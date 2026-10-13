//! Local NJU-seq runner: explicit argv, sample isolation, checked exit status.
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

pub fn json(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c < ' ' => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
fn p(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}
pub struct Config {
    pub samples: PathBuf,
    pub from: String,
    pub exons: PathBuf,
    pub out: PathBuf,
    pub bac_index: Option<PathBuf>,
    pub rrna_index: Option<PathBuf>,
    pub mrna_index: Option<PathBuf>,
    pub adapter_r1: String,
    pub adapter_r2: String,
    pub threads: usize,
    pub jobs: usize,
    pub exe: PathBuf,
    pub storage: String,
    pub rrna_filter: String,
    pub qc: bool,
}

impl Config {
    fn aligns_to(&self, mode: &str) -> bool {
        match mode {
            "mrna" => true,
            "bac" => matches!(self.rrna_filter.as_str(), "both" | "bac"),
            "rrna" => matches!(self.rrna_filter.as_str(), "both" | "host"),
            _ => false,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum UmiDesign {
    None,
    R1ThreePrime8,
}

impl UmiDesign {
    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "none" => Ok(Self::None),
            "r1_3p8" => Ok(Self::R1ThreePrime8),
            _ => Err(format!(
                "Invalid manifest umi value {value:?}; expected none or r1_3p8"
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::R1ThreePrime8 => "r1_3p8",
        }
    }

    fn enabled(self) -> bool {
        self != Self::None
    }

    fn count_modes(self) -> &'static [&'static str] {
        if self.enabled() {
            &["reads", "umi"]
        } else {
            &["reads"]
        }
    }
}

struct Sample {
    name: String,
    role: String,
    r1: PathBuf,
    r2: PathBuf,
    alignment: PathBuf,
    umi: UmiDesign,
}

fn umi_score_status(samples: &[Sample]) -> &'static str {
    let enabled = samples.iter().filter(|sample| sample.umi.enabled()).count();
    if enabled == 0 {
        "disabled_no_umi"
    } else if enabled == samples.len() {
        "enabled"
    } else {
        "disabled_mixed_libraries"
    }
}

fn score_modes(samples: &[Sample]) -> &'static [&'static str] {
    if umi_score_status(samples) == "enabled" {
        &["reads", "umi"]
    } else {
        &["reads"]
    }
}
fn absolute(path: &Path) -> Result<PathBuf, String> {
    if path.is_absolute() {
        Ok(path.to_owned())
    } else {
        Ok(std::env::current_dir()
            .map_err(|e| e.to_string())?
            .join(path))
    }
}
fn samples(config: &Config) -> Result<Vec<Sample>, String> {
    let path = absolute(&config.samples)?;
    let base = path.parent().unwrap();
    let text = fs::read_to_string(&path).map_err(|e| e.to_string())?;
    let mut rows = text.lines();
    let has_umi = match rows.next() {
        Some("sample\trole\tr1\tr2\talignment") => false,
        Some("sample\trole\tr1\tr2\talignment\tumi") => true,
        _ => return Err("Expected manifest header: sample\trole\tr1\tr2\talignment[\tumi]".into()),
    };
    let mut names = std::collections::HashSet::new();
    let mut result = vec![];
    for row in rows {
        if row.trim().is_empty() {
            continue;
        }
        let f: Vec<_> = row.split('\t').collect();
        if f.len() != if has_umi { 6 } else { 5 }
            || f[0].is_empty()
            || matches!(f[0], "." | "..")
            || !f[0]
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
            || !["NC", "TR"].contains(&f[1])
            || !names.insert(f[0].to_owned())
        {
            return Err(format!("Invalid/duplicate manifest row: {row}"));
        }
        let resolve = |s: &str| {
            if Path::new(s).is_absolute() {
                PathBuf::from(s)
            } else {
                base.join(s)
            }
        };
        let s = Sample {
            name: f[0].into(),
            role: f[1].into(),
            r1: resolve(f[2]),
            r2: resolve(f[3]),
            alignment: resolve(f[4]),
            // Legacy five-column manifests describe ordinary libraries. UMI
            // extraction/deduplication always requires explicit declaration.
            umi: if has_umi {
                UmiDesign::parse(f[5])?
            } else {
                UmiDesign::None
            },
        };
        let mut required = vec![&s.r1];
        if ["raw", "prepared"].contains(&config.from.as_str()) {
            required.push(&s.r2);
        } else {
            required.push(&s.alignment);
        }
        for f in required {
            if !f.is_file() {
                return Err(format!("Input missing: {}", f.display()));
            }
        }
        result.push(s);
    }
    if result.iter().filter(|s| s.role == "NC").count() != 1
        || !result.iter().any(|s| s.role == "TR")
    {
        return Err("Exactly one NC and at least one TR required".into());
    }
    Ok(result)
}
#[derive(Clone)]
struct Step {
    name: String,
    program: String,
    args: Vec<String>,
}
impl Step {
    fn json(&self) -> String {
        format!(
            "{{\"stage\":{},\"program\":{},\"argv\":[{}]}}",
            json(&self.name),
            json(&self.program),
            self.args
                .iter()
                .map(|s| json(s))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
    fn run(&self, logs: &Path, journal: &Path) -> Result<(), String> {
        let log_path = logs.join(format!("{}.log", self.name));
        let out = File::create(&log_path).map_err(|e| e.to_string())?;
        let start = Instant::now();
        let status = Command::new(&self.program)
            .args(&self.args)
            .stdout(Stdio::from(out.try_clone().map_err(|e| e.to_string())?))
            .stderr(Stdio::from(out))
            .status();
        let elapsed = start.elapsed().as_secs_f64();
        let code = status.as_ref().ok().and_then(|s| s.code()).unwrap_or(-1);
        let mut journal = OpenOptions::new()
            .create(true)
            .append(true)
            .open(journal)
            .map_err(|e| e.to_string())?;
        writeln!(
            journal,
            "{}\t{elapsed:.6}\t{code}\t{}",
            self.name,
            self.json()
        )
        .map_err(|e| e.to_string())?;
        match status {
            Ok(s) if s.success() => Ok(()),
            Ok(s) => Err(format!(
                "Stage {} failed ({s}); {}",
                self.name,
                log_path.display()
            )),
            Err(e) => Err(format!("Stage {} could not start: {e}", self.name)),
        }
    }
}
fn step(name: &str, program: &str, args: Vec<String>) -> Step {
    Step {
        name: name.into(),
        program: program.into(),
        args,
    }
}
fn strings(s: &[&str]) -> Vec<String> {
    s.iter().map(|s| (*s).into()).collect()
}
fn native(name: &str, c: &Config, args: Vec<String>) -> Step {
    step(name, &p(&c.exe), args)
}
fn legacy_sample_steps(c: &Config, s: &Sample, root: &Path) -> Vec<Step> {
    let dir = root.join(&s.name);
    let f = |n: &str| p(&dir.join(n));
    let mut steps = vec![];
    let mut r1 = p(&s.r1);
    let mut alignment = p(&s.alignment);
    if c.from == "raw" {
        let mut cutadapt = strings(&["-a", &c.adapter_r1, "-A", &c.adapter_r2]);
        if s.umi.enabled() {
            cutadapt.extend(strings(&["-U", "8"]));
        }
        cutadapt.extend(strings(&[
            "-O",
            "6",
            "-m",
            if s.umi.enabled() { "18:10" } else { "10:10" },
            "-e",
            "0.1",
            "--discard-untrimmed",
            "-j",
            &c.threads.to_string(),
            "-o",
            &f("R1.cut.fq.gz"),
            "-p",
            &f("R2.cut.fq.gz"),
            &p(&s.r1),
            &p(&s.r2),
        ]));
        steps.push(step("01_cutadapt", "cutadapt", cutadapt));
        if s.umi.enabled() {
            steps.push(step(
                "02_extract",
                "umi_tools",
                strings(&[
                    "extract",
                    "--extract-method=string",
                    "--bc-pattern=NNNNNNNN",
                    "--3prime",
                    "--stdin",
                    &f("R1.cut.fq.gz"),
                    "--stdout",
                    &f("R1.origin.fq.gz"),
                    "--read2-in",
                    &f("R2.cut.fq.gz"),
                    "--read2-out",
                    &f("R2.origin.fq.gz"),
                    "--log",
                    &f("umi_extract.log"),
                ]),
            ));
        }
        let mut input = if s.umi.enabled() { "origin" } else { "cut" };
        for (number, mode, index, sam, output) in [
            (3, "bac", c.bac_index.as_ref(), "bac_rna.raw.sam", "filter"),
            (
                6,
                "rrna",
                c.rrna_index.as_ref(),
                "rrna.filter.sam",
                "filter.mrna",
            ),
            (9, "mrna", c.mrna_index.as_ref(), "mrna.filter.sam", ""),
        ] {
            if !c.aligns_to(mode) {
                continue;
            }
            let mut a = strings(&[
                "--end-to-end",
                "-p",
                &c.threads.to_string(),
                "-t",
                "--no-unal",
                "--no-mixed",
                "--no-discordant",
                "-D",
                "20",
                "-R",
                "3",
                "-N",
                "0",
                "-L",
                "10",
                "-i",
                "S,1,0.75",
                "--xeq",
            ]);
            if mode == "mrna" {
                a.push("-a".into());
            } else {
                a.extend(strings(&["-k", "5"]));
            }
            if mode == "rrna" {
                a.extend(strings(&["--np", "0"]));
            } else {
                a.extend(strings(&[
                    "--maxins",
                    "120",
                    "--score-min",
                    if mode == "bac" {
                        "L,-0.8,-0.8"
                    } else {
                        "L,0.8,-0.8"
                    },
                    "--rdg",
                    "8,8",
                    "--rfg",
                    "8,8",
                    "--np",
                    "12",
                    "--mp",
                    "12,12",
                    "--ignore-quals",
                ]));
            }
            a.extend(strings(&[
                "-x",
                &p(index.unwrap()),
                "-1",
                &f(&format!("R1.{input}.fq.gz")),
                "-2",
                &f(&format!("R2.{input}.fq.gz")),
                "-S",
                &f(sam),
            ]));
            steps.push(step(&format!("{number:02}_{mode}_align"), "bowtie2", a));
            if mode != "mrna" {
                let selected = f(&format!("{mode}.selected.tsv"));
                steps.push(native(
                    &format!("{:02}_{mode}_select", number + 1),
                    c,
                    strings(&[
                        "filter",
                        "--mode",
                        mode,
                        "--format",
                        "sam",
                        "--input",
                        &f(sam),
                        "--output",
                        &selected,
                    ]),
                ));
                steps.push(native(
                    &format!("{:02}_{mode}_remove", number + 2),
                    c,
                    strings(&[
                        "delete-pairs",
                        "--names",
                        &selected,
                        "--r1",
                        &f(&format!("R1.{input}.fq.gz")),
                        "--r2",
                        &f(&format!("R2.{input}.fq.gz")),
                        "--out1",
                        &f(&format!("R1.{output}.fq.gz")),
                        "--out2",
                        &f(&format!("R2.{output}.fq.gz")),
                    ]),
                ));
                input = output;
            }
        }
        r1 = f(&format!("R1.{input}.fq.gz"));
        alignment = f("mrna.filter.sam");
    }
    steps.push(native(
        "10_select",
        c,
        strings(&[
            "filter",
            "--mode",
            "mrna",
            "--format",
            if c.from == "filtered" { "tsv" } else { "sam" },
            "--input",
            &alignment,
            "--exons",
            &p(&c.exons),
            "--output",
            &f("mrna.dedup.filter.tmp"),
        ]),
    ));
    steps.push(native(
        "11_almostunique",
        c,
        strings(&[
            "almostunique",
            "--input",
            &f("mrna.dedup.filter.tmp"),
            "--r1",
            &r1,
            "--output",
            &f("mrna.almostuniquematch.filter.tmp"),
        ]),
    ));
    if s.umi.enabled() {
        steps.push(native(
            "12_umi",
            c,
            strings(&[
                "umi",
                "--input",
                &f("mrna.almostuniquematch.filter.tmp"),
                "--exons",
                &p(&c.exons),
                "--sample",
                &s.name,
                "--output",
                &f("mrna.umi.filter.tmp"),
                "--families",
                &f("mrna.umi.families.tsv.gz"),
                "--summary",
                &f("mrna.umi.summary.json"),
            ]),
        ));
    }
    for &mode in s.umi.count_modes() {
        let input = if mode == "reads" {
            "mrna.almostuniquematch.filter.tmp"
        } else {
            "mrna.umi.filter.tmp"
        };
        steps.push(native(
            &format!("13_{mode}_count"),
            c,
            strings(&[
                "count",
                "--input",
                &f(input),
                "--output",
                &f(&format!("{mode}/mrna.count.filter.tmp")),
            ]),
        ));
        steps.push(native(
            &format!("14_{mode}_merge"),
            c,
            strings(&[
                "merge",
                "--input",
                &f(&format!("{mode}/mrna.count.filter.tmp")),
                "--exons",
                &p(&c.exons),
                "--output",
                &f(&format!("{mode}/mrna.filter.tsv")),
            ]),
        ));
    }
    steps
}
fn sample_steps(c: &Config, s: &Sample, root: &Path) -> Vec<Step> {
    let old = legacy_sample_steps(c, s, root);
    if c.storage == "all" && !c.qc && !["prepared", "unique"].contains(&c.from.as_str()) {
        return old;
    }
    let dir = root.join(&s.name);
    let f = |name: &str| p(&dir.join(name));
    let mut steps = vec![];
    let mut one = p(&s.r1);
    let mut two = p(&s.r2);
    let qc_pair = |stage: &str, one: &str, two: &str| {
        native(
            &format!("qc_{stage}"),
            c,
            strings(&[
                "qc",
                "--r1",
                one,
                "--r2",
                two,
                "--sample",
                &s.name,
                "--stage",
                stage,
                "--out",
                &f(&format!("qc_{stage}")),
            ]),
        )
    };
    let mut unique = f("mrna.almostuniquematch.filter.tmp");
    if ["raw", "prepared"].contains(&c.from.as_str()) {
        if c.from == "raw" {
            if c.qc {
                steps.push(qc_pair("raw", &one, &two));
            }
            for step in old
                .iter()
                .filter(|step| step.name.starts_with("01_") || step.name.starts_with("02_"))
            {
                steps.push(step.clone());
                if c.qc && step.name.starts_with("01_") {
                    steps.push(qc_pair("cut", &f("R1.cut.fq.gz"), &f("R2.cut.fq.gz")));
                }
            }
            let suffix = if s.umi.enabled() { "origin" } else { "cut" };
            one = f(&format!("R1.{suffix}.fq.gz"));
            two = f(&format!("R2.{suffix}.fq.gz"));
        }
        if c.qc {
            steps.push(qc_pair("prepared", &one, &two));
        }
        for (mode, number, index, suffix) in [
            ("bac", 3, c.bac_index.as_ref(), "filter"),
            ("rrna", 6, c.rrna_index.as_ref(), "filter.mrna"),
            ("mrna", 9, c.mrna_index.as_ref(), ""),
        ] {
            if !c.aligns_to(mode) {
                continue;
            }
            let output = if mode == "mrna" {
                unique.clone()
            } else {
                f(&format!("{mode}.removed.ids"))
            };
            let mut args = strings(&[
                "align-select",
                "--mode",
                mode,
                "--index",
                &p(index.unwrap()),
                "--r1",
                &one,
                "--r2",
                &two,
                "--threads",
                &c.threads.to_string(),
                "--output",
                &output,
                "--log",
                &f(&format!("{mode}.bowtie2.log")),
                "--stats",
                &f(&format!("{mode}.selection.json")),
            ]);
            if mode == "mrna" {
                args.extend(strings(&["--exons", &p(&c.exons)]));
            }
            if c.storage == "all" {
                args.extend(strings(&["--sam-output", &f(&format!("{mode}.sam"))]));
            }
            steps.push(native(&format!("{number:02}_{mode}_stream"), c, args));
            if mode != "mrna" {
                let next_one = f(&format!("R1.{suffix}.fq.gz"));
                let next_two = f(&format!("R2.{suffix}.fq.gz"));
                steps.push(native(
                    &format!("{:02}_{mode}_remove", number + 2),
                    c,
                    strings(&[
                        "delete-pairs",
                        "--names",
                        &output,
                        "--r1",
                        &one,
                        "--r2",
                        &two,
                        "--out1",
                        &next_one,
                        "--out2",
                        &next_two,
                    ]),
                ));
                one = next_one;
                two = next_two;
                if c.qc {
                    steps.push(qc_pair(
                        if mode == "bac" {
                            "after_bac"
                        } else {
                            "after_host"
                        },
                        &one,
                        &two,
                    ));
                }
            }
        }
    } else if c.from == "unique" {
        unique = p(&s.alignment);
    } else {
        steps.extend(
            old.iter()
                .filter(|step| step.name.starts_with("10_") || step.name.starts_with("11_"))
                .cloned(),
        );
    }
    if c.qc {
        steps.push(native(
            "qc_mrna",
            c,
            strings(&[
                "qc",
                "--input",
                &unique,
                "--sample",
                &s.name,
                "--stage",
                "mrna_unique",
                "--out",
                &f("qc_mrna"),
            ]),
        ));
    }
    for mut step in old.into_iter().filter(|step| step.name.as_str() >= "12_") {
        for arg in &mut step.args {
            if *arg == f("mrna.almostuniquematch.filter.tmp") {
                *arg = unique.clone();
            }
        }
        let is_umi = step.name == "12_umi";
        steps.push(step);
        if is_umi && c.qc {
            steps.push(native(
                "qc_umi",
                c,
                strings(&[
                    "qc",
                    "--input",
                    &f("mrna.umi.filter.tmp"),
                    "--sample",
                    &s.name,
                    "--stage",
                    "umi_molecules",
                    "--unit",
                    "molecules",
                    "--out",
                    &f("qc_umi"),
                ]),
            ));
        }
    }
    steps
}
fn score_steps(c: &Config, samples: &[Sample], root: &Path) -> Vec<Step> {
    score_modes(samples)
        .iter()
        .map(|mode| {
            let mut args = strings(&["score", "--mode", "mrna", "--format", "unified"]);
            for role in ["NC", "TR"] {
                for s in samples.iter().filter(|s| s.role == role) {
                    args.push(format!("--{}", role.to_lowercase()));
                    args.push(p(&root
                        .join("score_inputs")
                        .join(mode)
                        .join(&s.name)
                        .join("mrna.filter.tsv")));
                }
            }
            args.extend(strings(&[
                "--output",
                &p(&root.join(format!("mrna.{mode}.scored.tsv"))),
                "--summary",
                &p(&root.join(format!("mrna.{mode}.score_summary.json"))),
            ]));
            native(&format!("15_{mode}_score"), c, args)
        })
        .collect()
}
pub fn execute(c: Config, dry_run: bool) -> Result<(), String> {
    for adapter in [&c.adapter_r1, &c.adapter_r2] {
        if adapter.is_empty() || !adapter.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
            return Err("Adapters must be nonempty uppercase IUPAC DNA sequences".into());
        }
    }
    if !["raw", "prepared", "sam", "filtered", "unique"].contains(&c.from.as_str())
        || c.threads == 0
        || c.jobs == 0
    {
        return Err("from=raw|prepared|sam|filtered|unique; threads/jobs must be positive".into());
    }
    if !["compact", "all"].contains(&c.storage.as_str())
        || !["both", "bac", "host", "none"].contains(&c.rrna_filter.as_str())
    {
        return Err("storage=compact|all; rrna-filter=both|bac|host|none".into());
    }
    if !c.exons.is_file() {
        return Err("Missing exon annotation".into());
    }
    if ["raw", "prepared"].contains(&c.from.as_str()) {
        for (mode, prefix) in [
            ("bac", &c.bac_index),
            ("rrna", &c.rrna_index),
            ("mrna", &c.mrna_index),
        ] {
            if !c.aligns_to(mode) {
                continue;
            }
            let prefix = prefix.as_ref().ok_or_else(|| {
                format!("{} entry requires --{mode}-index or --index DIR", c.from)
            })?;
            let complete = |extension: &str| {
                ["1", "2", "3", "4", "rev.1", "rev.2"].iter().all(|part| {
                    PathBuf::from(format!("{}.{part}.{extension}", p(prefix))).is_file()
                })
            };
            if !complete("bt2") && !complete("bt2l") {
                return Err(format!(
                    "Index requires six components of one format (.bt2 or .bt2l): {}",
                    prefix.display()
                ));
            }
        }
    }
    let samples = samples(&c)?;
    let root = absolute(&c.out)?;
    let groups: Vec<_> = samples.iter().map(|s| sample_steps(&c, s, &root)).collect();
    let scores = score_steps(&c, &samples, &root);
    let mut plan = format!(
        "{{\"version\":\"0.5.0\",\"from\":{},\"threads_per_alignment\":{},\"sample_jobs\":{},\"out\":{},\"umi_score_status\":{},\"samples\":[",
        json(&c.from),
        c.threads,
        c.jobs,
        json(&p(&root)),
        json(umi_score_status(&samples))
    );
    plan.push_str(
        &samples
            .iter()
            .zip(&groups)
            .map(|(s, steps)| {
                format!(
                    "{{\"sample\":{},\"role\":{},\"umi\":{},\"steps\":[{}]}}",
                    json(&s.name),
                    json(&s.role),
                    json(s.umi.as_str()),
                    steps.iter().map(Step::json).collect::<Vec<_>>().join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(","),
    );
    plan.push_str(&format!(
        "],\"storage\":{},\"rrna_filter\":{},\"qc\":{},\"score_steps\":[{}]}}\n",
        json(&c.storage),
        json(&c.rrna_filter),
        c.qc,
        scores.iter().map(Step::json).collect::<Vec<_>>().join(",")
    ));
    if dry_run {
        print!("{plan}");
        return Ok(());
    }
    fs::create_dir(&root).map_err(|e| {
        format!(
            "Refusing existing/uncreatable output {}: {e}",
            root.display()
        )
    })?;
    fs::write(root.join("plan.json"), plan).map_err(|e| e.to_string())?;
    for s in &samples {
        for sub in std::iter::once("logs").chain(s.umi.count_modes().iter().copied()) {
            fs::create_dir_all(root.join(&s.name).join(sub)).map_err(|e| e.to_string())?;
        }
    }
    let index = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let errors = Mutex::new(vec![]);
    std::thread::scope(|scope| {
        for _ in 0..c.jobs.min(samples.len()) {
            scope.spawn(|| {
                loop {
                    if failed.load(Ordering::SeqCst) {
                        break;
                    }
                    let i = index.fetch_add(1, Ordering::SeqCst);
                    if i >= samples.len() {
                        break;
                    }
                    let dir = root.join(&samples[i].name);
                    eprintln!("START {}", samples[i].name);
                    for step in &groups[i] {
                        if failed.load(Ordering::SeqCst) {
                            break;
                        }
                        if let Err(e) = step.run(&dir.join("logs"), &dir.join("stages.tsv")) {
                            failed.store(true, Ordering::SeqCst);
                            errors.lock().unwrap().push(e);
                            break;
                        }
                    }
                    if !failed.load(Ordering::SeqCst) {
                        eprintln!("DONE {}", samples[i].name);
                    }
                }
            });
        }
    });
    if failed.load(Ordering::SeqCst) {
        let e = errors.lock().unwrap().join("\n");
        fs::write(root.join("FAILED"), &e).map_err(|e| e.to_string())?;
        return Err(e);
    }
    for mode in score_modes(&samples) {
        for s in &samples {
            let dir = root.join("score_inputs").join(mode).join(&s.name);
            fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
            std::os::unix::fs::symlink(
                root.join(&s.name).join(mode).join("mrna.filter.tsv"),
                dir.join("mrna.filter.tsv"),
            )
            .map_err(|e| e.to_string())?;
        }
    }
    fs::create_dir(root.join("logs")).map_err(|e| e.to_string())?;
    for step in scores {
        if let Err(e) = step.run(&root.join("logs"), &root.join("stages.tsv")) {
            fs::write(root.join("FAILED"), &e).map_err(|e| e.to_string())?;
            return Err(e);
        }
    }
    if c.qc
        && let Err(e) = crate::report::run(&root, &root.join("qc.html"))
    {
        fs::write(root.join("FAILED"), &e).map_err(|e| e.to_string())?;
        return Err(e);
    }
    // Only names created by this runner inside its new output directory are
    // eligible for cleanup. Input FASTQ and caller-owned alignments are untouched.
    if c.storage == "compact" {
        let mut audit = String::from("path\tbytes\n");
        for s in &samples {
            for name in [
                "R1.cut.fq.gz",
                "R2.cut.fq.gz",
                "R1.origin.fq.gz",
                "R2.origin.fq.gz",
                "R1.filter.fq.gz",
                "R2.filter.fq.gz",
                "R1.filter.mrna.fq.gz",
                "R2.filter.mrna.fq.gz",
                "mrna.dedup.filter.tmp",
                "reads/mrna.count.filter.tmp",
                "umi/mrna.count.filter.tmp",
            ] {
                let path = root.join(&s.name).join(name);
                if let Ok(meta) = fs::symlink_metadata(&path) {
                    if !meta.file_type().is_file() {
                        return Err("Unexpected cleanup file type".into());
                    }
                    audit.push_str(&format!("{}/{}\t{}\n", s.name, name, meta.len()));
                    fs::remove_file(path).map_err(|e| e.to_string())?;
                }
            }
        }
        fs::write(root.join("cleaned_intermediates.tsv"), audit).map_err(|e| e.to_string())?;
    }
    fs::write(
        root.join("SUCCESS"),
        format!(
            "All sample stages and {} joint score branch(es) completed. UMI score status: {}.\n",
            score_modes(&samples).join(" + "),
            umi_score_status(&samples)
        ),
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicU64;

    fn config(from: &str) -> Config {
        Config {
            samples: "/fixture/samples.tsv".into(),
            from: from.into(),
            exons: "/fixture/exon.info".into(),
            out: "/fixture/output".into(),
            bac_index: Some("/fixture/bac".into()),
            rrna_index: Some("/fixture/rrna".into()),
            mrna_index: Some("/fixture/mrna".into()),
            adapter_r1: "AGATCGGAAGAGCACA".into(),
            adapter_r2: "GATCGTCGGACTGTAG".into(),
            threads: 1,
            jobs: 1,
            exe: "/fixture/nju-seq".into(),
            storage: "all".into(),
            rrna_filter: "both".into(),
            qc: false,
        }
    }

    fn sample(name: &str, role: &str, umi: UmiDesign) -> Sample {
        Sample {
            name: name.into(),
            role: role.into(),
            r1: format!("/fixture/{name}.R1.fq.gz").into(),
            r2: format!("/fixture/{name}.R2.fq.gz").into(),
            alignment: format!("/fixture/{name}.sam").into(),
            umi,
        }
    }

    fn argument<'a>(step: &'a Step, option: &str) -> &'a str {
        &step.args[step.args.iter().position(|arg| arg == option).unwrap() + 1]
    }

    #[test]
    fn ordinary_raw_libraries_do_not_trim_or_extract_umi_bases() {
        let c = config("raw");
        let s = sample("NC", "NC", UmiDesign::None);
        let steps = sample_steps(&c, &s, &c.out);
        let cut = &steps[0];
        assert_eq!(cut.name, "01_cutadapt");
        assert_eq!(argument(cut, "-m"), "10:10");
        assert!(!cut.args.iter().any(|arg| arg == "-U"));
        assert!(
            !steps
                .iter()
                .any(|step| step.name == "02_extract" || step.name == "12_umi")
        );
        assert!(
            !steps
                .iter()
                .any(|step| step.name == "13_umi_count" || step.name == "14_umi_merge")
        );
        let bac = steps
            .iter()
            .find(|step| step.name == "03_bac_align")
            .unwrap();
        assert_eq!(argument(bac, "-1"), "/fixture/output/NC/R1.cut.fq.gz");
        assert_eq!(argument(bac, "-2"), "/fixture/output/NC/R2.cut.fq.gz");
        let removal = steps
            .iter()
            .find(|step| step.name == "05_bac_remove")
            .unwrap();
        assert_eq!(argument(removal, "--r1"), "/fixture/output/NC/R1.cut.fq.gz");
        assert!(steps.iter().any(|step| step.name == "13_reads_count"));
        assert_eq!(s.umi.count_modes(), &["reads"]);
    }

    #[test]
    fn explicit_r1_3p8_keeps_current_trimming_and_both_count_branches() {
        let c = config("raw");
        let s = sample("TR", "TR", UmiDesign::R1ThreePrime8);
        let steps = sample_steps(&c, &s, &c.out);
        assert_eq!(argument(&steps[0], "-U"), "8");
        assert_eq!(argument(&steps[0], "-m"), "18:10");
        let extract = steps.iter().find(|step| step.name == "02_extract").unwrap();
        assert!(extract.args.iter().any(|arg| arg == "--3prime"));
        assert!(
            extract
                .args
                .iter()
                .any(|arg| arg == "--bc-pattern=NNNNNNNN")
        );
        let bac = steps
            .iter()
            .find(|step| step.name == "03_bac_align")
            .unwrap();
        assert_eq!(argument(bac, "-1"), "/fixture/output/TR/R1.origin.fq.gz");
        for name in [
            "12_umi",
            "13_reads_count",
            "14_reads_merge",
            "13_umi_count",
            "14_umi_merge",
        ] {
            assert!(steps.iter().any(|step| step.name == name), "missing {name}");
        }
        assert_eq!(s.umi.count_modes(), &["reads", "umi"]);
    }

    #[test]
    fn non_raw_optional_umi_does_not_change_read_filtering() {
        for from in ["sam", "filtered"] {
            let c = config(from);
            let ordinary = sample_steps(&c, &sample("NC", "NC", UmiDesign::None), &c.out);
            let tagged = sample_steps(&c, &sample("NC", "NC", UmiDesign::R1ThreePrime8), &c.out);
            assert_eq!(ordinary[0].name, "10_select");
            assert_eq!(ordinary[0].args, tagged[0].args);
            assert_eq!(ordinary[1].args, tagged[1].args);
            assert!(!ordinary.iter().any(|step| step.name == "12_umi"));
            assert!(tagged.iter().any(|step| step.name == "12_umi"));
        }
    }

    #[test]
    fn joint_umi_requires_every_library_without_dropping_treatments() {
        let c = config("filtered");
        for designs in [
            [UmiDesign::None, UmiDesign::None, UmiDesign::None],
            [
                UmiDesign::R1ThreePrime8,
                UmiDesign::None,
                UmiDesign::R1ThreePrime8,
            ],
            [
                UmiDesign::R1ThreePrime8,
                UmiDesign::R1ThreePrime8,
                UmiDesign::R1ThreePrime8,
            ],
        ] {
            let samples = vec![
                sample("NC", "NC", designs[0]),
                sample("TR1", "TR", designs[1]),
                sample("TR2", "TR", designs[2]),
            ];
            let steps = score_steps(&c, &samples, &c.out);
            assert_eq!(steps[0].name, "15_reads_score");
            assert_eq!(steps[0].args.iter().filter(|arg| *arg == "--nc").count(), 1);
            assert_eq!(steps[0].args.iter().filter(|arg| *arg == "--tr").count(), 2);
            for name in ["NC", "TR1", "TR2"] {
                assert!(steps[0].args.contains(&format!(
                    "/fixture/output/score_inputs/reads/{name}/mrna.filter.tsv"
                )));
            }
            let enabled = designs.iter().filter(|design| design.enabled()).count();
            let status = match enabled {
                0 => "disabled_no_umi",
                3 => "enabled",
                _ => "disabled_mixed_libraries",
            };
            assert_eq!(umi_score_status(&samples), status);
            assert_eq!(steps.len(), if enabled == 3 { 2 } else { 1 });
            if enabled == 3 {
                assert_eq!(steps[1].name, "15_umi_score");
            }
        }
    }

    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);
    struct ManifestFixture(PathBuf);
    impl ManifestFixture {
        fn new() -> Self {
            let root = std::env::temp_dir().join(format!(
                "nju-umi-manifest-{}-{}",
                std::process::id(),
                NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&root).unwrap();
            for name in ["R1.fq", "R2.fq"] {
                fs::write(root.join(name), "").unwrap();
            }
            Self(root)
        }
        fn parse(&self, text: &str) -> Result<Vec<Sample>, String> {
            let mut c = config("raw");
            c.samples = self.0.join("samples.tsv");
            fs::write(&c.samples, text).unwrap();
            samples(&c)
        }
    }
    impl Drop for ManifestFixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn five_column_manifest_defaults_to_no_umi() {
        let fixture = ManifestFixture::new();
        let parsed = fixture.parse("sample\trole\tr1\tr2\talignment\nNC\tNC\tR1.fq\tR2.fq\t-\nTR\tTR\tR1.fq\tR2.fq\t-\n").unwrap();
        assert!(parsed.iter().all(|sample| sample.umi == UmiDesign::None));
        assert_eq!(umi_score_status(&parsed), "disabled_no_umi");
    }

    #[test]
    fn sixth_column_declares_design_and_rejects_unknown_or_missing_values() {
        let fixture = ManifestFixture::new();
        let text = "sample\trole\tr1\tr2\talignment\tumi\nNC\tNC\tR1.fq\tR2.fq\t-\tnone\nTR\tTR\tR1.fq\tR2.fq\t-\tr1_3p8\n";
        let parsed = fixture.parse(text).unwrap();
        assert_eq!(parsed[0].umi.as_str(), "none");
        assert_eq!(parsed[1].umi.as_str(), "r1_3p8");
        assert_eq!(umi_score_status(&parsed), "disabled_mixed_libraries");
        for invalid in [
            text.replace("r1_3p8", "yes"),
            text.replace("r1_3p8", ""),
            text.replace("\t-\tr1_3p8", "\t-"),
        ] {
            assert!(fixture.parse(&invalid).is_err());
        }
    }
}
