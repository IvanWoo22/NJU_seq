use nju_seq::{
    ablation, align, counts, fastq, filter, model::ExonIndex, proof, qc, report, scoring, umi,
    workflow,
};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

const HELP: &str = "nju-seq 0.5.0 — NJU-seq Rust workflow, human and Arabidopsis\n\
Commands:\n\
  prove-impossible-short --input PATH --format fastq|tsv --reference INDEX_TRANSCRIPT_FASTA --output IDS --stats JSON\n\
  length-select --input ACCEPTED_TSV --min-length NT --output TSV (explicit sensitivity test)\n\
  align-select --mode bac|rrna|mrna --index PREFIX --r1 FQ --r2 FQ [--exons PATH] [--threads 4] --output TSV --log LOG --stats JSON [--sam-output SAM]\n\
      optional experimental --screen-multimap K [--screen-kind paired|exact-pair] emits rejection IDs; unresolved reads still require -a.\n\
  qc --r1 FQ [--r2 FQ] --sample NAME --stage NAME --out DIR\n\
  qc --input ALIGNMENT_TSV --sample NAME --stage NAME [--unit reads|molecules] --out DIR\n\
  qc-report --root DIR --output HTML\n\
  partition-pairs --r1 FQ --r2 FQ --bac-list IDS --host-list IDS --sample NAME --out DIR\n\
  ablation-assemble --baseline TSV --bac TSV --host TSV --bac-host-list IDS --out DIR\n\
  filter --input PATH --format sam|tsv --mode mrna|bac|rrna [--exons PATH] --output PATH\n\
  almostunique --input TSV --r1 FASTQ[.gz] --output TSV\n\
  delete-pairs --names ID_OR_ALIGNMENT_TSV --r1 FQ --r2 FQ --out1 FQ[.gz] --out2 FQ[.gz]\n\
  umi --input TSV --exons PATH --sample NAME --output TSV --families TSV[.gz] --summary JSON [--umi-length 8]\n\
  count --input TSV --output TSV\n\
  merge --input TSV --exons PATH --output TSV\n\
  score --mode mrna|rrna --nc TSV --tr TSV [--tr TSV ...] --output TSV [--summary JSON]\n\
      [--format unified|legacy] [--rows candidates|all] [--reference-id ID (required for rrna)]\n\
      defaults: mode=mrna, format=unified; rows=candidates for mrna, all for rrna.\n\
  plan|run --samples TSV --from raw|prepared|sam|filtered|unique --exons PATH --out DIR [--threads 4] [--jobs 1]\n\
      [--storage compact|all] [--rrna-filter both|bac|host|none] [--qc yes|no]\n\
      [--adapter-r1 DNA] [--adapter-r2 DNA]; --index DIR uses bac / rrna / mrna prefixes\n\
      raw/prepared require an mRNA index and only the enabled rRNA indices (--index DIR or individual prefixes).\n\
      Individual prefixes: --bac-index PREFIX, --rrna-index PREFIX, --mrna-index PREFIX.\n\
      rrna-filter: both (default), host (skip bacterial alignment/removal), bac (skip host), none (skip both).\n\
      compact (default): streamed SAM selection, automatic QC, remove generated scratch FASTQ only after success.\n\
      prepared: adapters and optional UMI already removed; unique: accepted one-alignment-per-read five-column TSV.\n\
Manifest columns: sample role r1 r2 alignment [umi] (tabs; one NC, >=1 TR).\n\
umi: none (default without column) or r1_3p8 (R1 3prime 8nt; R2 5prime same tag).\n\
Reads scoring is always enabled; UMI joint scoring requires every library to have UMI.\n\
Relative manifest input paths are relative to its directory. Output must be new.\n\
SAM/5-column inputs must keep each read's alignments contiguous. Optional exact UMI only.\n\
Exon tables: Arabidopsis Parent=transcript:AT... or human Parent=ENST... with gene_id=ENSG...; PAR_Y IDs are preserved.\n\
Core run does not resume; scripts/align_chunked.py supports verified mRNA batch checkpoints.\n";

struct Args(HashMap<String, Vec<String>>);
impl Args {
    fn parse() -> Result<(String, Self), String> {
        let mut raw = std::env::args().skip(1);
        let verb = raw.next().unwrap_or_else(|| "help".into());
        let mut map = HashMap::<String, Vec<String>>::new();
        while let Some(key) = raw.next() {
            if !key.starts_with("--") {
                return Err(format!("Expected --option, got {key}"));
            }
            let value = raw
                .next()
                .ok_or_else(|| format!("Missing value for {key}"))?;
            if value.starts_with("--") {
                return Err(format!("Missing value for {key}"));
            }
            map.entry(key).or_default().push(value);
        }
        Ok((verb, Self(map)))
    }
    fn check(&self, allowed: &[&str]) -> Result<(), String> {
        for (k, v) in &self.0 {
            if !allowed.contains(&k.as_str()) || (v.len() > 1 && k != "--tr") {
                return Err(format!("Unknown or repeated option: {k}"));
            }
        }
        Ok(())
    }
    fn get(&self, k: &str) -> Result<&str, String> {
        self.0
            .get(k)
            .map(|v| v[0].as_str())
            .ok_or_else(|| format!("Missing {k}"))
    }
    fn optional(&self, k: &str) -> Option<&str> {
        self.0.get(k).map(|v| v[0].as_str())
    }
    fn path(&self, k: &str) -> Result<PathBuf, String> {
        Ok(PathBuf::from(self.get(k)?))
    }
    fn number(&self, k: &str, default: usize) -> Result<usize, String> {
        match self.optional(k) {
            Some(s) => s.parse().map_err(|_| format!("Invalid integer {k}: {s}")),
            None => Ok(default),
        }
    }
}

/// Publish complete outputs using no-overwrite hard links; failed stages leave no
/// public output. This is local-filesystem evaluation, not distributed locking.
struct Outputs {
    targets: Vec<PathBuf>,
    temp: Vec<PathBuf>,
}
impl Outputs {
    fn new(targets: Vec<PathBuf>) -> Result<Self, String> {
        let mut keys = HashSet::new();
        let mut temp = vec![];
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|e| e.to_string())?
            .as_nanos();
        for (i, p) in targets.iter().enumerate() {
            if fs::symlink_metadata(p).is_ok() {
                return Err(format!("Refusing to overwrite {}", p.display()));
            }
            let parent = p
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let name = p.file_name().ok_or("Output filename required")?;
            let key = parent
                .canonicalize()
                .map_err(|e| format!("Output parent unavailable: {e}"))?
                .join(name);
            if !keys.insert(key) {
                return Err("Output paths must be different".into());
            }
            temp.push(parent.join(format!(
                ".nju_eval.{}.{nonce}.{i}.partial",
                std::process::id()
            )));
        }
        Ok(Self { targets, temp })
    }
    fn commit(&self) -> Result<(), String> {
        let mut published = vec![];
        for (source, target) in self.temp.iter().zip(&self.targets) {
            if let Err(e) = fs::hard_link(source, target) {
                for p in published {
                    let _ = fs::remove_file(p);
                }
                return Err(format!("Publish {} failed: {e}", target.display()));
            }
            published.push(target);
        }
        Ok(())
    }
}
impl Drop for Outputs {
    fn drop(&mut self) {
        for p in &self.temp {
            let _ = fs::remove_file(p);
        }
    }
}

fn main_impl() -> Result<(), String> {
    let (verb, a) = Args::parse()?;
    match verb.as_str() {
        "help" | "--help" | "-h" => {
            print!("{HELP}");
            Ok(())
        }
        "--version" => {
            println!("nju-seq 0.5.0");
            Ok(())
        }
        "filter" => {
            a.check(&["--input", "--format", "--mode", "--exons", "--output"])?;
            let format = match a.get("--format")? {
                "sam" => filter::InputFormat::Sam,
                "tsv" => filter::InputFormat::AlignmentTsv,
                _ => return Err("format must be sam or tsv".into()),
            };
            let mode = match a.get("--mode")? {
                "mrna" => filter::Strategy::Mrna,
                "bac" => filter::Strategy::Bac,
                "rrna" => filter::Strategy::Rrna,
                _ => return Err("mode must be mrna, bac or rrna".into()),
            };
            let index = a
                .optional("--exons")
                .map(|p| ExonIndex::load(Path::new(p)))
                .transpose()?;
            let out = Outputs::new(vec![a.path("--output")?])?;
            let mut reader = fastq::Input::open(&a.path("--input")?)?;
            let mut writer = BufWriter::new(File::create(&out.temp[0]).map_err(|e| e.to_string())?);
            let stats = filter::run(&mut reader, &mut writer, format, mode, index.as_ref())?;
            reader.finish()?;
            writer.flush().map_err(|e| e.to_string())?;
            drop(writer);
            out.commit()?;
            eprintln!("{stats:?}");
            Ok(())
        }
        "align-select" => {
            a.check(&[
                "--mode",
                "--index",
                "--r1",
                "--r2",
                "--exons",
                "--threads",
                "--output",
                "--log",
                "--stats",
                "--sam-output",
                "--screen-multimap",
                "--screen-kind",
            ])?;
            let mut targets = vec![a.path("--output")?, a.path("--log")?, a.path("--stats")?];
            if let Some(s) = a.optional("--sam-output") {
                targets.push(PathBuf::from(s));
            }
            let out = Outputs::new(targets)?;
            let result = align::run(align::Options {
                mode: a.get("--mode")?,
                index: &a.path("--index")?,
                r1: &a.path("--r1")?,
                r2: &a.path("--r2")?,
                exons: a.optional("--exons").map(Path::new),
                threads: a.number("--threads", 4)?,
                output: &out.temp[0],
                log: &out.temp[1],
                sam: out.temp.get(3).map(PathBuf::as_path),
                screen_multimap: a
                    .optional("--screen-multimap")
                    .map(|s| {
                        s.parse()
                            .map_err(|_| "Invalid multimap screen limit".to_string())
                    })
                    .transpose()?,
                exact_pair_screen: match a.optional("--screen-kind").unwrap_or("paired") {
                    "paired" => false,
                    "exact-pair" => true,
                    _ => return Err("screen-kind=paired|exact-pair".into()),
                },
            });
            match result {
                Ok(stats) => {
                    fs::write(&out.temp[2], stats).map_err(|e| e.to_string())?;
                    out.commit()
                }
                Err(e) => {
                    if let Ok(log) = fs::read_to_string(&out.temp[1]) {
                        eprintln!("{log}");
                    }
                    Err(e)
                }
            }
        }
        "prove-impossible-short" => {
            a.check(&["--input", "--format", "--reference", "--output", "--stats"])?;
            let tsv = match a.get("--format")? {
                "tsv" => true,
                "fastq" => false,
                _ => return Err("format=fastq|tsv".into()),
            };
            let out = Outputs::new(vec![a.path("--output")?, a.path("--stats")?])?;
            let stats = proof::impossible_short(
                &a.path("--input")?,
                tsv,
                &a.path("--reference")?,
                &out.temp[0],
            )?;
            fs::write(&out.temp[1], stats).map_err(|e| e.to_string())?;
            out.commit()
        }
        "length-select" => {
            a.check(&["--input", "--min-length", "--output"])?;
            let out = Outputs::new(vec![a.path("--output")?])?;
            let stats = qc::length_select(
                &a.path("--input")?,
                &out.temp[0],
                a.number("--min-length", 0)?,
            )?;
            out.commit()?;
            eprintln!("input={} retained={}", stats.0, stats.1);
            Ok(())
        }
        "ablation-assemble" => {
            a.check(&["--baseline", "--bac", "--host", "--bac-host-list", "--out"])?;
            ablation::assemble(
                &a.path("--baseline")?,
                &a.path("--bac")?,
                &a.path("--host")?,
                &a.path("--bac-host-list")?,
                &a.path("--out")?,
            )
        }
        "qc-report" => {
            a.check(&["--root", "--output"])?;
            report::run(&a.path("--root")?, &a.path("--output")?)
        }
        "qc" => {
            a.check(&[
                "--r1", "--r2", "--input", "--sample", "--stage", "--unit", "--out",
            ])?;
            if let Some(input) = a.optional("--input") {
                qc::alignment_qc(
                    Path::new(input),
                    &a.path("--out")?,
                    a.get("--sample")?,
                    a.get("--stage")?,
                    a.optional("--unit").unwrap_or("reads"),
                )
            } else {
                qc::fastq_qc(
                    &a.path("--r1")?,
                    a.optional("--r2").map(Path::new),
                    &a.path("--out")?,
                    a.get("--sample")?,
                    a.get("--stage")?,
                )
            }
        }
        "partition-pairs" => {
            a.check(&[
                "--r1",
                "--r2",
                "--bac-list",
                "--host-list",
                "--sample",
                "--out",
            ])?;
            qc::partition(
                &a.path("--r1")?,
                &a.path("--r2")?,
                &a.path("--bac-list")?,
                &a.path("--host-list")?,
                &a.path("--out")?,
                a.get("--sample")?,
            )
        }
        "almostunique" => {
            a.check(&["--input", "--r1", "--output"])?;
            let out = Outputs::new(vec![a.path("--output")?])?;
            let (n, k) = fastq::almost_unique(&a.path("--input")?, &a.path("--r1")?, &out.temp[0])?;
            out.commit()?;
            eprintln!("input_alignments={n} retained_reads={k}");
            Ok(())
        }
        "delete-pairs" => {
            a.check(&["--names", "--r1", "--r2", "--out1", "--out2"])?;
            let targets = vec![a.path("--out1")?, a.path("--out2")?];
            let gzip = targets[0].extension().is_some_and(|e| e == "gz");
            if targets[1].extension().is_some_and(|e| e == "gz") != gzip {
                return Err("Both output mates must use the same compression".into());
            }
            let out = Outputs::new(targets)?;
            let (n, k) = fastq::delete_pairs(
                &a.path("--names")?,
                &a.path("--r1")?,
                &a.path("--r2")?,
                &out.temp[0],
                &out.temp[1],
                gzip,
            )?;
            out.commit()?;
            eprintln!("input_pairs={n} retained_pairs={k}");
            Ok(())
        }
        "umi" => {
            a.check(&[
                "--input",
                "--exons",
                "--sample",
                "--output",
                "--families",
                "--summary",
                "--umi-length",
            ])?;
            let index = ExonIndex::load(&a.path("--exons")?)?;
            let out = Outputs::new(vec![
                a.path("--output")?,
                a.path("--families")?,
                a.path("--summary")?,
            ])?;
            let mut reader = fastq::Input::open(&a.path("--input")?)?;
            let mut writer = BufWriter::new(File::create(&out.temp[0]).map_err(|e| e.to_string())?);
            let compressed = out.targets[1].extension().is_some_and(|e| e == "gz");
            let mut compressor = None;
            let mut family_writer: Box<dyn Write> = if compressed {
                let file = File::create(&out.temp[1]).map_err(|e| e.to_string())?;
                let mut child = Command::new("gzip")
                    .args(["-c", "-3"])
                    .stdin(Stdio::piped())
                    .stdout(Stdio::from(file))
                    .spawn()
                    .map_err(|e| e.to_string())?;
                let writer = BufWriter::new(child.stdin.take().unwrap());
                compressor = Some(child);
                Box::new(writer)
            } else {
                Box::new(BufWriter::new(
                    File::create(&out.temp[1]).map_err(|e| e.to_string())?,
                ))
            };
            let result = umi::run(
                &mut reader,
                &index,
                &mut writer,
                &mut family_writer,
                a.get("--sample")?,
                a.number("--umi-length", 8)?,
            );
            let flushed = family_writer.flush().map_err(|e| e.to_string());
            drop(family_writer);
            let compressor_status = compressor
                .as_mut()
                .map(|p| p.wait().map_err(|e| e.to_string()))
                .transpose()?;
            let summary = result?;
            flushed?;
            if compressor_status.is_some_and(|s| !s.success()) {
                return Err("gzip families output failed".into());
            }
            reader.finish()?;
            writer.flush().map_err(|e| e.to_string())?;
            drop(writer);
            fs::write(&out.temp[2], summary.to_json() + "\n").map_err(|e| e.to_string())?;
            out.commit()?;
            eprintln!(
                "input_reads={} output_molecules={}",
                summary.input_reads, summary.output_molecules
            );
            Ok(())
        }
        "count" | "merge" => {
            a.check(if verb == "count" {
                &["--input", "--output"]
            } else {
                &["--input", "--exons", "--output"]
            })?;
            let out = Outputs::new(vec![a.path("--output")?])?;
            if verb == "count" {
                counts::count(&a.path("--input")?, &out.temp[0])?;
            } else {
                counts::merge(
                    &a.path("--input")?,
                    &ExonIndex::load(&a.path("--exons")?)?,
                    &out.temp[0],
                )?;
            }
            out.commit()
        }
        "score" => {
            a.check(&[
                "--nc",
                "--tr",
                "--output",
                "--mode",
                "--format",
                "--rows",
                "--reference-id",
                "--summary",
            ])?;
            let mut inputs = vec![a.path("--nc")?];
            for t in a.0.get("--tr").ok_or("At least one --tr required")? {
                inputs.push(PathBuf::from(t));
            }
            let mode = scoring::Mode::parse(a.optional("--mode").unwrap_or("mrna"))?;
            let format = match a.optional("--format").unwrap_or("unified") {
                "unified" => scoring::Format::Unified,
                "legacy" => scoring::Format::Legacy,
                _ => return Err("score --format must be unified|legacy".into()),
            };
            let rows = match a
                .optional("--rows")
                .unwrap_or(if mode == scoring::Mode::Mrna {
                    "candidates"
                } else {
                    "all"
                }) {
                "candidates" => scoring::Rows::Candidates,
                "all" => scoring::Rows::All,
                _ => return Err("score --rows must be candidates|all".into()),
            };
            let mut targets = vec![a.path("--output")?];
            if let Some(path) = a.optional("--summary") {
                targets.push(PathBuf::from(path));
            }
            let out = Outputs::new(targets)?;
            scoring::run(
                &inputs,
                &out.temp[0],
                out.temp.get(1).map(PathBuf::as_path),
                &scoring::Options {
                    mode,
                    format,
                    rows,
                    reference: a.optional("--reference-id").map(String::from),
                },
            )?;
            out.commit()
        }
        "plan" | "run" => {
            a.check(&[
                "--samples",
                "--from",
                "--exons",
                "--out",
                "--index",
                "--bac-index",
                "--rrna-index",
                "--mrna-index",
                "--threads",
                "--jobs",
                "--storage",
                "--rrna-filter",
                "--adapter-r1",
                "--adapter-r2",
                "--qc",
            ])?;
            let prefix = |option: &str, default: &str| {
                a.optional(option)
                    .map(PathBuf::from)
                    .or_else(|| a.optional("--index").map(|d| Path::new(d).join(default)))
            };
            workflow::execute(
                workflow::Config {
                    samples: a.path("--samples")?,
                    from: a.get("--from")?.into(),
                    exons: a.path("--exons")?,
                    out: a.path("--out")?,
                    bac_index: prefix("--bac-index", "bac"),
                    rrna_index: prefix("--rrna-index", "rrna"),
                    mrna_index: prefix("--mrna-index", "mrna"),
                    adapter_r1: a
                        .optional("--adapter-r1")
                        .unwrap_or("AGATCGGAAGAGCACA")
                        .into(),
                    adapter_r2: a
                        .optional("--adapter-r2")
                        .unwrap_or("GATCGTCGGACTGTAG")
                        .into(),
                    threads: a.number("--threads", 4)?,
                    jobs: a.number("--jobs", 1)?,
                    exe: std::env::current_exe().map_err(|e| e.to_string())?,
                    storage: a.optional("--storage").unwrap_or("compact").into(),
                    rrna_filter: a.optional("--rrna-filter").unwrap_or("both").into(),
                    qc: match a.optional("--qc").unwrap_or(
                        if a.optional("--storage") == Some("all") {
                            "no"
                        } else {
                            "yes"
                        },
                    ) {
                        "yes" => true,
                        "no" => false,
                        _ => return Err("qc=yes|no".into()),
                    },
                },
                verb == "plan",
            )
        }
        _ => Err(format!("Unknown command: {verb}\n{HELP}")),
    }
}
fn main() {
    if let Err(e) = main_impl() {
        eprintln!("ERROR: {e}");
        std::process::exit(1);
    }
}
