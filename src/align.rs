//! Bowtie2 -> existing native selection, with no intermediate SAM by default.
//! --reorder makes the FASTQ join a streaming merge instead of a whole-library hash.
use crate::{
    fastq, filter,
    model::{Alignment, ExonIndex},
};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};

pub struct Options<'a> {
    pub mode: &'a str,
    pub index: &'a Path,
    pub r1: &'a Path,
    pub r2: &'a Path,
    pub exons: Option<&'a Path>,
    pub threads: usize,
    pub output: &'a Path,
    pub log: &'a Path,
    pub sam: Option<&'a Path>,
    pub screen_multimap: Option<usize>,
    pub exact_pair_screen: bool,
}

pub fn arguments(mode: &str, threads: usize) -> Result<Vec<String>, String> {
    if !["bac", "rrna", "mrna"].contains(&mode) || threads == 0 {
        return Err("align-select: mode=bac|rrna|mrna; threads must be positive".into());
    }
    let mut args: Vec<String> = [
        "--end-to-end",
        "--reorder",
        "--seed",
        "0",
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
    ]
    .map(String::from)
    .to_vec();
    args.extend(["-p".into(), threads.to_string()]);
    if mode == "mrna" {
        args.push("-a".into());
    } else {
        args.extend(["-k".into(), "5".into()]);
    }
    if mode == "rrna" {
        args.extend(["--np".into(), "0".into()]);
    } else {
        args.extend(
            [
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
            ]
            .map(String::from),
        );
    }
    Ok(args)
}

/// Receives already quality/multimatch/coordinate-filtered lines. A group is
/// accepted only if exactly one alignment has SEQ identical to original R1.
struct SelectedWriter<W: Write> {
    out: W,
    fq: Option<fastq::Input>,
    pending: Vec<u8>,
    current: String,
    original_sequence: String,
    first_match: Option<String>,
    matches: u64,
    pub kept: u64,
    reject_multiple: bool,
    fq2: Option<fastq::Input>,
    valid_complement: bool,
}
impl<W: Write> SelectedWriter<W> {
    fn new(out: W, r1: Option<&Path>) -> Result<Self, String> {
        Ok(Self {
            out,
            fq: r1.map(fastq::Input::open).transpose()?,
            pending: vec![],
            current: String::new(),
            original_sequence: String::new(),
            first_match: None,
            matches: 0,
            kept: 0,
            reject_multiple: false,
            fq2: None,
            valid_complement: true,
        })
    }
    fn flush_group(&mut self) -> Result<(), String> {
        if self.reject_multiple {
            if self.matches >= 2 && self.valid_complement {
                writeln!(self.out, "{}", self.current).map_err(|e| e.to_string())?;
                self.kept += 1;
            }
            self.first_match = None;
            self.matches = 0;
            return Ok(());
        }
        if self.fq.is_some() && self.matches == 1 {
            writeln!(self.out, "{}", self.first_match.take().unwrap())
                .map_err(|e| e.to_string())?;
            self.kept += 1;
        }
        self.first_match = None;
        self.matches = 0;
        Ok(())
    }
    fn line(&mut self, line: &str) -> Result<(), String> {
        let a = Alignment::parse(line)?;
        let id = fastq::fragment_id(&a.read);
        if id != self.current {
            self.flush_group()?;
            if let Some(fq) = &mut self.fq {
                loop {
                    let r = fastq::record(fq)?.ok_or_else(|| format!(
                        "Ordered FASTQ join could not find {id}; SAM must come from --reorder with the same R1"))?;
                    if let Some(fq2) = &mut self.fq2 {
                        let b = fastq::record(fq2)?.ok_or("Truncated R2 in exact-pair proof")?;
                        if fastq::fragment_id(r.id()) != fastq::fragment_id(b.id()) {
                            return Err("Exact-pair proof requires synchronized mates".into());
                        }
                        self.valid_complement = r.sequence.len() <= 120
                            && r.sequence.len() == b.sequence.len()
                            && r.sequence.bytes().zip(b.sequence.bytes().rev()).all(
                                |(a, b)| match a {
                                    b'A' => b == b'T',
                                    b'T' => b == b'A',
                                    b'C' => b == b'G',
                                    b'G' => b == b'C',
                                    _ => false,
                                },
                            );
                    }
                    if fastq::fragment_id(r.id()) == id {
                        self.original_sequence = r.sequence;
                        break;
                    }
                }
            } else {
                writeln!(self.out, "{id}").map_err(|e| e.to_string())?;
                self.kept += 1;
            }
            self.current = id.to_owned();
        }
        if self.fq.is_some() && a.sequence == self.original_sequence {
            self.matches += 1;
            if self.matches == 1 {
                self.first_match = Some(line.to_owned());
            }
        }
        Ok(())
    }
    fn finish(mut self) -> Result<u64, String> {
        if !self.pending.is_empty() {
            return Err("Unterminated selected alignment".into());
        }
        self.flush_group()?;
        if let Some(mut fq) = self.fq.take() {
            while let Some(a) = fastq::record(&mut fq)? {
                if let Some(two) = &mut self.fq2 {
                    let b = fastq::record(two)?.ok_or("Truncated R2 in exact-pair proof")?;
                    if fastq::fragment_id(a.id()) != fastq::fragment_id(b.id()) {
                        return Err("Exact-pair proof mate mismatch".into());
                    }
                }
            }
            fq.finish()?;
        }
        if let Some(mut two) = self.fq2.take() {
            if fastq::record(&mut two)?.is_some() {
                return Err("Excess R2 in exact-pair proof".into());
            }
            two.finish()?;
        }
        self.out.flush().map_err(|e| e.to_string())?;
        Ok(self.kept)
    }
}
impl<W: Write> Write for SelectedWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        for &b in bytes {
            if b == b'\n' {
                let line = String::from_utf8(std::mem::take(&mut self.pending))
                    .map_err(io::Error::other)?;
                self.line(&line).map_err(io::Error::other)?;
            } else {
                self.pending.push(b);
            }
        }
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        self.out.flush()
    }
}

fn select<R: BufRead>(reader: R, options: &Options<'_>) -> Result<String, String> {
    let strategy = match options.mode {
        "bac" => filter::Strategy::Bac,
        "rrna" => filter::Strategy::Rrna,
        "mrna" => filter::Strategy::Mrna,
        _ => return Err("Unknown alignment mode".into()),
    };
    let index = if options.mode == "mrna" {
        Some(ExonIndex::load(
            options.exons.ok_or("mRNA requires --exons")?,
        )?)
    } else {
        None
    };
    let out = BufWriter::new(File::create(options.output).map_err(|e| e.to_string())?);
    let mut selected = SelectedWriter::new(
        out,
        if options.mode == "mrna" {
            Some(options.r1)
        } else {
            None
        },
    )?;
    selected.reject_multiple = options.screen_multimap.is_some();
    if options.exact_pair_screen {
        selected.fq2 = Some(fastq::Input::open(options.r2)?);
    }
    let stats = filter::run(
        reader,
        &mut selected,
        if options.exact_pair_screen {
            filter::InputFormat::ExactSingleEndSam
        } else {
            filter::InputFormat::Sam
        },
        strategy,
        index.as_ref(),
    )?;
    let kept = selected.finish()?;
    let stats_json = format!(
        "{{\"sam_records\":{},\"sam_selected\":{},\"quality_passed\":{},\"reference_unique\":{},\"coordinate_unique\":{},\"read_groups\":{},\"selected_reads\":{},\"sam_materialized\":{}}}\n",
        stats.input_records,
        stats.sam_selected,
        stats.quality_passed,
        stats.reference_unique,
        stats.coordinate_unique,
        stats.read_groups,
        kept,
        options.sam.is_some()
    );
    Ok(stats_json.trim_end().trim_end_matches('}').to_owned()
        + &format!(
            ",\"screen_multimap\":{},\"screen_kind\":\"{}\"}}\n",
            options
                .screen_multimap
                .map(|x| x.to_string())
                .unwrap_or("null".into()),
            if options.screen_multimap.is_none() {
                "none"
            } else if options.exact_pair_screen {
                "exact-pair"
            } else {
                "paired"
            }
        ))
}

pub fn run(options: Options<'_>) -> Result<String, String> {
    let mut args = arguments(options.mode, options.threads)?;
    if let Some(k) = options.screen_multimap {
        if options.mode != "mrna" || k < 2 {
            return Err("Multimap screen requires mrna and k>=2".into());
        }
        args.retain(|arg| arg != "-a");
        args.extend(["-k".into(), k.to_string()]);
    }
    if options.exact_pair_screen {
        if options.screen_multimap.is_none() {
            return Err("Exact-pair proof requires --screen-multimap".into());
        }
        let i = args.iter().position(|x| x == "--score-min").unwrap();
        args[i + 1] = "L,0,0".into();
        args.push("--norc".into());
    }
    let mut command = Command::new("bowtie2");
    command
        .args(args)
        .arg("-x")
        .arg(options.index)
        .stderr(Stdio::from(
            File::create(options.log).map_err(|e| e.to_string())?,
        ));
    if options.exact_pair_screen {
        command.arg("-U").arg(options.r1);
    } else {
        command.arg("-1").arg(options.r1).arg("-2").arg(options.r2);
    }
    if let Some(sam) = options.sam {
        let status = command
            .arg("-S")
            .arg(sam)
            .stdout(Stdio::null())
            .status()
            .map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("Bowtie2 failed: {status}"));
        }
        return select(
            BufReader::with_capacity(1024 * 1024, File::open(sam).map_err(|e| e.to_string())?),
            &options,
        );
    }
    let mut child = command
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    let result = select(
        BufReader::with_capacity(1024 * 1024, child.stdout.take().unwrap()),
        &options,
    );
    if result.is_err() {
        let _ = child.kill();
    }
    let status = child.wait().map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!(
            "Bowtie2 failed: {status}; native selection: {result:?}"
        ));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parameters_keep_scientific_settings_and_order() {
        for mode in ["bac", "rrna", "mrna"] {
            let a = arguments(mode, 2).unwrap();
            assert!(a.iter().any(|x| x == "--reorder"));
            assert_eq!(a.iter().any(|x| x == "-a"), mode == "mrna");
            assert_eq!(a.iter().any(|x| x == "-k"), mode != "mrna");
        }
    }
    #[test]
    fn rrna_multiple_alignments_emit_one_fragment_id() {
        let mut bytes = vec![];
        let mut w = SelectedWriter::new(&mut bytes, None).unwrap();
        write!(
            w,
            "a/1\tt\t1\t3=\tAAA\na/1\tu\t1\t3=\tAAA\nb\tt\t1\t3=\tAAA\n"
        )
        .unwrap();
        assert_eq!(w.finish().unwrap(), 2);
        assert_eq!(bytes, b"a\nb\n");
    }
    #[test]
    fn ordered_join_retains_exactly_one_sequence_match_and_checks_remaining_fastq() {
        let path = std::env::temp_dir().join(format!("nju-ordered-{}.fq", std::process::id()));
        std::fs::write(
            &path,
            "@a/1\nAAA\n+\nIII\n@b/1\nCCC\n+\nIII\n@c/1\nGGG\n+\nIII\n@d/1\nATA\n+\nIII\n",
        )
        .unwrap();
        let mut bytes = vec![];
        let mut writer = SelectedWriter::new(&mut bytes, Some(&path)).unwrap();
        write!(
            writer,
            "a\tt\t1\t3=\tAAA\na\tu\t1\t3=\tAAA\nb\tt\t1\t3=\tTTT\nc\tt\t1\t3=\tGGG\n"
        )
        .unwrap();
        assert_eq!(writer.finish().unwrap(), 1);
        assert_eq!(bytes, b"c\tt\t1\t3=\tGGG\n");
        let mut writer = SelectedWriter::new(Vec::new(), Some(&path)).unwrap();
        assert!(write!(writer, "c\tt\t1\t3=\tGGG\nb\tt\t1\t3=\tCCC\n").is_err());
        std::fs::remove_file(path).unwrap();
    }
}
