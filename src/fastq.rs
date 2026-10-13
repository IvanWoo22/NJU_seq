//! Strict streaming FASTQ I/O. gzip is a checked child process, never a shell string.
use crate::model::Alignment;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::{Child, Command, Stdio};

pub struct Input {
    reader: BufReader<Box<dyn Read>>,
    child: Option<Child>,
}
impl Input {
    pub fn open(path: &Path) -> Result<Self, String> {
        if path.extension().is_some_and(|s| s == "gz") {
            let mut child = Command::new("gzip")
                .arg("-dc")
                .arg("--")
                .arg(path)
                .stdout(Stdio::piped())
                .spawn()
                .map_err(|e| e.to_string())?;
            let reader = BufReader::with_capacity(
                1024 * 1024,
                Box::new(child.stdout.take().unwrap()) as Box<dyn Read>,
            );
            Ok(Self {
                reader,
                child: Some(child),
            })
        } else {
            Ok(Self {
                reader: BufReader::with_capacity(
                    1024 * 1024,
                    Box::new(File::open(path).map_err(|e| e.to_string())?),
                ),
                child: None,
            })
        }
    }
    pub fn finish(mut self) -> Result<(), String> {
        if let Some(mut child) = self.child.take() {
            let status = child.wait().map_err(|e| e.to_string())?;
            if !status.success() {
                return Err(format!("gzip input failed: {status}"));
            }
        }
        Ok(())
    }
}
impl Read for Input {
    fn read(&mut self, b: &mut [u8]) -> io::Result<usize> {
        self.reader.read(b)
    }
}
impl BufRead for Input {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.reader.fill_buf()
    }
    fn consume(&mut self, n: usize) {
        self.reader.consume(n)
    }
}
impl Drop for Input {
    fn drop(&mut self) {
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

pub(crate) struct Output {
    writer: Option<BufWriter<Box<dyn Write>>>,
    child: Option<Child>,
}
impl Output {
    pub(crate) fn open(path: &Path, gzip: bool) -> Result<Self, String> {
        let file = File::create(path).map_err(|e| e.to_string())?;
        if gzip {
            let mut child = Command::new("gzip")
                .args(["-c", "-3"])
                .stdin(Stdio::piped())
                .stdout(Stdio::from(file))
                .spawn()
                .map_err(|e| e.to_string())?;
            let writer = BufWriter::with_capacity(
                1024 * 1024,
                Box::new(child.stdin.take().unwrap()) as Box<dyn Write>,
            );
            Ok(Self {
                writer: Some(writer),
                child: Some(child),
            })
        } else {
            Ok(Self {
                writer: Some(BufWriter::new(Box::new(file))),
                child: None,
            })
        }
    }
    pub(crate) fn record(&mut self, r: &Record) -> Result<(), String> {
        writeln!(
            self.writer.as_mut().unwrap(),
            "{}\n{}\n{}\n{}",
            r.header,
            r.sequence,
            r.plus,
            r.quality
        )
        .map_err(|e| e.to_string())
    }
    pub(crate) fn finish(mut self) -> Result<(), String> {
        if let Some(mut writer) = self.writer.take() {
            writer.flush().map_err(|e| e.to_string())?;
            drop(writer);
        }
        if let Some(mut child) = self.child.take() {
            let status = child.wait().map_err(|e| e.to_string())?;
            if !status.success() {
                return Err(format!("gzip output failed: {status}"));
            }
        }
        Ok(())
    }
}
impl Drop for Output {
    fn drop(&mut self) {
        self.writer.take();
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

#[derive(Debug)]
pub struct Record {
    pub header: String,
    pub sequence: String,
    pub plus: String,
    pub quality: String,
}
impl Record {
    pub fn id(&self) -> &str {
        self.header[1..].split_whitespace().next().unwrap_or("")
    }
}

/// Normalize only the conventional terminal mate suffix. Underscores and UMI
/// strings remain part of the read name; they never imply a library design.
pub(crate) fn fragment_id(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

pub fn record<R: BufRead>(input: &mut R) -> Result<Option<Record>, String> {
    let mut fields = Vec::with_capacity(4);
    for i in 0..4 {
        let mut s = String::new();
        let n = input.read_line(&mut s).map_err(|e| e.to_string())?;
        if n == 0 {
            if i == 0 {
                return Ok(None);
            }
            return Err("Truncated FASTQ record".into());
        }
        fields.push(s.trim_end_matches(['\n', '\r']).to_owned());
    }
    let mut fields = fields.into_iter();
    let r = Record {
        header: fields.next().unwrap(),
        sequence: fields.next().unwrap(),
        plus: fields.next().unwrap(),
        quality: fields.next().unwrap(),
    };
    if !r.header.starts_with('@')
        || r.id().is_empty()
        || !r.plus.starts_with('+')
        || r.sequence.is_empty()
        || r.sequence.len() != r.quality.len()
    {
        return Err(format!("Invalid FASTQ record: {}", r.header));
    }
    Ok(Some(r))
}

pub fn almost_unique(input: &Path, fastq: &Path, output: &Path) -> Result<(u64, u64), String> {
    let mut fq = Input::open(fastq)?;
    let mut sequences = HashMap::new();
    while let Some(r) = record(&mut fq)? {
        let id = fragment_id(r.id()).to_owned();
        if id.is_empty() {
            return Err("Empty FASTQ fragment ID after mate suffix normalization".into());
        }
        if sequences.insert(id.clone(), r.sequence).is_some() {
            return Err(format!(
                "Repeated FASTQ fragment ID after mate suffix normalization: {id}"
            ));
        }
    }
    fq.finish()?;
    let mut matches: HashMap<String, u64> = HashMap::new();
    let mut n = 0;
    let mut stream = Input::open(input)?;
    let mut line = String::new();
    while stream.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        let a = Alignment::parse(line.trim_end_matches(['\n', '\r']))?;
        n += 1;
        let id = fragment_id(&a.read);
        if sequences.get(id) == Some(&a.sequence) {
            *matches.entry(id.to_owned()).or_default() += 1;
        }
        line.clear();
    }
    stream.finish()?;
    let mut stream = Input::open(input)?;
    let mut out = BufWriter::new(File::create(output).map_err(|e| e.to_string())?);
    let mut kept = 0;
    while stream.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        let a = Alignment::parse(line.trim_end_matches(['\n', '\r']))?;
        let id = fragment_id(&a.read);
        if matches.get(id) == Some(&1) && sequences.get(id) == Some(&a.sequence) {
            writeln!(out, "{}", a.to_tsv()).map_err(|e| e.to_string())?;
            kept += 1;
        }
        line.clear();
    }
    stream.finish()?;
    out.flush().map_err(|e| e.to_string())?;
    Ok((n, kept))
}

pub fn delete_pairs(
    names: &Path,
    r1: &Path,
    r2: &Path,
    out1: &Path,
    out2: &Path,
    gzip: bool,
) -> Result<(u64, u64), String> {
    let mut names_input = Input::open(names)?;
    let mut names_set = HashSet::new();
    let mut line = String::new();
    while names_input
        .read_line(&mut line)
        .map_err(|e| e.to_string())?
        > 0
    {
        let id = line
            .trim_end_matches(['\n', '\r'])
            .split('\t')
            .next()
            .unwrap_or("");
        if !id.is_empty() {
            let id = fragment_id(id);
            if id.is_empty() {
                return Err("Empty deletion fragment ID after mate suffix normalization".into());
            }
            names_set.insert(id.to_owned());
        }
        line.clear();
    }
    names_input.finish()?;
    let mut one = Input::open(r1)?;
    let mut two = Input::open(r2)?;
    let mut first = Output::open(out1, gzip)?;
    let mut second = Output::open(out2, gzip)?;
    let (mut total, mut kept) = (0, 0);
    loop {
        match (record(&mut one)?, record(&mut two)?) {
            (None, None) => break,
            (Some(a), Some(b)) => {
                let id = fragment_id(a.id());
                if id.is_empty() || id != fragment_id(b.id()) {
                    return Err(format!("Pair ID mismatch: {} / {}", a.id(), b.id()));
                }
                total += 1;
                if !names_set.contains(id) {
                    first.record(&a)?;
                    second.record(&b)?;
                    kept += 1;
                }
            }
            _ => return Err("R1/R2 record counts differ".into()),
        }
    }
    one.finish()?;
    two.finish()?;
    first.finish()?;
    second.finish()?;
    Ok((total, kept))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);
    struct Fixture(PathBuf);
    impl Fixture {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "nju-fastq-mate-{}-{}",
                std::process::id(),
                NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn path(&self, name: &str) -> PathBuf {
            self.0.join(name)
        }
        fn write(&self, name: &str, content: &str) -> PathBuf {
            let path = self.path(name);
            fs::write(&path, content).unwrap();
            path
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    #[test]
    fn fastq_truncation_and_lengths() {
        assert!(record(&mut &b"@a\nAC\n+\n"[..]).is_err());
        assert!(record(&mut &b"@a\nAC\n+\nI\n"[..]).is_err());
        assert!(record(&mut &b""[..]).unwrap().is_none());
    }
    #[test]
    fn fastq_header_whitespace() {
        let r = record(&mut &b"@a desc\nAC\n+\nII\n"[..]).unwrap().unwrap();
        assert_eq!(r.id(), "a");
    }

    #[test]
    fn only_terminal_mate_suffix_is_normalized() {
        assert_eq!(fragment_id("read/1"), "read");
        assert_eq!(fragment_id("read/2"), "read");
        assert_eq!(fragment_id("read_AAAAAAAA/1"), "read_AAAAAAAA");
        for name in ["read_AAAAAAAA", "read/1_AAAAAAAA", "read/3", "read_name"] {
            assert_eq!(fragment_id(name), name);
        }
    }

    #[test]
    fn deletion_normalizes_sam_or_mate_names_and_preserves_original_records() {
        let fixture = Fixture::new();
        let kept_r1 = "@keep_AAAAAAAA/1 description\nACGT\n+\nIJKL\n";
        let kept_r2 = "@keep_AAAAAAAA/2 description\nTGCA\n+\nLKJI\n";
        let one = fixture.write("R1.fq", &format!("@drop/1\nACGT\n+\nIIII\n{kept_r1}"));
        let two = fixture.write("R2.fq", &format!("@drop/2\nACGT\n+\nIIII\n{kept_r2}"));
        for (number, name) in ["drop", "drop/1", "drop/2"].iter().enumerate() {
            let names = fixture.write("names.tsv", &format!("{name}\tignored_fields\n"));
            let out1 = fixture.path(&format!("{number}.R1.fq"));
            let out2 = fixture.path(&format!("{number}.R2.fq"));
            assert_eq!(
                delete_pairs(&names, &one, &two, &out1, &out2, false).unwrap(),
                (2, 1)
            );
            assert_eq!(fs::read_to_string(out1).unwrap(), kept_r1);
            assert_eq!(fs::read_to_string(out2).unwrap(), kept_r2);
        }
    }

    #[test]
    fn almost_unique_normalizes_queries_and_preserves_alignment_read_ids() {
        let fixture = Fixture::new();
        let fq = fixture.write(
            "R1.fq",
            "@keep_AAAAAAAA/1 description\nACGT\n+\nIIII\n@multi/1\nACGT\n+\nIIII\n",
        );
        let row = "keep_AAAAAAAA\tAT1.1\t1\t4=\tACGT\n";
        // Two aliases of the same read are still two compatible alignments.
        let input = fixture.write(
            "input.tsv",
            &format!("{row}multi\tAT1.1\t1\t4=\tACGT\nmulti/1\tAT1.1\t2\t4=\tACGT\n"),
        );
        let out = fixture.path("output.tsv");
        assert_eq!(almost_unique(&input, &fq, &out).unwrap(), (3, 1));
        assert_eq!(fs::read_to_string(out).unwrap(), row);
        let fq = fixture.write("plain.fq", "@keep_AAAAAAAA\nACGT\n+\nIIII\n");
        let row = "keep_AAAAAAAA/1\tAT1.1\t1\t4=\tACGT\n";
        let input = fixture.write("reverse.tsv", row);
        let out = fixture.path("reverse.output.tsv");
        assert_eq!(almost_unique(&input, &fq, &out).unwrap(), (1, 1));
        assert_eq!(fs::read_to_string(out).unwrap(), row);
    }

    #[test]
    fn almost_unique_rejects_normalization_collisions_in_fastq() {
        let fixture = Fixture::new();
        let one = fixture.write("R1.fq", "@same\nACGT\n+\nIIII\n@same/1\nACGT\n+\nIIII\n");
        let alignments = fixture.write("input.tsv", "same\tAT1.1\t1\t4=\tACGT\n");
        let error = almost_unique(&alignments, &one, &fixture.path("out.tsv")).unwrap_err();
        assert!(error.contains("Repeated FASTQ fragment ID"), "{error}");
    }
}
