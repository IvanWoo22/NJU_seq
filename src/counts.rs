//! Native counterparts of the current mRNA count, merge and score scripts.
//!
//! These intentionally retain the existing end-minus-one candidate coordinate,
//! first-observed sequence context and NC-first, all-TR intersection scoring.
use crate::model::{Alignment, ExonIndex, parse_cigar};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

type Result<T> = std::result::Result<T, String>;

fn reader(path: &Path) -> Result<BufReader<File>> {
    File::open(path)
        .map(BufReader::new)
        .map_err(|e| format!("{}: {e}", path.display()))
}

fn writer(path: &Path) -> Result<BufWriter<File>> {
    File::create(path)
        .map(BufWriter::new)
        .map_err(|e| format!("{}: {e}", path.display()))
}

fn integer(value: &str, label: &str, positive: bool) -> Result<u64> {
    if value.is_empty() || !value.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Invalid {label}: {value}"));
    }
    let n = value.parse::<u64>().map_err(|e| format!("{label}: {e}"))?;
    if positive && n == 0 {
        return Err(format!("{label} must be positive"));
    }
    Ok(n)
}

fn add(a: u64, b: u64, label: &str) -> Result<u64> {
    a.checked_add(b).ok_or_else(|| format!("{label} overflow"))
}

/// Last three reference positions: insertions consume only query bases;
/// deletions contribute N. This matches count.pl's repaired ALIGNMENT_TAIL.
fn alignment_tail(cigar: &str, sequence: &str) -> Result<(u64, [u8; 3])> {
    if cigar.is_empty() || !sequence.is_ascii() {
        return Err("Empty CIGAR or non-ASCII sequence".into());
    }
    let mut query = 0usize;
    let mut reference = 0u64;
    let mut tail = [0u8; 3];
    let mut tail_len = 0usize;
    let sequence = sequence.as_bytes();
    for (length, operation) in parse_cigar(cigar)? {
        let number = usize::try_from(length).map_err(|_| "CIGAR length overflow")?;
        let ch = operation as u8;
        if !matches!(ch, b'=' | b'X' | b'I' | b'D') {
            return Err(format!("Invalid --xeq CIGAR: {cigar}"));
        }
        let query_end = if ch != b'D' {
            let end = query.checked_add(number).ok_or("Query length overflow")?;
            if end > sequence.len() {
                return Err("CIGAR/sequence length mismatch".into());
            }
            end
        } else {
            query
        };
        if ch != b'I' {
            reference = add(reference, number as u64, "Reference length")?;
            // Only three positions are needed, even for very long operations.
            let keep = number.min(3);
            let start = number - keep;
            if number >= 3 {
                tail_len = 0;
            }
            for offset in start..number {
                let base = if ch == b'D' {
                    b'N'
                } else {
                    sequence[query + offset]
                };
                if tail_len < 3 {
                    tail[tail_len] = base;
                    tail_len += 1;
                } else {
                    tail.rotate_left(1);
                    tail[2] = base;
                }
            }
        }
        query = query_end;
    }
    if query != sequence.len() || reference < 3 || tail_len != 3 {
        return Err("CIGAR/sequence length mismatch or fragment too short".into());
    }
    Ok((reference, tail))
}

#[derive(Debug)]
struct TranscriptCount {
    context: [u8; 3],
    ends: u64,
}

/// Seven columns, lexicographic transcript order and numeric position order.
pub fn count(input: &Path, output: &Path) -> Result<()> {
    let mut ends: HashMap<(String, u64), TranscriptCount> = HashMap::new();
    let mut starts: HashMap<(String, u64), u64> = HashMap::new();
    for (line_no, line) in reader(input)?.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let a = Alignment::parse(&line)
            .map_err(|e| format!("{}:{}: {e}", input.display(), line_no + 1))?;
        let (length, context) = alignment_tail(&a.cigar, &a.sequence)
            .map_err(|e| format!("{}:{}: {e}", input.display(), line_no + 1))?;
        let position = add(a.start, length - 2, "Candidate position")?;
        let start = starts.entry((a.reference.clone(), a.start)).or_default();
        *start = add(*start, 1, "StartCount")?;
        let entry = ends
            .entry((a.reference, position))
            .or_insert(TranscriptCount { context, ends: 0 });
        entry.ends = add(entry.ends, 1, "EndCount")?;
    }
    let mut rows: Vec<_> = ends.into_iter().collect();
    rows.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut out = writer(output)?;
    for ((transcript, position), row) in rows {
        let start = starts
            .get(&(transcript.clone(), position))
            .copied()
            .unwrap_or(0);
        writeln!(
            out,
            "{transcript}\t{position}\t{}\t{}\t{}\t{start}\t{}",
            row.context[0] as char, row.context[1] as char, row.context[2] as char, row.ends
        )
        .map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct Site {
    chrom: String,
    position: u64,
    strand: char,
}

#[derive(Debug)]
struct GenomicCount {
    context: [String; 3],
    genes: Vec<String>,
    starts: u64,
    ends: u64,
}

/// Map transcript counts onto genomic positions and sum isoform contributions.
/// Bases remain transcript-oriented on both strands, as in merge.pl.
pub fn merge(input: &Path, index: &ExonIndex, output: &Path) -> Result<()> {
    let mut sites: HashMap<Site, GenomicCount> = HashMap::new();
    for (line_no, line) in reader(input)?.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 7 || fields[2..5].iter().any(|v| v.len() != 1) {
            return Err(format!(
                "{}:{}: expected seven-column count row",
                input.display(),
                line_no + 1
            ));
        }
        // UMI handling and the Perl merge stage accept a reference identifier
        // followed by a pipe-delimited description. Preserve the original text
        // in count output but resolve its identifier for annotation lookup.
        let transcript_id = fields[0]
            .split(|c: char| c == '|' || c.is_whitespace())
            .next()
            .unwrap_or("");
        let transcript = index
            .transcripts
            .get(fields[0])
            .or_else(|| index.transcripts.get(transcript_id))
            .ok_or_else(|| format!("Transcript absent from annotation: {}", fields[0]))?;
        let position = integer(fields[1], "Transcript position", true)?;
        let starts = integer(fields[5], "StartCount", false)?;
        let ends = integer(fields[6], "EndCount", true)?;
        let site = Site {
            chrom: transcript.chrom.clone(),
            position: transcript.map_position(position)?,
            strand: transcript.strand,
        };
        if let Some(row) = sites.get_mut(&site) {
            if !row.genes.contains(&transcript.gene) {
                row.genes.push(transcript.gene.clone());
            }
            row.starts = add(row.starts, starts, "Merged StartCount")?;
            row.ends = add(row.ends, ends, "Merged EndCount")?;
        } else {
            sites.insert(
                site,
                GenomicCount {
                    context: [fields[2].into(), fields[3].into(), fields[4].into()],
                    genes: vec![transcript.gene.clone()],
                    starts,
                    ends,
                },
            );
        }
    }
    let mut rows: Vec<_> = sites.into_iter().collect();
    rows.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut out = writer(output)?;
    for (site, row) in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            site.chrom,
            site.position,
            site.strand,
            row.context[0],
            row.context[1],
            row.context[2],
            row.genes.join("/"),
            row.starts,
            row.ends
        )
        .map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}

/// Compatibility adapter for callers using the historical mRNA table schema.
pub fn score(inputs: &[PathBuf], output: &Path) -> Result<()> {
    crate::scoring::run(
        inputs,
        output,
        None,
        &crate::scoring::Options {
            mode: crate::scoring::Mode::Mrna,
            format: crate::scoring::Format::Legacy,
            rows: crate::scoring::Rows::Candidates,
            reference: None,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Transcript;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);
    struct Fixture(PathBuf);
    impl Fixture {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "nju-rust-counts-{}-{}",
                std::process::id(),
                NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn file(&self, name: &str, text: &str) -> PathBuf {
            let path = self.0.join(name);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(&path, text).unwrap();
            path
        }
        fn output(&self) -> PathBuf {
            self.0.join("output.tsv")
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn tail_tracks_reference_positions_through_indels() {
        assert_eq!(alignment_tail("3=1I2=", "ACGTAC").unwrap(), (5, *b"GAC"));
        assert_eq!(alignment_tail("3=1D1=", "ACGT").unwrap(), (5, *b"GNT"));
        assert_eq!(alignment_tail("3=3D", "ACG").unwrap(), (6, *b"NNN"));
        assert_eq!(alignment_tail("4=2I", "ACGTAA").unwrap(), (4, *b"CGT"));
        assert_eq!(alignment_tail("3=1D2=", "ACGTA").unwrap(), (6, *b"NTA"));
    }

    #[test]
    fn tail_rejects_non_xeq_cigars_and_truncation() {
        for (cigar, seq) in [
            ("3M", "ACG"),
            ("03=", "ACG"),
            ("2=", "AC"),
            ("3=1", "ACG"),
            ("4=", "ACG"),
            ("3=", "ACGA"),
        ] {
            assert!(
                alignment_tail(cigar, seq).is_err(),
                "accepted {cigar}/{seq}"
            );
        }
    }

    #[test]
    fn count_preserves_first_context_and_collects_start_count() {
        let f = Fixture::new();
        let input = f.file("input.tsv", "r1\tAT1.1\t2\t3=\tACG\nr2\tAT1.1\t2\t3=\tTTT\nr3\tAT1.1\t3\t3=\tGCA\nr4\tAT1.1\t10\t3=\tATC\n");
        count(&input, &f.output()).unwrap();
        assert_eq!(
            fs::read_to_string(f.output()).unwrap(),
            "AT1.1\t3\tA\tC\tG\t1\t2\nAT1.1\t4\tG\tC\tA\t0\t1\nAT1.1\t11\tA\tT\tC\t0\t1\n"
        );
    }

    #[test]
    fn merge_maps_exon_junctions_both_strands_and_sums_isoforms() {
        let f = Fixture::new();
        let mut index = ExonIndex {
            transcripts: HashMap::new(),
        };
        index.transcripts.insert(
            "AT1.1".into(),
            Transcript::new("1".into(), '+', "AT1".into(), vec![(100, 102), (200, 202)]).unwrap(),
        );
        index.transcripts.insert(
            "AT1.2".into(),
            Transcript::new("1".into(), '+', "AT1".into(), vec![(100, 102), (200, 202)]).unwrap(),
        );
        index.transcripts.insert(
            "AT10.1".into(),
            Transcript::new("1".into(), '+', "AT10".into(), vec![(200, 202)]).unwrap(),
        );
        index.transcripts.insert(
            "AT2.1".into(),
            Transcript::new("1".into(), '-', "AT2".into(), vec![(100, 102), (200, 202)]).unwrap(),
        );
        let input = f.file("input.tsv", "AT1.1\t4\tA\tC\tG\t1\t2\nAT1.2\t4\tT\tT\tT\t2\t3\nAT10.1\t1\tG\tG\tG\t3\t4\nAT2.1\t4\tG\tC\tA\t0\t1\nAT2.1\t3\tA\tT\tG\t1\t2\n");
        merge(&input, &index, &f.output()).unwrap();
        assert_eq!(
            fs::read_to_string(f.output()).unwrap(),
            "1\t102\t-\tG\tC\tA\tAT2\t0\t1\n1\t200\t+\tA\tC\tG\tAT1/AT10\t6\t9\n1\t200\t-\tA\tT\tG\tAT2\t1\t2\n"
        );
    }

    #[test]
    fn merge_resolves_reference_descriptions() {
        let f = Fixture::new();
        let index = ExonIndex {
            transcripts: HashMap::from([(
                "AT1.1".into(),
                Transcript::new("1".into(), '+', "AT1".into(), vec![(100, 102)]).unwrap(),
            )]),
        };
        let input = f.file(
            "input.tsv",
            "AT1.1|description\t1\tA\tC\tG\t0\t2\nAT1.1 description\t1\tA\tC\tG\t1\t3\n",
        );
        merge(&input, &index, &f.output()).unwrap();
        assert_eq!(
            fs::read_to_string(f.output()).unwrap(),
            "1\t100\t+\tA\tC\tG\tAT1\t1\t5\n"
        );
    }

    #[test]
    fn score_has_multi_tr_intersection_and_historical_thresholds() {
        let f = Fixture::new();
        let nc = f.file(
            "NC/counts.tsv",
            "1\t10\t+\tA\tC\tG\tAT1\t0\t1\n1\t30\t+\tA\tC\tG\tAT1\t0\t100\n",
        );
        let tr1 = f.file("TR1/counts.tsv", "1\t10\t+\tA\tC\tG\tAT1\t0\t40\n1\t20\t+\tA\tC\tG\tAT1\t0\t80\n1\t30\t+\tA\tC\tG\tAT1\t0\t1\n");
        let tr2 = f.file(
            "TR2/counts.tsv",
            "1\t10\t+\tA\tC\tG\tAT1\t0\t50\n1\t30\t+\tA\tC\tG\tAT1\t0\t1\n",
        );
        score(&[nc, tr1, tr2], &f.output()).unwrap();
        let output = fs::read_to_string(f.output()).unwrap();
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].split('\t').count(), 18);
        assert_eq!(
            lines[1],
            "1\t10\t+\tA\tC\tG\tAT1\t1\t40\t34\t39\t50\t100\t49\t134\t67\t88\t44"
        );
    }

    #[test]
    fn empty_tr_is_valid_and_duplicate_sites_are_rejected() {
        let f = Fixture::new();
        let nc = f.file("NC/counts.tsv", "");
        let tr = f.file("TR/counts.tsv", "");
        score(&[nc.clone(), tr], &f.output()).unwrap();
        assert_eq!(fs::read_to_string(f.output()).unwrap().lines().count(), 1);
        let duplicate = f.file(
            "dup.tsv",
            "1\t10\t+\tA\tC\tG\tAT1\t0\t1\n1\t10\t+\tA\tC\tG\tAT1\t0\t1\n",
        );
        assert!(
            score(&[nc, duplicate], &f.output())
                .unwrap_err()
                .contains("duplicate")
        );
    }
}
