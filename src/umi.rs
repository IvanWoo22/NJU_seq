//! Exact UMI families after almost-unique mapping, before site counting.
//! UMI matching is exact; directional correction is not applied.

use crate::model::{ExonIndex, parse_cigar};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Context {
    chrom: String,
    strand: char,
    blocks: Vec<(u64, u64)>,
}

#[derive(Debug)]
struct Family {
    count: u64,
    errors: u64,
    representative: String,
    row: String,
    first_transcript: String,
    first_sequence: String,
    multiple_transcripts: bool,
    multiple_sequences: bool,
}

#[derive(Debug)]
struct Bucket {
    site: u64,
    families: BTreeMap<String, Family>,
}

#[derive(Debug)]
struct Ambiguous {
    context: Context,
    site: u64,
    umi: String,
    row: String,
}

#[derive(Debug, Default)]
pub struct UmiSummary {
    pub sample: String,
    pub umi_length: usize,
    pub input_reads: u64,
    pub output_molecules: u64,
    pub removed_reads: u64,
    pub exact_acgt_families: u64,
    pub multi_transcript_families: u64,
    pub mixed_insert_sequence_families: u64,
    pub umi_corrected_families: u64,
    pub acgt_fragment_contexts: u64,
    pub ambiguous_umi_reads_kept: u64,
    pub exact_molecules_including_ambiguous: u64,
    pub additional_molecules_collapsed_by_correction: u64,
    pub max_exact_umis_in_one_context: u64,
    pub family_size_histogram: BTreeMap<u64, u64>,
    pub observed_umis_per_family_histogram: BTreeMap<u64, u64>,
    pub umi_base_counts_by_position: Vec<BTreeMap<char, u64>>,
    pub top_umis: Vec<(String, u64)>,
    pub observed_umi_strings: u64,
}

impl UmiSummary {
    /// Business metrics only. The evaluation runner records file hashes and timing.
    pub fn to_json(&self) -> String {
        fn histogram(value: &BTreeMap<u64, u64>) -> String {
            format!(
                "{{{}}}",
                value
                    .iter()
                    .map(|(k, v)| format!("\"{k}\":{v}"))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        }
        let numeric = [
            ("input_reads", self.input_reads),
            ("output_molecules", self.output_molecules),
            ("removed_reads", self.removed_reads),
            ("exact_acgt_families", self.exact_acgt_families),
            ("multi_transcript_families", self.multi_transcript_families),
            (
                "mixed_insert_sequence_families",
                self.mixed_insert_sequence_families,
            ),
            ("umi_corrected_families", self.umi_corrected_families),
            ("acgt_fragment_contexts", self.acgt_fragment_contexts),
            ("ambiguous_umi_reads_kept", self.ambiguous_umi_reads_kept),
            (
                "exact_molecules_including_ambiguous",
                self.exact_molecules_including_ambiguous,
            ),
            (
                "additional_molecules_collapsed_by_correction",
                self.additional_molecules_collapsed_by_correction,
            ),
            (
                "max_exact_umis_in_one_context",
                self.max_exact_umis_in_one_context,
            ),
            ("observed_umi_strings", self.observed_umi_strings),
        ]
        .iter()
        .map(|(key, value)| format!("  \"{key}\": {value}"))
        .collect::<Vec<_>>()
        .join(",\n");
        let base_counts = self
            .umi_base_counts_by_position
            .iter()
            .map(|counts| {
                format!(
                    "{{{}}}",
                    counts
                        .iter()
                        .map(|(base, count)| format!("\"{base}\":{count}"))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let top = self
            .top_umis
            .iter()
            .map(|(umi, reads)| format!("{{\"umi\":\"{umi}\",\"reads\":{reads}}}"))
            .collect::<Vec<_>>()
            .join(",");
        // sample and UMI strings have already passed a restrictive ASCII check.
        format!(
            concat!(
                "{{\n  \"sample\": \"{}\",\n  \"method\": \"unique\",\n  \"umi_length\": {},\n",
                "{},\n  \"family_size_histogram\": {},\n  \"observed_umis_per_family_histogram\": {},\n",
                "  \"umi_base_counts_by_position\": [{}],\n  \"top_umis\": [{}],\n",
                "  \"key\": \"library + chromosome + strand + complete R1 genomic exon path + UMI\",\n",
                "  \"ambiguous_umi_policy\": \"retain each N-containing read separately\",\n",
                "  \"correction_distance\": 0\n}}\n"
            ),
            self.sample,
            self.umi_length,
            numeric,
            histogram(&self.family_size_histogram),
            histogram(&self.observed_umis_per_family_histogram),
            base_counts,
            top
        )
    }
}

/// The caller opens a single physical library and owns compression/publication.
/// All validation completes before either output stream receives result records.
pub fn run<R: BufRead, W: Write, F: Write>(
    mut input: R,
    annotation: &ExonIndex,
    mut output: W,
    mut families: F,
    sample: &str,
    umi_length: usize,
) -> Result<UmiSummary, String> {
    if umi_length == 0
        || sample.is_empty()
        || !sample
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-'))
    {
        return Err("positive UMI length and a simple, nonempty sample identifier required".into());
    }
    let mut summary = UmiSummary {
        sample: sample.into(),
        umi_length,
        umi_base_counts_by_position: vec![BTreeMap::new(); umi_length],
        ..Default::default()
    };
    let mut buckets: HashMap<Context, Bucket> = HashMap::new();
    let mut ambiguous: Vec<Ambiguous> = Vec::new();
    let mut seen = HashSet::new();
    // Counter.most_common in the reference breaks frequency ties by first sighting.
    let mut umi_counts: HashMap<String, (u64, u64)> = HashMap::new();
    let mut line = String::new();
    let mut number = 0;
    loop {
        line.clear();
        if input.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
            break;
        }
        number += 1;
        // Python's text reader normalizes CRLF before storing the original row.
        let row = line.trim_end_matches(['\r', '\n']);
        let parsed = (|| -> Result<(), String> {
            let fields: Vec<_> = row.split('\t').collect();
            if fields.len() != 5 {
                return Err("expected five fields".into());
            }
            let (read, reference, start, cigar, sequence) =
                (fields[0], fields[1], fields[2], fields[3], fields[4]);
            if !seen.insert(read.to_owned()) {
                return Err(format!(
                    "read ID repeated: {read}; input must be almostuniquematch output"
                ));
            }
            let umi = read.rsplit_once('_').ok_or("read ID has no _UMI suffix")?.1;
            if umi.len() != umi_length
                || !umi
                    .bytes()
                    .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
            {
                return Err(format!("invalid {umi_length}-nt UMI {umi:?}"));
            }
            let tx = reference
                .split(|c: char| c == '|' || c.is_whitespace())
                .next()
                .filter(|s| !s.is_empty())
                .ok_or("empty transcript ID")?;
            let model = annotation
                .transcripts
                .get(tx)
                .ok_or_else(|| format!("transcript missing from annotation: {tx}"))?;
            let (mut reference_length, mut query_length, mut errors) = (0u64, 0u64, 0u64);
            for (length, op) in parse_cigar(cigar)? {
                if op != 'I' {
                    reference_length = reference_length
                        .checked_add(length)
                        .ok_or("CIGAR reference length overflow")?;
                }
                if op != 'D' {
                    query_length = query_length
                        .checked_add(length)
                        .ok_or("CIGAR query length overflow")?;
                }
                if op != '=' {
                    errors = errors
                        .checked_add(length)
                        .ok_or("CIGAR error count overflow")?;
                }
            }
            if reference_length < 2 || sequence.len() < 3 || sequence.len() as u64 != query_length {
                return Err("CIGAR/sequence length mismatch or fragment too short".into());
            }
            let start = start
                .parse::<u64>()
                .map_err(|e| format!("invalid start: {e}"))?;
            let blocks = model.blocks(start, reference_length)?;
            let context = Context {
                chrom: model.chrom.clone(),
                strand: model.strand,
                blocks,
            };
            let site = model.map_position(start + (reference_length - 2))?;
            let counter = umi_counts
                .entry(umi.to_owned())
                .or_insert((0, summary.input_reads));
            counter.0 += 1;
            for (i, base) in umi.chars().enumerate() {
                *summary.umi_base_counts_by_position[i]
                    .entry(base)
                    .or_default() += 1;
            }
            summary.input_reads += 1;
            if umi.contains('N') {
                ambiguous.push(Ambiguous {
                    context,
                    site,
                    umi: umi.into(),
                    row: row.into(),
                });
                return Ok(());
            }
            let bucket = buckets.entry(context).or_insert_with(|| Bucket {
                site,
                families: BTreeMap::new(),
            });
            if bucket.site != site {
                return Err("same fragment path maps to conflicting stop sites".into());
            }
            if let Some(family) = bucket.families.get_mut(umi) {
                family.count += 1;
                family.multiple_transcripts |= tx != family.first_transcript;
                family.multiple_sequences |= sequence != family.first_sequence;
                if (errors, read, row)
                    < (
                        family.errors,
                        family.representative.as_str(),
                        family.row.as_str(),
                    )
                {
                    family.errors = errors;
                    family.representative = read.into();
                    family.row = row.into();
                }
            } else {
                bucket.families.insert(
                    umi.into(),
                    Family {
                        count: 1,
                        errors,
                        representative: read.into(),
                        row: row.into(),
                        first_transcript: tx.into(),
                        first_sequence: sequence.into(),
                        multiple_transcripts: false,
                        multiple_sequences: false,
                    },
                );
            }
            Ok(())
        })();
        parsed.map_err(|e| format!("input line {number}: {e}"))?;
    }
    // These are no longer needed while writing; reduce live memory for sorting.
    drop(seen);
    summary.observed_umi_strings = umi_counts.len() as u64;
    let mut ordered_umis: Vec<_> = umi_counts.into_iter().collect();
    ordered_umis.sort_unstable_by(|a, b| b.1.0.cmp(&a.1.0).then(a.1.1.cmp(&b.1.1)));
    summary.top_umis = ordered_umis
        .into_iter()
        .take(20)
        .map(|(umi, (count, _))| (umi, count))
        .collect();
    summary.acgt_fragment_contexts = buckets.len() as u64;
    summary.ambiguous_umi_reads_kept = ambiguous.len() as u64;
    writeln!(families, "FamilyID\tSample\tChromosome\tStrand\tGenomicBlocks\tCandidatePosition\tRootUMI\tReads\tObservedUMIs\tRepresentativeRead\tStatus").map_err(|e| e.to_string())?;
    let mut ordered_buckets: Vec<_> = buckets.into_iter().collect();
    ordered_buckets.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    for (context, bucket) in ordered_buckets {
        summary.max_exact_umis_in_one_context = summary
            .max_exact_umis_in_one_context
            .max(bucket.families.len() as u64);
        summary.exact_acgt_families += bucket.families.len() as u64;
        for (umi, family) in bucket.families {
            summary.multi_transcript_families += u64::from(family.multiple_transcripts);
            summary.mixed_insert_sequence_families += u64::from(family.multiple_sequences);
            emit(
                &mut output,
                &mut families,
                &mut summary,
                &context,
                bucket.site,
                &umi,
                &family.row,
                family.count,
                "acgt",
            )?;
        }
    }
    ambiguous.sort_unstable_by(|a, b| (&a.context, &a.row).cmp(&(&b.context, &b.row)));
    for item in ambiguous {
        emit(
            &mut output,
            &mut families,
            &mut summary,
            &item.context,
            item.site,
            &item.umi,
            &item.row,
            1,
            "ambiguous_umi_kept",
        )?;
    }
    summary.exact_molecules_including_ambiguous =
        summary.exact_acgt_families + summary.ambiguous_umi_reads_kept;
    summary.removed_reads = summary.input_reads - summary.output_molecules;
    output.flush().map_err(|e| e.to_string())?;
    families.flush().map_err(|e| e.to_string())?;
    Ok(summary)
}

#[allow(clippy::too_many_arguments)]
fn emit<W: Write, F: Write>(
    out: &mut W,
    report: &mut F,
    summary: &mut UmiSummary,
    context: &Context,
    site: u64,
    umi: &str,
    row: &str,
    count: u64,
    status: &str,
) -> Result<(), String> {
    summary.output_molecules += 1;
    *summary.family_size_histogram.entry(count).or_default() += 1;
    *summary
        .observed_umis_per_family_histogram
        .entry(1)
        .or_default() += 1;
    let blocks = context
        .blocks
        .iter()
        .map(|(a, b)| format!("{a}-{b}"))
        .collect::<Vec<_>>()
        .join(",");
    let read = row.split('\t').next().unwrap();
    writeln!(out, "{row}").map_err(|e| e.to_string())?;
    writeln!(
        report,
        "{}:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        summary.sample,
        summary.output_molecules,
        summary.sample,
        context.chrom,
        context.strand,
        blocks,
        site,
        umi,
        count,
        umi,
        read,
        status
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Transcript;

    fn annotation() -> ExonIndex {
        let mut transcripts = HashMap::new();
        for (id, strand, exons) in [
            ("AT1G00010.1", '+', vec![(101, 110), (201, 220)]),
            ("AT1G00010.2", '+', vec![(51, 60), (101, 110), (201, 220)]),
            ("AT1G00010.3", '+', vec![(101, 112), (203, 220)]),
            ("AT1G00020.1", '-', vec![(101, 110), (201, 220)]),
        ] {
            transcripts.insert(
                id.into(),
                Transcript::new("1".into(), strand, "AT1".into(), exons).unwrap(),
            );
        }
        ExonIndex { transcripts }
    }

    fn row(name: &str, umi: &str, tx: &str, start: u64, cigar: &str, seq: &str) -> String {
        format!("{name}_{umi}\t{tx}\t{start}\t{cigar}\t{seq}\n")
    }

    fn basic(name: &str, umi: &str) -> String {
        row(name, umi, "AT1G00010.1", 8, "8=", "ACGTACGT")
    }

    fn execute(rows: &[String]) -> Result<(String, String, UmiSummary), String> {
        let mut out = Vec::new();
        let mut families = Vec::new();
        let summary = run(
            rows.concat().as_bytes(),
            &annotation(),
            &mut out,
            &mut families,
            "sample",
            8,
        )?;
        Ok((
            String::from_utf8(out).unwrap(),
            String::from_utf8(families).unwrap(),
            summary,
        ))
    }

    #[test]
    fn equivalent_isoforms_share_a_family_and_representative() {
        let (out, families, stats) = execute(&[
            basic("z", "AAAAAAAA"),
            row("a", "AAAAAAAA", "AT1G00010.2", 18, "8=", "ACGTACGT"),
            basic("b", "CCCCCCCC"),
        ])
        .unwrap();
        assert_eq!(stats.input_reads, 3);
        assert_eq!(stats.output_molecules, 2);
        assert_eq!(stats.multi_transcript_families, 1);
        assert!(out.starts_with("a_AAAAAAAA\t"));
        assert!(!out.contains("z_AAAAAAAA"));
        assert!(families.contains(
            "sample:1\tsample\t1\t+\t108-110,201-205\t204\tAAAAAAAA\t2\tAAAAAAAA\ta_AAAAAAAA\tacgt"
        ));
    }

    #[test]
    fn endpoints_strands_and_complete_paths_remain_distinct() {
        let (_, _, stats) = execute(&[
            basic("original", "AAAAAAAA"),
            row("start", "AAAAAAAA", "AT1G00010.1", 7, "9=", "AACGTACGT"),
            row("end", "AAAAAAAA", "AT1G00010.1", 8, "9=", "ACGTACGTA"),
            row("strand", "AAAAAAAA", "AT1G00020.1", 16, "8=", "ACGTACGT"),
            row("path", "AAAAAAAA", "AT1G00010.3", 8, "8=", "ACGTACGT"),
        ])
        .unwrap();
        assert_eq!(stats.output_molecules, 5);
    }

    #[test]
    fn unknown_umis_are_independent_and_output_last() {
        let (out, families, stats) = execute(&[
            basic("z", "AAAANAAA"),
            basic("a", "AAAANAAA"),
            basic("x", "AAAAAAAA"),
        ])
        .unwrap();
        assert_eq!(stats.output_molecules, 3);
        assert_eq!(stats.ambiguous_umi_reads_kept, 2);
        assert!(out.starts_with("x_AAAAAAAA"));
        assert!(out.find("a_AAAANAAA").unwrap() < out.find("z_AAAANAAA").unwrap());
        assert_eq!(families.matches("ambiguous_umi_kept").count(), 2);
    }

    #[test]
    fn fewer_error_bases_win_and_indels_use_reference_span() {
        let (out, _, stats) = execute(&[
            row("a", "AAAAAAAA", "AT1G00010.1", 8, "7=1X", "ACGTACGA"),
            row("b", "AAAAAAAA", "AT1G00010.1", 8, "4=1I4=", "ACGTAACGT"),
            basic("z", "AAAAAAAA"),
        ])
        .unwrap();
        assert_eq!(stats.output_molecules, 1);
        assert_eq!(stats.mixed_insert_sequence_families, 1);
        assert!(out.starts_with("z_AAAAAAAA"));
        assert_eq!(stats.family_size_histogram[&3], 1);
    }

    #[test]
    fn unique_output_is_input_order_independent() {
        let mut rows = vec![
            basic("a", "AAAAAAAA"),
            basic("b", "AAAAAAAA"),
            basic("x", "CCCCCCCC"),
            basic("n", "NAAAAAAA"),
        ];
        let first = execute(&rows).unwrap();
        rows.reverse();
        let second = execute(&rows).unwrap();
        assert_eq!(first.0, second.0);
        assert_eq!(first.1, second.1);
        assert_eq!(first.2.output_molecules, second.2.output_molecules);
    }

    #[test]
    fn malformed_input_fails_before_writing_results() {
        for rows in [
            vec![basic("a", "AAA")],
            vec![basic("a", "AAAAAAAZ")],
            vec![basic("a", "AAAAAAAA"), basic("a", "AAAAAAAA")],
            vec![row("a", "AAAAAAAA", "AT9G99999.1", 8, "8=", "ACGTACGT")],
            vec![row("a", "AAAAAAAA", "AT1G00010.1", 28, "8=", "ACGTACGT")],
            vec![row("a", "AAAAAAAA", "AT1G00010.1", 8, "8M", "ACGTACGT")],
            vec![row("a", "AAAAAAAA", "AT1G00010.1", 8, "7=", "ACGTACGT")],
            vec!["not five fields\n".into()],
        ] {
            let mut output = Vec::new();
            let mut families = Vec::new();
            assert!(
                run(
                    rows.concat().as_bytes(),
                    &annotation(),
                    &mut output,
                    &mut families,
                    "sample",
                    8
                )
                .is_err()
            );
            assert!(output.is_empty());
            assert!(families.is_empty());
        }
    }

    #[test]
    fn empty_input_has_explicit_zero_metrics_and_only_a_family_header() {
        let (out, families, stats) = execute(&[]).unwrap();
        assert!(out.is_empty());
        assert_eq!(families.lines().count(), 1);
        assert_eq!(stats.output_molecules, 0);
        assert_eq!(stats.umi_corrected_families, 0);
        assert!(stats.to_json().contains("\"output_molecules\": 0"));
        assert!(stats.top_umis.is_empty());
    }

    #[test]
    fn crlf_input_has_the_same_normalized_result() {
        let row = basic("a", "AAAAAAAA");
        let unix = execute(std::slice::from_ref(&row)).unwrap();
        let windows = execute(&[row.replace('\n', "\r\n")]).unwrap();
        assert_eq!(unix.0, windows.0);
        assert_eq!(unix.1, windows.1);
    }

    #[test]
    fn business_totals_conserve_reads_and_counter_ties_follow_input() {
        let (_, _, stats) = execute(&[
            basic("a", "CCCCCCCC"),
            basic("b", "AAAAAAAA"),
            basic("c", "AAAAAAAA"),
            basic("d", "CCCCCCCC"),
            basic("n", "NAAAAAAA"),
        ])
        .unwrap();
        assert_eq!(stats.top_umis[0], ("CCCCCCCC".into(), 2));
        assert_eq!(
            stats
                .family_size_histogram
                .iter()
                .map(|(size, n)| size * n)
                .sum::<u64>(),
            5
        );
        assert_eq!(stats.umi_base_counts_by_position[0][&'A'], 2);
        assert_eq!(
            stats.exact_molecules_including_ambiguous,
            stats.output_molecules
        );
        assert_eq!(stats.removed_reads, 2);
    }
}
