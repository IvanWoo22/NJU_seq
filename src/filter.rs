//! Streaming equivalent of the current SAM -> quality -> multimatch -> location
//! filtering stages. Input must keep every read's records contiguous, as required
//! by the existing multimatch scripts. This module deliberately does not change
//! their selection rules and never performs PCR/UMI deduplication.

use crate::model::{Alignment, ExonIndex, parse_cigar};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InputFormat {
    Sam,
    /// Proof-only exact forward R1; the consumer also verifies paired sequences.
    ExactSingleEndSam,
    AlignmentTsv,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strategy {
    Mrna,
    Bac,
    Rrna,
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct FilterStats {
    pub input_records: u64,
    pub sam_selected: u64,
    pub quality_passed: u64,
    pub reference_unique: u64,
    pub coordinate_unique: u64,
    pub read_groups: u64,
}

#[derive(Debug, Default)]
struct CigarSummary {
    reference_len: u64,
    matches: u64,
    differences: u64,
    insertions: u64,
}

fn checked_add(target: &mut u64, amount: u64) -> Result<(), String> {
    *target = target.checked_add(amount).ok_or("CIGAR length overflow")?;
    Ok(())
}

fn cigar_summary(cigar: &str) -> Result<CigarSummary, String> {
    let mut summary = CigarSummary::default();
    for (length, operation) in parse_cigar(cigar)? {
        if length == 0 {
            return Err(format!("Zero-length CIGAR operation: {cigar}"));
        }
        match operation {
            '=' => {
                checked_add(&mut summary.reference_len, length)?;
                checked_add(&mut summary.matches, length)?;
            }
            'X' | 'D' => {
                checked_add(&mut summary.reference_len, length)?;
                checked_add(&mut summary.differences, length)?;
            }
            'I' => {
                checked_add(&mut summary.insertions, length)?;
                checked_add(&mut summary.differences, length)?;
            }
            _ => {
                return Err(format!(
                    "Expected end-to-end --xeq CIGAR (= X I D), got: {cigar}"
                ));
            }
        }
    }
    Ok(summary)
}

pub(crate) fn quality_passes(cigar: &str, strategy: Strategy) -> Result<bool, String> {
    if strategy == Strategy::Bac {
        // The bacterial stage does not run matchquality_judge.pl.
        return Ok(true);
    }
    let summary = cigar_summary(cigar)?;
    let first_boundary = if strategy == Strategy::Mrna { 18 } else { 20 };
    let permitted = if summary.reference_len < first_boundary {
        0
    } else if summary.reference_len < 30 {
        1
    } else if summary.reference_len < 40 {
        2
    } else {
        3
    };
    Ok(summary.differences <= permitted)
}

fn parse_sam(line: &str) -> Result<Option<Alignment>, String> {
    if line.starts_with('@') {
        return Ok(None);
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return Err("SAM record requires at least 11 fields".into());
    }
    let flag = fields[1]
        .parse::<u16>()
        .map_err(|_| format!("Invalid SAM FLAG: {}", fields[1]))?;
    // Equivalent to samtools view -F 128 | samtools view -F 16. Secondary
    // alignments are retained; the existing pipeline uses Bowtie2 -a or -k.
    if flag & 144 != 0 || fields[5] == "*" || fields[6] != "=" {
        return Ok(None);
    }
    let start = fields[3]
        .parse::<u64>()
        .map_err(|_| format!("Invalid SAM POS: {}", fields[3]))?;
    let mate_start = fields[7]
        .parse::<u64>()
        .map_err(|_| format!("Invalid SAM PNEXT: {}", fields[7]))?;
    if start != mate_start {
        return Ok(None);
    }
    Alignment::parse(&format!(
        "{}\t{}\t{}\t{}\t{}",
        fields[0], fields[2], start, fields[5], fields[9]
    ))
    .map(Some)
}

fn parse_exact_single_end_sam(line: &str) -> Result<Option<Alignment>, String> {
    let f: Vec<_> = line.split('\t').collect();
    if f.len() < 11 {
        return Err("SAM record requires at least 11 fields".into());
    }
    let flag = f[1].parse::<u16>().map_err(|e| e.to_string())?;
    if flag & (4 | 16 | 128) != 0 || f[5] == "*" {
        return Ok(None);
    }
    let a = Alignment::parse(&format!("{}\t{}\t{}\t{}\t{}", f[0], f[2], f[3], f[5], f[9]))?;
    if parse_cigar(&a.cigar)?.iter().any(|(_, op)| *op != '=') {
        return Ok(None);
    }
    Ok(Some(a))
}

/// The original dedup.pl resolves the leading `\w+\.\d+` transcript ID;
/// reference descriptions following that identifier do not affect the lookup.
fn legacy_transcript_id(reference: &str) -> Result<&str, String> {
    if reference.starts_with("ENST") {
        let id = reference
            .split(|c: char| c == '|' || c.is_whitespace())
            .next()
            .unwrap();
        crate::model::human_stable_id(id, "ENST")?;
        return Ok(id);
    }
    let bytes = reference.as_bytes();
    let mut offset = 0;
    while offset < bytes.len() && (bytes[offset].is_ascii_alphanumeric() || bytes[offset] == b'_') {
        offset += 1;
    }
    if offset == 0 || bytes.get(offset) != Some(&b'.') {
        return Err(format!("Invalid transcript ID: {reference}"));
    }
    offset += 1;
    let digits_start = offset;
    while offset < bytes.len() && bytes[offset].is_ascii_digit() {
        offset += 1;
    }
    if digits_start == offset {
        return Err(format!("Invalid transcript ID: {reference}"));
    }
    Ok(&reference[..offset])
}

fn emit<W: Write>(alignment: &Alignment, writer: &mut W) -> Result<(), String> {
    writeln!(writer, "{}", alignment.to_tsv()).map_err(|error| error.to_string())
}

fn mrna_group<W: Write>(
    alignments: &[Alignment],
    writer: &mut W,
    index: Option<&ExonIndex>,
    stats: &mut FilterStats,
) -> Result<(), String> {
    let mut reference_seen: HashSet<(&str, u64)> = HashSet::new();
    let mut genomic_seen: HashSet<(&str, char, u64)> = HashSet::new();
    for alignment in alignments {
        if !reference_seen.insert((&alignment.reference, alignment.start)) {
            continue;
        }
        stats.reference_unique += 1;
        if let Some(index) = index {
            let transcript_id = legacy_transcript_id(&alignment.reference)?;
            let transcript = index
                .transcripts
                .get(transcript_id)
                .ok_or_else(|| format!("Transcript absent from annotation: {transcript_id}"))?;
            let genomic_start = transcript.map_position(alignment.start)?;
            if !genomic_seen.insert((transcript.chrom.as_str(), transcript.strand, genomic_start)) {
                continue;
            }
        }
        emit(alignment, writer)?;
        stats.coordinate_unique += 1;
    }
    Ok(())
}

fn rrna_group<W: Write>(
    alignments: &[Alignment],
    writer: &mut W,
    stats: &mut FilterStats,
) -> Result<(), String> {
    if alignments.len() == 1 {
        emit(&alignments[0], writer)?;
        stats.reference_unique += 1;
        stats.coordinate_unique += 1;
        return Ok(());
    }
    // The literal tab delimiter and lexical tie break exactly mirror the
    // corrected Perl implementation, including lexical (not numeric) sites.
    let mut candidates: HashMap<String, (usize, u64)> = HashMap::new();
    let mut perfect_keys = Vec::new();
    for (position, alignment) in alignments.iter().enumerate() {
        let summary = cigar_summary(&alignment.cigar)?;
        let key = format!("{}\t{}", alignment.reference, alignment.start);
        let replace = candidates
            .get(&key)
            .is_none_or(|(_, matches)| summary.matches > *matches);
        if replace {
            candidates.insert(key.clone(), (position, summary.matches));
            if summary.reference_len == summary.matches && summary.insertions == 0 {
                perfect_keys.push(key);
            }
        }
    }
    let mut selected = Vec::new();
    if perfect_keys.is_empty() {
        if let Some((_, (position, _))) =
            candidates
                .iter()
                .max_by(|(key_a, (_, match_a)), (key_b, (_, match_b))| {
                    match_a.cmp(match_b).then_with(|| key_a.cmp(key_b))
                })
        {
            selected.push(*position);
        }
    } else {
        let mut printed = HashSet::new();
        for key in perfect_keys {
            if printed.insert(key.clone()) {
                // Compatibility detail: Perl remembers that a key was perfect
                // even if a later, longer imperfect alignment replaced it.
                selected.push(candidates[&key].0);
            }
        }
    }
    for position in selected {
        emit(&alignments[position], writer)?;
        stats.reference_unique += 1;
        stats.coordinate_unique += 1;
    }
    Ok(())
}

fn flush_group<W: Write>(
    alignments: &[Alignment],
    writer: &mut W,
    strategy: Strategy,
    index: Option<&ExonIndex>,
    stats: &mut FilterStats,
) -> Result<(), String> {
    match strategy {
        Strategy::Mrna => mrna_group(alignments, writer, index, stats),
        Strategy::Bac | Strategy::Rrna => rrna_group(alignments, writer, stats),
    }
}

/// Filter SAM records or five-column alignment rows. `index` optionally enables
/// mRNA coordinate filtering after quality and transcript/start filtering.
///
/// The caller must supply contiguous read groups. Neither the legacy
/// multimatch scripts nor this streaming replacement accepts arbitrary row
/// order. Memory is proportional to the largest read group plus exon index.
pub fn run<R: BufRead, W: Write>(
    mut reader: R,
    writer: &mut W,
    format: InputFormat,
    strategy: Strategy,
    index: Option<&ExonIndex>,
) -> Result<FilterStats, String> {
    if strategy != Strategy::Mrna && index.is_some() {
        return Err("Coordinate deduplication is only defined for mRNA filtering".into());
    }
    let mut stats = FilterStats::default();
    let mut line = String::new();
    let mut line_number = 0_u64;
    let mut current_read: Option<String> = None;
    let mut group = Vec::new();
    loop {
        line.clear();
        if reader
            .read_line(&mut line)
            .map_err(|error| error.to_string())?
            == 0
        {
            break;
        }
        line_number += 1;
        let record = line.trim_end_matches(['\r', '\n']);
        if format != InputFormat::AlignmentTsv && record.starts_with('@') {
            continue;
        }
        stats.input_records += 1;
        let alignment = match format {
            InputFormat::Sam => parse_sam(record),
            InputFormat::ExactSingleEndSam => parse_exact_single_end_sam(record),
            InputFormat::AlignmentTsv => Alignment::parse(record).map(Some),
        }
        .map_err(|error| format!("Input line {line_number}: {error}"))?;
        let Some(alignment) = alignment else { continue };
        stats.sam_selected += 1;
        if current_read.as_deref() != Some(alignment.read.as_str()) {
            flush_group(&group, writer, strategy, index, &mut stats)?;
            group.clear();
            current_read = Some(alignment.read.clone());
            stats.read_groups += 1;
        }
        if quality_passes(&alignment.cigar, strategy)
            .map_err(|error| format!("Input line {line_number}: {error}"))?
        {
            stats.quality_passed += 1;
            group.push(alignment);
        }
    }
    flush_group(&group, writer, strategy, index, &mut stats)?;
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Transcript;
    use std::io::Cursor;

    fn filtered(input: &str, strategy: Strategy) -> String {
        let mut output = Vec::new();
        run(
            Cursor::new(input),
            &mut output,
            InputFormat::AlignmentTsv,
            strategy,
            None,
        )
        .unwrap();
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn quality_boundaries_and_indels() {
        assert!(!quality_passes("16=1X", Strategy::Mrna).unwrap());
        assert!(quality_passes("17=1X", Strategy::Mrna).unwrap());
        assert!(!quality_passes("17=1X", Strategy::Rrna).unwrap());
        assert!(quality_passes("19=1X", Strategy::Rrna).unwrap());
        assert!(quality_passes("28=1D1X", Strategy::Mrna).unwrap());
        assert!(!quality_passes("28=1D1I1X", Strategy::Mrna).unwrap());
        assert!(quality_passes("37=3X", Strategy::Mrna).unwrap());
        assert!(!quality_passes("36=4X", Strategy::Mrna).unwrap());
        assert!(quality_passes("20M", Strategy::Mrna).is_err());
        assert!(quality_passes("0=", Strategy::Mrna).is_err());
    }

    #[test]
    fn mrna_keeps_first_and_does_not_conflate_keys_or_qnames() {
        let input = "q.1\ttx1\t23\t3=\tAAA\nq.1\ttx1\t23\t4=\tAAAA\nq.1\ttx12\t3\t3=\tCCC\nqX1\ttx1\t23\t3=\tGGG\n";
        assert_eq!(
            filtered(input, Strategy::Mrna),
            "q.1\ttx1\t23\t3=\tAAA\nq.1\ttx12\t3\t3=\tCCC\nqX1\ttx1\t23\t3=\tGGG\n"
        );
    }

    #[test]
    fn genomic_start_dedup_retains_strands_and_first_isoform() {
        let index = ExonIndex {
            transcripts: HashMap::from([
                (
                    "a.1".into(),
                    Transcript::new("1".into(), '+', "a".into(), vec![(10, 12)]).unwrap(),
                ),
                (
                    "b.1".into(),
                    Transcript::new("1".into(), '+', "b".into(), vec![(10, 12)]).unwrap(),
                ),
                (
                    "c.1".into(),
                    Transcript::new("1".into(), '-', "c".into(), vec![(8, 10)]).unwrap(),
                ),
            ]),
        };
        let input = "q\ta.1|description\t1\t3=\tAAA\nq\tb.1\t1\t3=\tAAA\nq\tc.1\t1\t3=\tTTT\n";
        let mut output = Vec::new();
        let stats = run(
            Cursor::new(input),
            &mut output,
            InputFormat::AlignmentTsv,
            Strategy::Mrna,
            Some(&index),
        )
        .unwrap();
        assert_eq!(stats.reference_unique, 3);
        assert_eq!(stats.coordinate_unique, 2);
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "q\ta.1|description\t1\t3=\tAAA\nq\tc.1\t1\t3=\tTTT\n"
        );
    }

    #[test]
    fn missing_and_out_of_range_transcripts_fail() {
        let index = ExonIndex {
            transcripts: HashMap::from([(
                "a.1".into(),
                Transcript::new("1".into(), '+', "a".into(), vec![(10, 12)]).unwrap(),
            )]),
        };
        for input in ["q\tz.1\t1\t3=\tAAA\n", "q\ta.1\t4\t3=\tAAA\n"] {
            assert!(
                run(
                    Cursor::new(input),
                    &mut Vec::new(),
                    InputFormat::AlignmentTsv,
                    Strategy::Mrna,
                    Some(&index)
                )
                .is_err()
            );
        }
    }

    #[test]
    fn human_par_y_does_not_collapse_onto_x() {
        let exons = concat!(
            "chrX\t101\t120\t+\tParent=ENST00001.1;gene_id=ENSG00001.1\n",
            "chrY\t101\t120\t+\tParent=ENST00001.1_PAR_Y;gene_id=ENSG00001.1\n",
        );
        let index = ExonIndex::from_reader(exons.as_bytes(), "auto", "").unwrap();
        let input = concat!(
            "q\tENST00001.1|ENSG00001.1|description\t1\t3=\tAAA\n",
            "q\tENST00001.1_PAR_Y|ENSG00001.1_PAR_Y|description\t1\t3=\tAAA\n",
        );
        let mut out = Vec::new();
        let stats = run(
            Cursor::new(input),
            &mut out,
            InputFormat::AlignmentTsv,
            Strategy::Mrna,
            Some(&index),
        )
        .unwrap();
        assert_eq!(stats.coordinate_unique, 2);
        assert_eq!(String::from_utf8(out).unwrap(), input);
    }

    #[test]
    fn chromosome_prefix_is_removed_only_once_by_annotation_loader() {
        let annotation = "chrchr1\t10\t12\t+\tParent=transcript:AT1.1\nchr1\t10\t12\t+\tParent=transcript:AT2.1\n";
        let index =
            ExonIndex::from_reader(annotation.as_bytes(), "Parent=transcript:", "AT").unwrap();
        let input = "q\tAT1.1\t1\t3=\tAAA\nq\tAT2.1\t1\t3=\tAAA\n";
        let mut output = Vec::new();
        let stats = run(
            Cursor::new(input),
            &mut output,
            InputFormat::AlignmentTsv,
            Strategy::Mrna,
            Some(&index),
        )
        .unwrap();
        assert_eq!(stats.coordinate_unique, 2);
        assert_eq!(String::from_utf8(output).unwrap(), input);
    }

    #[test]
    fn sam_selection_keeps_forward_r1_secondary() {
        let input = "@HD\tVN:1.6\nq\t99\ta.1\t1\t255\t3=\t=\t1\t3\tAAA\tIII\nq\t355\tb.1\t2\t255\t3=\t=\t2\t3\tAAA\tIII\nq\t147\ta.1\t1\t255\t3=\t=\t1\t-3\tTTT\tIII\nr\t83\ta.1\t1\t255\t3=\t=\t1\t3\tAAA\tIII\ns\t99\ta.1\t1\t255\t3=\t=\t2\t3\tAAA\tIII\nt\t99\ta.1\t1\t255\t3=\tb.1\t1\t3\tAAA\tIII\n";
        let mut output = Vec::new();
        let stats = run(
            Cursor::new(input),
            &mut output,
            InputFormat::Sam,
            Strategy::Mrna,
            None,
        )
        .unwrap();
        assert_eq!(stats.input_records, 6);
        assert_eq!(stats.sam_selected, 2);
        assert_eq!(stats.read_groups, 1);
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "q\ta.1\t1\t3=\tAAA\nq\tb.1\t2\t3=\tAAA\n"
        );
    }

    #[test]
    fn rrna_keeps_all_perfect_keys_once() {
        let input = "q\ta\t1\t3=\tAAA\nq\ta\t1\t4=\tAAAA\nq\tb\t1\t2=1I\tAAA\nq\tc\t1\t3=\tCCC\n";
        assert_eq!(
            filtered(input, Strategy::Bac),
            "q\ta\t1\t4=\tAAAA\nq\tc\t1\t3=\tCCC\n"
        );
    }

    #[test]
    fn rrna_best_nonperfect_uses_lexical_key_tie_break() {
        let input = "q\ta\t10\t2=1X\tAAA\nq\ta\t2\t2=1X\tCCC\n";
        assert_eq!(filtered(input, Strategy::Bac), "q\ta\t2\t2=1X\tCCC\n");
    }

    #[test]
    fn rrna_retains_legacy_perfect_history_rule() {
        let input = "q\ta\t1\t2=\tAA\nq\ta\t1\t3=1X\tAAAT\nq\tb\t1\t3=\tCCC\n";
        assert_eq!(
            filtered(input, Strategy::Bac),
            "q\ta\t1\t3=1X\tAAAT\nq\tb\t1\t3=\tCCC\n"
        );
    }

    #[test]
    fn empty_input_and_fully_filtered_groups_are_valid() {
        assert_eq!(filtered("", Strategy::Mrna), "");
        assert_eq!(
            filtered("q\ta\t1\t2=1X\tAAA\nr\tb\t2\t3=\tCCC\n", Strategy::Mrna),
            "r\tb\t2\t3=\tCCC\n"
        );
    }
}
