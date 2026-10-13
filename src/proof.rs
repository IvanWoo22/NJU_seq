//! Conservative absence proof, not an alternative mapper. Query length <=16
//! can only pass the existing reference-span quality thresholds with all '='.
//! Reaching span 18 needs >=2 deletions but allows only one. Never apply to 17.
use crate::{fastq, model::Alignment};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
fn encoded(sequence: &str) -> Option<(usize, u64)> {
    if sequence.is_empty() || sequence.len() > 16 {
        return None;
    }
    let mut value = 0;
    for b in sequence.bytes() {
        value = (value << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some((sequence.len(), value))
}
fn records(
    path: &Path,
    tsv: bool,
    mut visit: impl FnMut(&str, &str) -> Result<(), String>,
) -> Result<(), String> {
    let mut input = fastq::Input::open(path)?;
    if tsv {
        let mut line = String::new();
        while input.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
            let a = Alignment::parse(line.trim_end())?;
            visit(&a.read, &a.sequence)?;
            line.clear();
        }
    } else {
        while let Some(r) = fastq::record(&mut input)? {
            visit(r.id(), &r.sequence)?;
        }
    }
    input.finish()
}
pub fn impossible_short(
    input: &Path,
    tsv: bool,
    reference: &Path,
    output: &Path,
) -> Result<String, String> {
    let mut absent: [HashSet<u64>; 17] = std::array::from_fn(|_| HashSet::new());
    let (mut total, mut eligible) = (0u64, 0u64);
    records(input, tsv, |_, s| {
        total += 1;
        if let Some((len, value)) = encoded(s) {
            eligible += 1;
            absent[len].insert(value);
        }
        Ok(())
    })?;
    let unique: usize = absent.iter().map(HashSet::len).sum();
    let lengths: Vec<_> = (1..17).filter(|i| !absent[*i].is_empty()).collect();
    let masks: Vec<_> = (0..17).map(|length| (1u64 << (2 * length)) - 1).collect();
    let file = BufReader::new(File::open(reference).map_err(|e| e.to_string())?);
    let (mut value, mut valid, mut transcripts, mut bases) = (0u64, 0usize, 0u64, 0u64);
    for line in file.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('>') {
            transcripts += 1;
            value = 0;
            valid = 0;
            continue;
        }
        if line.is_empty() {
            continue;
        }
        if transcripts == 0 {
            return Err("FASTA sequence before header".into());
        }
        for b in line.trim_end().bytes() {
            bases += 1;
            let v = match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    value = 0;
                    valid = 0;
                    continue;
                }
            };
            value = ((value << 2) | v) & masks[16];
            valid = (valid + 1).min(16);
            for &len in &lengths {
                if len <= valid && !absent[len].is_empty() {
                    absent[len].remove(&(value & masks[len]));
                }
            }
        }
    }
    if transcripts == 0 {
        return Err("Empty reference FASTA".into());
    }
    let mut out = BufWriter::new(File::create(output).map_err(|e| e.to_string())?);
    let mut excluded = 0u64;
    records(input, tsv, |id, s| {
        if let Some((len, value)) = encoded(s)
            && absent[len].contains(&value)
        {
            writeln!(out, "{}", fastq::fragment_id(id)).map_err(|e| e.to_string())?;
            excluded += 1;
        }
        Ok(())
    })?;
    out.flush().map_err(|e| e.to_string())?;
    Ok(format!(
        "{{\"input_reads\":{total},\"eligible_reads\":{eligible},\"unique_short_sequences\":{unique},\"proven_impossible_reads\":{excluded},\"reference_transcripts\":{transcripts},\"reference_bases\":{bases},\"maximum_query_length\":16}}\n"
    ))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn length_17_and_unknown_bases_never_receive_absence_proof() {
        assert!(encoded(&"A".repeat(16)).is_some());
        assert!(encoded(&"A".repeat(17)).is_none());
        assert!(encoded("ACNG").is_none());
        assert_eq!(encoded("ACGT"), Some((4, 27)));
    }
    #[test]
    fn absence_bound_tracks_the_actual_reference_span_quality_policy() {
        use crate::filter::{Strategy, quality_passes};
        for q in 1..=16 {
            for d in 1..=50 {
                assert!(!quality_passes(&format!("{q}={d}D"), Strategy::Mrna).unwrap());
            }
        }
        assert!(quality_passes("17=1D", Strategy::Mrna).unwrap());
    }
}
