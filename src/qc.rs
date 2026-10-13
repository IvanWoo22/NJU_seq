//! Streaming short-fragment QC; counts records, never unique sequences.
use crate::{fastq, model::Alignment};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

#[derive(Default, Clone)]
pub struct Stats {
    pub reads: u64,
    pub bases: u64,
    pub lengths: BTreeMap<usize, u64>,
    pub base: [u64; 6],
    pub head: [u64; 6],
    pub tail: [u64; 6],
    pub quality_bases: u64,
    pub quality_sum: u64,
    pub q20: u64,
    pub q30: u64,
}
fn base_index(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        b'N' => 4,
        _ => 5,
    }
}
impl Stats {
    pub fn add(&mut self, seq: &str, quality: Option<&str>) -> Result<(), String> {
        if seq.is_empty() {
            return Err("QC cannot count an empty sequence".into());
        }
        self.reads += 1;
        self.bases += seq.len() as u64;
        *self.lengths.entry(seq.len()).or_default() += 1;
        for b in seq.bytes() {
            self.base[base_index(b)] += 1;
        }
        self.head[base_index(seq.as_bytes()[0])] += 1;
        self.tail[base_index(*seq.as_bytes().last().unwrap())] += 1;
        if let Some(q) = quality {
            if q.len() != seq.len() {
                return Err("QC quality/sequence length mismatch".into());
            }
            for q in q.bytes() {
                if !(33..=126).contains(&q) {
                    return Err("Invalid Phred+33 character".into());
                }
                let q = (q - 33) as u64;
                self.quality_bases += 1;
                self.quality_sum += q;
                self.q20 += u64::from(q >= 20);
                self.q30 += u64::from(q >= 30);
            }
        }
        Ok(())
    }
    pub fn quantile(&self, num: u64, den: u64) -> usize {
        let target = (self.reads * num).div_ceil(den);
        let mut total = 0;
        for (&length, &n) in &self.lengths {
            total += n;
            if total >= target {
                return length;
            }
        }
        0
    }
    pub fn below(&self, cut: usize) -> u64 {
        self.lengths.range(..cut).map(|(_, n)| n).sum()
    }
}
pub const SUMMARY_HEADER: &str = "sample\tstage\tunit\tmate\treads\tbases\tmin_nt\tp10_nt\tmedian_nt\tp90_nt\tmax_nt\tmean_nt\tlt10\tlt18\tlt30\tgc_fraction\tn_fraction\tmean_phred\tq20_fraction\tq30_fraction";
/// Explicit sensitivity-test gate on accepted R1 sequences; never enabled by default.
pub fn length_select(input: &Path, output: &Path, minimum: usize) -> Result<(u64, u64), String> {
    let mut input = fastq::Input::open(input)?;
    let mut writer = BufWriter::new(File::create(output).map_err(|e| e.to_string())?);
    let mut line = String::new();
    let (mut total, mut kept) = (0, 0);
    while input.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        let a = Alignment::parse(line.trim_end())?;
        total += 1;
        if a.sequence.len() >= minimum {
            writeln!(writer, "{}", line.trim_end()).map_err(|e| e.to_string())?;
            kept += 1;
        }
        line.clear();
    }
    input.finish()?;
    writer.flush().map_err(|e| e.to_string())?;
    Ok((total, kept))
}
pub fn write_metrics(
    root: &Path,
    sample: &str,
    stage: &str,
    unit: &str,
    stats: &[Stats],
) -> Result<(), String> {
    fs::create_dir(root).map_err(|e| format!("QC output must be new: {e}"))?;
    let mut summary =
        BufWriter::new(File::create(root.join("summary.tsv")).map_err(|e| e.to_string())?);
    let mut lengths =
        BufWriter::new(File::create(root.join("lengths.tsv")).map_err(|e| e.to_string())?);
    let mut bases =
        BufWriter::new(File::create(root.join("bases.tsv")).map_err(|e| e.to_string())?);
    writeln!(summary, "{SUMMARY_HEADER}").map_err(|e| e.to_string())?;
    writeln!(lengths, "sample\tstage\tunit\tmate\tlength_nt\tcount").map_err(|e| e.to_string())?;
    writeln!(
        bases,
        "sample\tstage\tmate\tregion\tA\tC\tG\tT\tN\tother\ttotal"
    )
    .map_err(|e| e.to_string())?;
    let ratio = |n: u64, d: u64| {
        if d == 0 {
            "NA".into()
        } else {
            format!("{:.8}", n as f64 / d as f64)
        }
    };
    for (i, s) in stats.iter().enumerate() {
        let mate = format!("R{}", i + 1);
        writeln!(summary,"{sample}\t{stage}\t{unit}\t{mate}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.reads,s.bases,s.lengths.keys().next().copied().unwrap_or(0),s.quantile(1,10),s.quantile(1,2),s.quantile(9,10),s.lengths.keys().next_back().copied().unwrap_or(0),ratio(s.bases,s.reads),s.below(10),s.below(18),s.below(30),ratio(s.base[1]+s.base[2],s.bases),ratio(s.base[4],s.bases),ratio(s.quality_sum,s.quality_bases),ratio(s.q20,s.quality_bases),ratio(s.q30,s.quality_bases)).map_err(|e|e.to_string())?;
        for (len, n) in &s.lengths {
            writeln!(lengths, "{sample}\t{stage}\t{unit}\t{mate}\t{len}\t{n}")
                .map_err(|e| e.to_string())?;
        }
        for (region, counts) in [("body", s.base), ("head", s.head), ("tail", s.tail)] {
            writeln!(
                bases,
                "{sample}\t{stage}\t{mate}\t{region}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                counts[0],
                counts[1],
                counts[2],
                counts[3],
                counts[4],
                counts[5],
                counts.iter().sum::<u64>()
            )
            .map_err(|e| e.to_string())?;
        }
    }
    summary.flush().map_err(|e| e.to_string())?;
    lengths.flush().map_err(|e| e.to_string())?;
    bases.flush().map_err(|e| e.to_string())?;
    Ok(())
}
pub fn fastq_qc(
    r1: &Path,
    r2: Option<&Path>,
    root: &Path,
    sample: &str,
    stage: &str,
) -> Result<(), String> {
    let mut one = fastq::Input::open(r1)?;
    let mut two = r2.map(fastq::Input::open).transpose()?;
    let mut stats = vec![Stats::default(); if two.is_some() { 2 } else { 1 }];
    while let Some(a) = fastq::record(&mut one)? {
        stats[0].add(&a.sequence, Some(&a.quality))?;
        if let Some(two) = &mut two {
            let b = fastq::record(two)?.ok_or("R2 shorter than R1")?;
            if fastq::fragment_id(a.id()) != fastq::fragment_id(b.id()) {
                return Err("QC pair ID mismatch".into());
            }
            stats[1].add(&b.sequence, Some(&b.quality))?;
        }
    }
    one.finish()?;
    if let Some(mut two) = two {
        if fastq::record(&mut two)?.is_some() {
            return Err("R2 longer than R1".into());
        }
        two.finish()?;
    }
    write_metrics(
        root,
        sample,
        stage,
        if r2.is_some() { "pairs" } else { "reads" },
        &stats,
    )
}
pub fn alignment_qc(
    input: &Path,
    root: &Path,
    sample: &str,
    stage: &str,
    unit: &str,
) -> Result<(), String> {
    let mut input = fastq::Input::open(input)?;
    let mut line = String::new();
    let mut stats = Stats::default();
    let mut previous = String::new();
    while input.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        let a = Alignment::parse(line.trim_end())?;
        if a.read != previous {
            stats.add(&a.sequence, None)?;
            previous = a.read;
        }
        line.clear();
    }
    input.finish()?;
    write_metrics(root, sample, stage, unit, &[stats])
}
fn names(path: &Path) -> Result<HashSet<String>, String> {
    let mut f = fastq::Input::open(path)?;
    let mut names = HashSet::new();
    let mut line = String::new();
    while f.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        let id = line.trim_end().split('\t').next().unwrap();
        if !id.is_empty() {
            names.insert(fastq::fragment_id(id).into());
        }
        line.clear();
    }
    f.finish()?;
    Ok(names)
}
pub fn partition(
    r1: &Path,
    r2: &Path,
    bac: &Path,
    host: &Path,
    out: &Path,
    sample: &str,
) -> Result<(), String> {
    let mut bac = names(bac)?;
    let mut host = names(host)?;
    if bac.intersection(&host).next().is_some() {
        return Err("Sequential bac/host removed sets overlap".into());
    }
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let mut writers = vec![];
    for class in ["removed_bac", "removed_host"] {
        for mate in ["R1", "R2"] {
            writers.push(fastq::Output::open(
                &out.join(format!("{class}.{mate}.fq.gz")),
                true,
            )?);
        }
    }
    let mut one = fastq::Input::open(r1)?;
    let mut two = fastq::Input::open(r2)?;
    let mut stats = vec![vec![Stats::default(); 2]; 5];
    while let Some(a) = fastq::record(&mut one)? {
        let b = fastq::record(&mut two)?.ok_or("Truncated R2")?;
        let id = fastq::fragment_id(a.id());
        if id != fastq::fragment_id(b.id()) {
            return Err("Partition pair ID mismatch".into());
        }
        let category = if bac.remove(id) {
            3
        } else if host.remove(id) {
            4
        } else {
            2
        };
        for (i, r) in [&a, &b].iter().enumerate() {
            stats[0][i].add(&r.sequence, Some(&r.quality))?;
            if category != 3 {
                stats[1][i].add(&r.sequence, Some(&r.quality))?;
            }
            stats[category][i].add(&r.sequence, Some(&r.quality))?;
        }
        if category >= 3 {
            writers[(category - 3) * 2].record(&a)?;
            writers[(category - 3) * 2 + 1].record(&b)?;
        }
    }
    if fastq::record(&mut two)?.is_some() {
        return Err("R2 longer than R1".into());
    }
    one.finish()?;
    two.finish()?;
    for w in writers {
        w.finish()?;
    }
    if !bac.is_empty() || !host.is_empty() {
        return Err(format!(
            "Removal IDs absent from FASTQ: bac={} host={}",
            bac.len(),
            host.len()
        ));
    }
    for (stage, s) in [
        "prepared",
        "after_bac",
        "after_host",
        "removed_bac",
        "removed_host",
    ]
    .iter()
    .zip(stats)
    {
        write_metrics(&out.join(format!("qc_{stage}")), sample, stage, "pairs", &s)?;
    }
    fs::write(
        out.join("SUCCESS"),
        "Partition and removal-list membership verified\n",
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn short_length_counts_weight_records_not_unique_sequences() {
        let mut s = Stats::default();
        s.add("ACNN", Some("II!!")).unwrap();
        s.add("ACNN", Some("IIII")).unwrap();
        s.add("T", None).unwrap();
        assert_eq!(s.reads, 3);
        assert_eq!(s.bases, 9);
        assert_eq!(s.quantile(1, 2), 4);
        assert_eq!(s.base, [2, 2, 0, 1, 4, 0]);
        assert_eq!(s.q30, 6);
        assert_eq!(s.quality_bases, 8);
    }
}
