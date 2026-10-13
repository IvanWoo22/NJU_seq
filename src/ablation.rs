//! Factorial reconstruction from disjoint, independently mapped read cohorts.
use crate::{fastq, model::Alignment};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
pub fn assemble(
    baseline: &Path,
    bac: &Path,
    host: &Path,
    bac_host: &Path,
    out: &Path,
) -> Result<(), String> {
    let mut excluded = HashSet::new();
    let mut input = fastq::Input::open(bac_host)?;
    let mut line = String::new();
    while input.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        excluded.insert(line.trim_end().to_owned());
        line.clear();
    }
    input.finish()?;
    fs::create_dir(out).map_err(|e| format!("Ablation output must be new: {e}"))?;
    let arms = ["both", "bac_only", "host_only", "neither"];
    let mut writers = arms
        .iter()
        .map(|arm| {
            File::create(out.join(format!("{arm}.tsv")))
                .map(BufWriter::new)
                .map_err(|e| e.to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut counts = [0u64; 4];
    let mut seen = HashSet::new();
    let mut cohorts = [0u64; 3];
    let mut overlap = 0;
    for (class, path) in [baseline, bac, host].iter().enumerate() {
        let mut input = fastq::Input::open(path)?;
        while input.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
            let a = Alignment::parse(line.trim_end())?;
            let id = fastq::fragment_id(&a.read);
            if !seen.insert(id.to_owned()) {
                return Err(format!(
                    "Ablation needs disjoint, one-alignment-per-read cohorts: {id}"
                ));
            }
            cohorts[class] += 1;
            let host_positive = excluded.contains(id);
            if class == 1 && host_positive {
                overlap += 1;
            }
            // Baseline S; bacterial-deleted B; subsequent host-deleted H.
            // host_only = S + (B minus host-positive B), NOT S+B.
            let include = match class {
                0 => [true, true, true, true],
                1 => [false, false, !host_positive, true],
                _ => [false, true, false, true],
            };
            for i in 0..4 {
                if include[i] {
                    writeln!(writers[i], "{}", line.trim_end()).map_err(|e| e.to_string())?;
                    counts[i] += 1;
                }
            }
            line.clear();
        }
        input.finish()?;
    }
    for mut w in writers {
        w.flush().map_err(|e| e.to_string())?;
    }
    let mut summary = String::from("arm\tmrna_unique_reads\n");
    for (a, n) in arms.iter().zip(counts) {
        summary.push_str(&format!("{a}\t{n}\n"));
    }
    fs::write(out.join("summary.tsv"), summary).map_err(|e| e.to_string())?;
    fs::write(out.join("cohorts.json"),format!("{{\"baseline\":{},\"rescued_bac\":{},\"rescued_host\":{},\"rescued_bac_host_positive\":{},\"host_positive_bac_reads\":{}}}\n",cohorts[0],cohorts[1],cohorts[2],overlap,excluded.len())).map_err(|e|e.to_string())?;
    fs::write(
        out.join("SUCCESS"),
        "Disjoint cohorts verified; baseline-first order, then B, then H.\n",
    )
    .map_err(|e| e.to_string())
}
