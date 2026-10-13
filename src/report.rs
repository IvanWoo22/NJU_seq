//! Self-contained, offline QC report. Numeric tables remain available for reuse.
use crate::workflow::json;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
fn collect(root: &Path, found: &mut Vec<PathBuf>, logs: &mut Vec<PathBuf>) -> Result<(), String> {
    for entry in fs::read_dir(root).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        if entry.file_type().map_err(|e| e.to_string())?.is_dir() {
            collect(&entry.path(), found, logs)?;
        } else if entry.file_name() == "summary.tsv"
            && entry.path().with_file_name("lengths.tsv").is_file()
        {
            found.push(entry.path());
        } else if entry
            .file_name()
            .to_string_lossy()
            .ends_with(".bowtie2.log")
        {
            logs.push(entry.path());
        }
    }
    Ok(())
}
fn stage_count(dir: &Path, names: &[&str]) -> Option<u64> {
    for name in names {
        if let Ok(text) = fs::read_to_string(dir.join(name).join("summary.tsv")) {
            for row in text.lines().skip(1) {
                let f: Vec<_> = row.split('\t').collect();
                if f.len() > 4 && f[3] == "R1" {
                    return f[4].parse().ok();
                }
            }
        }
    }
    None
}
fn alignment_summary(text: &str) -> Result<(Option<u64>, u64, String), String> {
    let (mut total, mut mapped) = (None, 0u64);
    let (mut summaries, mut timers, mut seconds) = (0, 0, 0u64);
    let mut block: Option<(u64, u64, [bool; 2])> = None;
    let valid = |b: Option<(u64, u64, [bool; 2])>| match b {
        None => Ok(()),
        Some((n, m, seen)) if m <= n && (seen == [true, true] || n == 0) => Ok(()),
        _ => Err("Incomplete or inconsistent Bowtie2 paired summary".to_string()),
    };
    for line in text.lines() {
        let number = line
            .split_whitespace()
            .next()
            .and_then(|x| x.replace(',', "").parse::<u64>().ok());
        if line.contains(" reads; of these:") || line.trim() == "0 reads" {
            valid(block)?;
            let n = number.ok_or("Invalid Bowtie2 read count")?;
            total = Some(
                total
                    .unwrap_or(0u64)
                    .checked_add(n)
                    .ok_or("Read count overflow")?,
            );
            block = Some((n, 0, [false; 2]));
            summaries += 1;
        }
        if line.contains("aligned concordantly exactly 1 time")
            || line.contains("aligned concordantly >1 times")
        {
            let n = number.ok_or("Invalid Bowtie2 concordant count")?;
            let b = block
                .as_mut()
                .ok_or("Concordant count without input count")?;
            let column = usize::from(line.contains("aligned concordantly >1 times"));
            if b.2[column] {
                return Err("Duplicate Bowtie2 concordant count".into());
            }
            b.2[column] = true;
            b.1 = b.1.checked_add(n).ok_or("Concordant count overflow")?;
            mapped = mapped.checked_add(n).ok_or("Concordant count overflow")?;
        }
        if let Some(t) = line.strip_prefix("Multiseed full-index search: ") {
            let parts = t
                .trim()
                .split(':')
                .map(str::parse::<u64>)
                .collect::<Result<Vec<_>, _>>()
                .map_err(|_| "Invalid Bowtie2 search duration")?;
            if parts.len() != 3 || parts[1] >= 60 || parts[2] >= 60 {
                return Err("Invalid Bowtie2 search duration".into());
            }
            let duration = parts[0]
                .checked_mul(3600)
                .and_then(|n| n.checked_add(parts[1] * 60 + parts[2]))
                .ok_or("Search duration overflow")?;
            seconds = seconds
                .checked_add(duration)
                .ok_or("Search duration overflow")?;
            timers += 1;
        }
    }
    valid(block)?;
    let elapsed = if summaries > 0 && timers == summaries {
        format!(
            "{:02}:{:02}:{:02}",
            seconds / 3600,
            seconds / 60 % 60,
            seconds % 60
        )
    } else {
        "NA".into()
    };
    Ok((total, mapped, elapsed))
}
fn alignment_table(path: &Path) -> Result<String, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let (mut total, mapped, elapsed) =
        alignment_summary(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let dir = path.parent().unwrap();
    let sample = dir.file_name().unwrap().to_string_lossy();
    let filename = path.file_name().unwrap().to_string_lossy();
    let stage = filename.trim_end_matches(".bowtie2.log");
    let prepared = stage_count(dir, &["qc_prepared"]);
    let after_bac = stage_count(dir, &["qc_after_bac"]);
    let after_host = stage_count(dir, &["qc_after_host"]);
    // A chunked run with zero input never launches Bowtie2; use the audited input QC.
    let input = match stage {
        "bac" => prepared,
        "rrna" => after_bac.or(prepared),
        "mrna" => after_host.or(after_bac).or(prepared),
        _ => None,
    };
    if total.is_none() && input == Some(0) {
        total = Some(0);
    }
    let selected = match stage {
        "bac" => prepared.zip(after_bac).and_then(|(a, b)| a.checked_sub(b)),
        "rrna" => after_bac
            .or(prepared)
            .zip(after_host)
            .and_then(|(a, b)| a.checked_sub(b)),
        "mrna" => stage_count(dir, &["qc_mrna", "qc_mrna_unique"]),
        _ => None,
    };
    let rate = total
        .filter(|n| *n > 0)
        .map(|n| format!("{:.8}", mapped as f64 / n as f64))
        .unwrap_or("NA".into());
    Ok(format!(
        "{{\"sample\":{},\"stage\":{},\"input_pairs\":{},\"concordant_pairs\":{},\"rate\":{},\"selected\":{},\"meaning\":{},\"search_time\":{}}}",
        json(&sample),
        json(stage),
        json(&total.map(|x| x.to_string()).unwrap_or("NA".into())),
        json(&mapped.to_string()),
        json(&rate),
        json(&selected.map(|x| x.to_string()).unwrap_or("NA".into())),
        json(if stage == "mrna" {
            "独定位 reads"
        } else {
            "实际删除 pairs"
        }),
        json(&elapsed)
    ))
}

fn table(path: &Path, source: &str) -> Result<Vec<String>, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = text.lines();
    let keys: Vec<_> = lines.next().ok_or("Empty QC table")?.split('\t').collect();
    let mut result = vec![];
    for line in lines {
        let values: Vec<_> = line.split('\t').collect();
        if values.len() != keys.len() {
            return Err(format!("Invalid QC table {}", path.display()));
        }
        let mut pairs = vec![format!("\"source\":{}", json(source))];
        pairs.extend(
            keys.iter()
                .zip(values)
                .map(|(k, v)| format!("{}:{}", json(k), json(v))),
        );
        result.push(format!("{{{}}}", pairs.join(",")));
    }
    Ok(result)
}
pub fn run(root: &Path, output: &Path) -> Result<(), String> {
    let mut found = vec![];
    let mut logs = vec![];
    collect(root, &mut found, &mut logs)?;
    found.sort();
    logs.sort();
    if found.is_empty() {
        return Err("No QC summary/length tables found".into());
    }
    let mut summaries = vec![];
    let mut lengths = vec![];
    let mut bases = vec![];
    for path in found {
        let folder = path.parent().unwrap();
        let source = folder.strip_prefix(root).unwrap().to_string_lossy();
        summaries.extend(table(&path, &source)?);
        lengths.extend(table(&folder.join("lengths.tsv"), &source)?);
        if folder.join("bases.tsv").is_file() {
            bases.extend(table(&folder.join("bases.tsv"), &source)?);
        }
    }
    let alignments = logs
        .iter()
        .map(|p| alignment_table(p))
        .collect::<Result<Vec<_>, _>>()?;
    let data = format!(
        "{{\"summary\":[{}],\"lengths\":[{}],\"bases\":[{}],\"alignment\":[{}]}}",
        summaries.join(","),
        lengths.join(","),
        bases.join(","),
        alignments.join(",")
    )
    .replace('<', "\\u003c");
    let html = include_str!("../assets/qc.html").replace("__DATA__", &data);
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .map_err(|e| e.to_string())?
        .write_all(html.as_bytes())
        .map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::alignment_summary;
    #[test]
    fn empty_chunk_uses_zero_input_qc() {
        let dir = std::env::temp_dir().join(format!(
            "nju-empty-qc-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log = dir.join("mrna.bowtie2.log");
        std::fs::write(&log, "").unwrap();
        for stage in ["qc_after_host", "qc_after_bac", "qc_prepared"] {
            let qc = dir.join(stage);
            std::fs::create_dir(&qc).unwrap();
            std::fs::write(
                qc.join("summary.tsv"),
                "sample\tstage\tunit\tmate\treads\nS\tinput\treads\tR1\t0\n",
            )
            .unwrap();
            let result = super::alignment_table(&log).unwrap();
            assert!(result.contains("\"input_pairs\":\"0\""), "{stage}");
            assert!(result.contains("\"rate\":\"NA\""));
            assert!(result.contains("\"search_time\":\"NA\""));
            std::fs::remove_dir_all(qc).unwrap();
        }
        std::fs::remove_dir_all(dir).unwrap();
    }
    fn batch(n: u64, unique: u64, multi: u64, time: &str) -> String {
        format!(
            "Multiseed full-index search: {time}\n{n} reads; of these:\n  {unique} (0%) aligned concordantly exactly 1 time\n  {multi} (0%) aligned concordantly >1 times\n"
        )
    }
    #[test]
    fn sums_batches_and_durations_over_one_day() {
        let log = batch(20000, 100, 900, "23:59:59") + &batch(10, 2, 3, "00:00:02");
        assert_eq!(
            alignment_summary(&log).unwrap(),
            (Some(20010), 1005, "24:00:01".into())
        );
    }
    #[test]
    fn single_zero_and_unavailable_timing() {
        assert_eq!(
            alignment_summary(&batch(10, 2, 3, "01:02:03")).unwrap(),
            (Some(10), 5, "01:02:03".into())
        );
        assert_eq!(
            alignment_summary("0 reads; of these:\n").unwrap(),
            (Some(0), 0, "NA".into())
        );
        assert_eq!(
            alignment_summary(
                "Multiseed full-index search: 00:00:01\n0 reads\n0.00% overall alignment rate\n"
            )
            .unwrap(),
            (Some(0), 0, "00:00:01".into())
        );
        assert_eq!(alignment_summary("").unwrap(), (None, 0, "NA".into()));
        let log = batch(10, 2, 3, "01:02:03") + "0 reads; of these:\n";
        assert_eq!(alignment_summary(&log).unwrap().2, "NA");
    }
    #[test]
    fn rejects_incomplete_and_impossible_summaries() {
        for log in [
            batch(1, 2, 0, "00:00:01"),
            batch(1, 0, 0, "00:99:00"),
            "5 reads; of these:\n".into(),
            "5 aligned concordantly exactly 1 time\n".into(),
            batch(10, 2, 3, "00:00:01") + "2 aligned concordantly exactly 1 time\n",
        ] {
            assert!(alignment_summary(&log).is_err(), "{log}");
        }
    }
}
