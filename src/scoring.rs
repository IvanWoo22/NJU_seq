//! Shared NJU-seq score engine with explicit mRNA/rRNA policies and one output schema.
use crate::workflow::json;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
type Result<T> = std::result::Result<T, String>;

fn reader(path: &Path) -> Result<BufReader<File>> {
    File::open(path)
        .map(BufReader::new)
        .map_err(|e| format!("{}: {e}", path.display()))
}
fn integer(s: &str, label: &str, positive: bool) -> Result<u64> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Invalid {label}: {s}"));
    }
    let n = s.parse::<u64>().map_err(|e| e.to_string())?;
    if positive && n == 0 {
        return Err(format!("{label} must be positive"));
    }
    Ok(n)
}
fn add(a: u64, b: u64, label: &str) -> Result<u64> {
    a.checked_add(b).ok_or_else(|| format!("{label} overflow"))
}
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct Site {
    chrom: String,
    position: u64,
    strand: char,
}
#[derive(Debug)]
struct ScoreInput {
    counts: HashMap<Site, u64>,
    info: HashMap<Site, String>,
}

fn read_score_input(path: &Path) -> Result<ScoreInput> {
    let mut result = ScoreInput {
        counts: HashMap::new(),
        info: HashMap::new(),
    };
    for (line_no, line) in reader(path)?.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let row: Vec<_> = line.split('\t').collect();
        if row.len() != 9 || !matches!(row[2], "+" | "-") || row[0].is_empty() {
            return Err(format!(
                "{}:{}: invalid nine-column count row",
                path.display(),
                line_no + 1
            ));
        }
        let site = Site {
            chrom: row[0].into(),
            position: integer(row[1], "Position", true)?,
            strand: row[2].chars().next().unwrap(),
        };
        let ends = integer(row[8], "EndCount", true)?;
        if result.counts.insert(site.clone(), ends).is_some() {
            return Err(format!(
                "{}: duplicate genomic position: {}:{}:{}",
                path.display(),
                site.chrom,
                site.position,
                site.strand
            ));
        }
        result.info.insert(site, row[3..7].join("\t"));
    }
    Ok(result)
}

fn sample_name(path: &Path) -> String {
    path.parent()
        .and_then(Path::file_name)
        .map(|n| n.to_string_lossy().into_owned())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| ".".into())
}

fn sum_counts<'a>(mut values: impl Iterator<Item = &'a u64>, label: &str) -> Result<u64> {
    values.try_fold(0u64, |a, &b| add(a, b, label))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mode {
    Mrna,
    Rrna,
}
impl Mode {
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "mrna" => Ok(Self::Mrna),
            "rrna" => Ok(Self::Rrna),
            _ => Err("score --mode must be mrna|rrna".into()),
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::Mrna => "mrna",
            Self::Rrna => "rrna",
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Format {
    Unified,
    Legacy,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Rows {
    Candidates,
    All,
}
pub struct Options {
    pub mode: Mode,
    pub format: Format,
    pub rows: Rows,
    pub reference: Option<String>,
}
#[derive(Clone, Copy, Debug)]
struct Local {
    signal: u64,
    adjacent: u64,
}
#[derive(Clone, Copy, Debug)]
struct Normalization {
    nc: u64,
    tr: u64,
}
impl Normalization {
    fn factor(self) -> Option<f64> {
        (self.nc > 0 && self.tr > 0).then(|| self.nc as f64 / self.tr as f64)
    }
}
#[derive(Debug)]
struct PairScore {
    corrected: f64,
    score: f64,
}
fn pair_score(tr: Local, nc: Local, norm: Normalization, mode: Mode) -> Option<PairScore> {
    norm.factor()?;
    let nc_signal = match mode {
        Mode::Mrna => nc.signal.max(1),
        Mode::Rrna => nc.signal,
    };
    Some(PairScore {
        // Preserve Perl's operation order; using raw * (NC/TR) changes ceil at ties.
        corrected: (tr.signal as f64 * norm.nc as f64 / norm.tr as f64).ceil(),
        score: tr.signal as f64 / tr.adjacent.max(1) as f64
            - nc_signal as f64 / nc.adjacent.max(1) as f64,
    })
}
#[derive(Debug)]
struct ScoredRow {
    meta: Vec<String>,
    nc: Local,
    tr: Vec<Local>,
    pairs: Vec<PairScore>,
    status: &'static str,
    total_score: Option<f64>,
    corrected: Option<f64>,
    fold: Option<f64>,
    fc_pass: bool,
    score_pass: bool,
}
impl ScoredRow {
    fn candidate(&self) -> bool {
        self.status == "evaluated" && self.fc_pass && self.score_pass
    }
}
fn evaluate(
    meta: Vec<String>,
    nc: Local,
    tr: Vec<Local>,
    totals: &[Normalization],
    mode: Mode,
    boundary: bool,
) -> ScoredRow {
    let status = if totals.iter().any(|n| n.factor().is_none()) {
        "normalization_unavailable"
    } else if boundary {
        "boundary_excluded"
    } else {
        "evaluated"
    };
    let pairs: Vec<_> = if status == "evaluated" {
        tr.iter()
            .zip(totals)
            .map(|(&t, &n)| pair_score(t, nc, n, mode).unwrap())
            .collect()
    } else {
        vec![]
    };
    let score = (status == "evaluated").then(|| pairs.iter().map(|p| p.score).sum::<f64>());
    let corrected = (status == "evaluated").then(|| pairs.iter().map(|p| p.corrected).sum::<f64>());
    let fold = corrected.map(|n| n / nc.signal.max(1) as f64);
    let n = tr.len() as f64;
    let fc_pass = corrected.is_some_and(|v| v > 0.0 && v >= 3.0 * n * nc.signal as f64);
    let score_pass = score.is_some_and(|v| v >= 30.0 * n);
    ScoredRow {
        meta,
        nc,
        tr,
        pairs,
        status,
        total_score: score,
        corrected,
        fold,
        fc_pass,
        score_pass,
    }
}
fn adjacent_site(site: &Site) -> Result<Site> {
    let mut adjacent = site.clone();
    adjacent.position = if site.strand == '+' {
        site.position - 1
    } else {
        site.position
            .checked_add(1)
            .ok_or("Adjacent position overflow")?
    };
    Ok(adjacent)
}
fn local(input: &ScoreInput, site: &Site, adjacent: &Site) -> Local {
    Local {
        signal: input.counts.get(site).copied().unwrap_or(0),
        adjacent: input.counts.get(adjacent).copied().unwrap_or(0),
    }
}
fn mrna_rows(
    inputs: &[PathBuf],
    which: Rows,
) -> Result<(Vec<ScoredRow>, Vec<Normalization>, usize)> {
    let nc = read_score_input(&inputs[0])?;
    let tr = inputs[1..]
        .iter()
        .map(|p| read_score_input(p))
        .collect::<Result<Vec<_>>>()?;
    let totals = tr
        .iter()
        .map(|t| {
            Ok(Normalization {
                tr: sum_counts(t.counts.values(), "TR total")?,
                nc: sum_counts(
                    t.counts.keys().filter_map(|k| nc.counts.get(k)),
                    "NC overlap total",
                )?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut shared: HashSet<_> = tr[0].counts.keys().collect();
    for t in &tr[1..] {
        shared.retain(|k| t.counts.contains_key(*k));
    }
    let considered = shared.len();
    let mut rows = vec![];
    for site in shared {
        let adjacent = adjacent_site(site)?;
        let mut meta = vec![
            site.chrom.clone(),
            site.position.to_string(),
            site.strand.to_string(),
        ];
        meta.extend(tr[0].info[site].split('\t').map(String::from));
        let row = evaluate(
            meta,
            local(&nc, site, &adjacent),
            tr.iter().map(|t| local(t, site, &adjacent)).collect(),
            &totals,
            Mode::Mrna,
            false,
        );
        if which == Rows::All || row.candidate() {
            rows.push(row);
        }
    }
    Ok((rows, totals, considered))
}
struct ReferenceCounts {
    bases: Vec<String>,
    ends: Vec<u64>,
}
fn read_reference_counts(path: &Path) -> Result<ReferenceCounts> {
    let mut counts = ReferenceCounts {
        bases: vec![],
        ends: vec![],
    };
    for line in reader(path)?.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let f: Vec<_> = line.split('\t').collect();
        if f.len() != 4 || f[1].len() != 1 || !f[1].bytes().all(|b| b.is_ascii_alphabetic()) {
            return Err(format!(
                "{}: expected four-column rRNA position/base/start/end counts",
                path.display()
            ));
        }
        let position = integer(f[0], "rRNA position", true)?;
        if position != counts.ends.len() as u64 + 1 {
            return Err("rRNA positions must be contiguous and 1-based".into());
        }
        integer(f[2], "rRNA start count", false)?;
        counts.ends.push(integer(f[3], "rRNA end count", false)?);
        counts.bases.push(f[1].into());
    }
    if counts.ends.is_empty() {
        return Err("Empty rRNA reference count table".into());
    }
    Ok(counts)
}
fn rrna_rows(
    inputs: &[PathBuf],
    reference: &str,
    which: Rows,
) -> Result<(Vec<ScoredRow>, Vec<Normalization>, usize)> {
    let counts = inputs
        .iter()
        .map(|p| read_reference_counts(p))
        .collect::<Result<Vec<_>>>()?;
    let nc = &counts[0];
    if counts[1..].iter().any(|t| t.bases != nc.bases) {
        return Err("rRNA reference lengths/bases differ between NC and TR".into());
    }
    let nc_total = sum_counts(nc.ends.iter(), "NC total")?;
    let totals = counts[1..]
        .iter()
        .map(|t| {
            Ok(Normalization {
                nc: nc_total,
                tr: sum_counts(t.ends.iter(), "TR total")?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut rows = vec![];
    for i in 0..nc.ends.len() - 1 {
        let position = i + 1;
        let meta = vec![
            reference.into(),
            position.to_string(),
            "+".into(),
            if i == 0 {
                "N".into()
            } else {
                nc.bases[i - 1].clone()
            },
            nc.bases[i].clone(),
            nc.bases[i + 1].clone(),
            reference.into(),
        ];
        let values: Vec<_> = counts
            .iter()
            .map(|t| Local {
                signal: t.ends[i + 1],
                adjacent: t.ends[i],
            })
            .collect();
        let row = evaluate(
            meta,
            values[0],
            values[1..].to_vec(),
            &totals,
            Mode::Rrna,
            position <= 20 || position > nc.ends.len().saturating_sub(20),
        );
        if which == Rows::All || row.candidate() {
            rows.push(row);
        }
    }
    Ok((rows, totals, nc.ends.len() - 1))
}
fn value(v: Option<f64>) -> String {
    v.map(|n| n.to_string()).unwrap_or_else(|| "NA".into())
}
fn criterion(row: &ScoredRow, pass: bool) -> &'static str {
    if row.status != "evaluated" {
        "not_evaluated"
    } else if pass {
        "meet"
    } else {
        "unmeet"
    }
}
fn common_header(names: &[String]) -> Vec<String> {
    let mut h = [
        "Chromosome",
        "Position",
        "Direction",
        "LeftBase",
        "Base",
        "RightBase",
        "GeneID",
    ]
    .map(String::from)
    .to_vec();
    h.push(format!("{}_count(NC)", names[0]));
    for s in &names[1..] {
        h.extend([
            format!("{s}_count"),
            format!("{s}_count_corr"),
            format!("{s}_score"),
        ]);
    }
    h.extend(["TotalFoldChange", "AveFC", "TotalScore", "AveScore"].map(String::from));
    h
}
fn common_values(row: &ScoredRow) -> Vec<String> {
    let mut fields = row.meta.clone();
    fields.push(row.nc.signal.to_string());
    for (i, tr) in row.tr.iter().enumerate() {
        fields.extend([
            tr.signal.to_string(),
            value(row.pairs.get(i).map(|p| p.corrected)),
            value(row.pairs.get(i).map(|p| p.score)),
        ]);
    }
    let n = row.tr.len() as f64;
    fields.extend([
        value(row.fold),
        value(row.fold.map(|v| v / n)),
        value(row.total_score),
        value(row.total_score.map(|v| v / n)),
    ]);
    fields
}
fn unified_header(names: &[String]) -> Vec<String> {
    let mut h = common_header(names);
    h.extend(
        [
            "ScoreMode",
            "CoordinateSystem",
            "EvaluationStatus",
            "Candidate",
            "FC_criteria",
            "Score_criteria",
            "TotalCorrectedCount",
        ]
        .map(String::from),
    );
    for suffix in ["adjacent_count", "score_numerator", "score_denominator"] {
        h.push(format!("{}_{suffix}", names[0]));
    }
    for s in &names[1..] {
        for suffix in [
            "adjacent_count",
            "score_denominator",
            "NC_TOTAL",
            "TR_TOTAL",
            "normalization_factor",
        ] {
            h.push(format!("{s}_{suffix}"));
        }
    }
    h
}
fn unified_values(row: &ScoredRow, totals: &[Normalization], mode: Mode) -> Vec<String> {
    let mut fields = common_values(row);
    fields.extend([
        mode.name().into(),
        if mode == Mode::Mrna {
            "genome_1based"
        } else {
            "rrna_reference_1based"
        }
        .into(),
        row.status.into(),
        row.candidate().to_string(),
        criterion(row, row.fc_pass).into(),
        criterion(row, row.score_pass).into(),
        value(row.corrected),
        row.nc.adjacent.to_string(),
        match mode {
            Mode::Mrna => row.nc.signal.max(1),
            Mode::Rrna => row.nc.signal,
        }
        .to_string(),
        row.nc.adjacent.max(1).to_string(),
    ]);
    for (tr, norm) in row.tr.iter().zip(totals) {
        fields.extend([
            tr.adjacent.to_string(),
            tr.adjacent.max(1).to_string(),
            norm.nc.to_string(),
            norm.tr.to_string(),
            value(norm.factor()),
        ]);
    }
    fields
}
fn legacy_rrna_header(names: &[String]) -> Vec<String> {
    let mut h = vec![
        "#Position".into(),
        "Base".into(),
        format!("{}_count", names[0]),
    ];
    for s in &names[1..] {
        h.extend([format!("{s}_count"), format!("{s}_score")]);
    }
    h.extend(["FoldChange", "FC_criteria", "TotalScore", "Score_criteria"].map(String::from));
    h
}
fn legacy_rrna_values(row: &ScoredRow) -> Vec<String> {
    let mut f = vec![
        row.meta[1].clone(),
        row.meta[4].clone(),
        row.nc.adjacent.to_string(),
    ];
    for (i, tr) in row.tr.iter().enumerate() {
        f.extend([
            tr.adjacent.to_string(),
            format!("{:.2}", row.pairs.get(i).map(|p| p.score).unwrap_or(0.0)),
        ]);
    }
    f.extend([
        format!("{:.2}", row.fold.unwrap_or(0.0)),
        if row.fc_pass { "meet" } else { "unmeet" }.into(),
        format!("{:.2}", row.total_score.unwrap_or(0.0)),
        if row.score_pass { "meet" } else { "unmeet" }.into(),
    ]);
    f
}

pub fn run(
    inputs: &[PathBuf],
    output: &Path,
    summary: Option<&Path>,
    options: &Options,
) -> Result<()> {
    if inputs.len() < 2 {
        return Err("score requires one NC and at least one TR".into());
    }
    if options.format == Format::Legacy && options.mode == Mode::Mrna && options.rows == Rows::All {
        return Err(
            "mRNA --rows all requires --format unified so unevaluable rows are explicit".into(),
        );
    }
    let reference = match (&options.reference, options.mode) {
        (Some(s), Mode::Rrna) if !s.is_empty() && !s.chars().any(char::is_whitespace) => s.as_str(),
        (None, Mode::Mrna) => "",
        (Some(_), Mode::Mrna) => return Err("--reference-id is only valid for rRNA mode".into()),
        _ => return Err("rRNA score requires a nonempty --reference-id without whitespace".into()),
    };
    let names: Vec<_> = inputs.iter().map(|p| sample_name(p)).collect();
    if names
        .iter()
        .any(|s| s.chars().any(|c| c == '\t' || c == '\n' || c == '\r'))
    {
        return Err("Sample names must not contain table separators".into());
    }
    let (mut rows, totals, considered) = match options.mode {
        Mode::Mrna => mrna_rows(inputs, options.rows)?,
        Mode::Rrna => rrna_rows(inputs, reference, options.rows)?,
    };
    let header = match (options.format, options.mode) {
        (Format::Unified, _) => unified_header(&names),
        (Format::Legacy, Mode::Mrna) => common_header(&names),
        (Format::Legacy, Mode::Rrna) => legacy_rrna_header(&names),
    };
    // Historical legacy callers may use identical parent names; unified output must be unambiguous.
    if options.format == Format::Unified
        && header.iter().collect::<HashSet<_>>().len() != header.len()
    {
        return Err("Duplicate sample names/column names; place NC and each TR in distinct sample directories".into());
    }
    if options.mode == Mode::Mrna {
        rows.sort_by(|a, b| {
            b.total_score
                .unwrap_or(f64::NEG_INFINITY)
                .total_cmp(&a.total_score.unwrap_or(f64::NEG_INFINITY))
                .then_with(|| {
                    common_values(a)
                        .join("\t")
                        .cmp(&common_values(b).join("\t"))
                })
        });
    }
    let mut out = BufWriter::new(File::create(output).map_err(|e| e.to_string())?);
    writeln!(out, "{}", header.join("\t")).map_err(|e| e.to_string())?;
    for row in &rows {
        let fields = match (options.format, options.mode) {
            (Format::Unified, _) => unified_values(row, &totals, options.mode),
            (Format::Legacy, Mode::Mrna) => common_values(row),
            (Format::Legacy, Mode::Rrna) => legacy_rrna_values(row),
        };
        writeln!(out, "{}", fields.join("\t")).map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())?;
    if let Some(path) = summary {
        let normal = totals.iter().enumerate().map(|(i,n)| format!(
            "{{\"sample\":{},\"NC_TOTAL\":{},\"TR_TOTAL\":{},\"normalization_factor\":{},\"status\":{}}}",
            json(&names[i+1]), n.nc, n.tr, n.factor().map(|x| x.to_string()).unwrap_or("null".into()),
            json(if n.tr == 0 { "zero_tr_total" } else if n.nc == 0 { "zero_nc_total" } else { "available" }))).collect::<Vec<_>>().join(",");
        let inputs_json = inputs
            .iter()
            .map(|p| json(&p.to_string_lossy()))
            .collect::<Vec<_>>()
            .join(",");
        let valid = totals.iter().all(|t| t.factor().is_some());
        let data = format!(
            concat!(
                "{{\"schema_version\":1,\"engine_version\":{},\"mode\":{},\"format\":{},\"rows\":{},\"reference_id\":{},",
                "\"normalization\":{},\"site_support\":{},\"nc_zero_score_numerator\":{},\"edge_mask_nt\":{},",
                "\"score_threshold_per_tr\":30,\"fold_threshold_per_tr\":3,\"replicate_rule\":\"sum, not individual pass\",",
                "\"count_meaning\":{},",
                "\"position_meaning\":\"candidate one nucleotide upstream of physical read end along transcript direction; mRNA is genome-mapped\",",
                "\"joint_evaluation\":{},\"considered_sites\":{},\"emitted_rows\":{},\"candidate_sites\":{},\"inputs\":[{}],\"samples\":[{}]}}\n"
            ),
            json(env!("CARGO_PKG_VERSION")),
            json(options.mode.name()),
            json(if options.format == Format::Unified {
                "unified"
            } else {
                "legacy"
            }),
            json(if options.rows == Rows::All {
                "all"
            } else {
                "candidates"
            }),
            if reference.is_empty() {
                "null".into()
            } else {
                json(reference)
            },
            json(if options.mode == Mode::Mrna {
                "nc_on_each_tr_observed_sites"
            } else {
                "whole_reference"
            }),
            json(if options.mode == Mode::Mrna {
                "all_tr_positive"
            } else {
                "full_reference_including_zero"
            }),
            if options.mode == Mode::Mrna { 1 } else { 0 },
            if options.mode == Mode::Rrna { 20 } else { 0 },
            json(
                if options.format == Format::Legacy && options.mode == Mode::Rrna {
                    "legacy count is the raw local-ratio denominator at candidate position; score signal is the next reference position"
                } else {
                    "signal end count; adjacent_count is the raw local-ratio denominator"
                }
            ),
            json(if valid {
                "available"
            } else {
                "normalization_unavailable"
            }),
            considered,
            rows.len(),
            rows.iter().filter(|r| r.candidate()).count(),
            inputs_json,
            normal
        );
        fs::write(path, data).map_err(|e| e.to_string())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_ratio_and_overlap_normalization() {
        let norm = Normalization { nc: 21, tr: 246 };
        let plus = pair_score(
            Local {
                signal: 80,
                adjacent: 2,
            },
            Local {
                signal: 10,
                adjacent: 5,
            },
            norm,
            Mode::Mrna,
        )
        .unwrap();
        let minus = pair_score(
            Local {
                signal: 160,
                adjacent: 4,
            },
            Local {
                signal: 4,
                adjacent: 2,
            },
            norm,
            Mode::Mrna,
        )
        .unwrap();
        assert_eq!((plus.score, minus.score), (38.0, 38.0));
        assert_eq!((plus.corrected, minus.corrected), (7.0, 14.0));
        let site = Site {
            chrom: "1".into(),
            position: 20,
            strand: '-',
        };
        assert_eq!(adjacent_site(&site).unwrap().position, 21);
        let site = Site {
            strand: '+',
            ..site
        };
        assert_eq!(adjacent_site(&site).unwrap().position, 19);
    }
    #[test]
    fn zero_totals_are_unevaluable_in_both_modes() {
        for mode in [Mode::Mrna, Mode::Rrna] {
            for norm in [
                Normalization { nc: 0, tr: 40 },
                Normalization { nc: 40, tr: 0 },
            ] {
                let row = evaluate(
                    vec![],
                    Local {
                        signal: 0,
                        adjacent: 0,
                    },
                    vec![Local {
                        signal: 100,
                        adjacent: 1,
                    }],
                    &[norm],
                    mode,
                    false,
                );
                assert_eq!(row.status, "normalization_unavailable");
                assert!(!row.candidate());
                assert!(row.fold.is_none() && row.total_score.is_none());
            }
        }
    }
    #[test]
    fn zero_current_nc_is_policy_specific_when_normalization_is_valid() {
        let t = Local {
            signal: 30,
            adjacent: 0,
        };
        let n = Local {
            signal: 0,
            adjacent: 0,
        };
        let norm = Normalization { nc: 100, tr: 100 };
        assert_eq!(pair_score(t, n, norm, Mode::Mrna).unwrap().score, 29.0);
        assert_eq!(pair_score(t, n, norm, Mode::Rrna).unwrap().score, 30.0);
    }
    #[test]
    fn one_unusable_replicate_invalidates_joint_call() {
        let row = evaluate(
            vec![],
            Local {
                signal: 1,
                adjacent: 1,
            },
            vec![
                Local {
                    signal: 100,
                    adjacent: 1
                };
                2
            ],
            &[
                Normalization { nc: 100, tr: 100 },
                Normalization { nc: 0, tr: 100 },
            ],
            Mode::Mrna,
            false,
        );
        assert!(!row.candidate());
        assert_eq!(row.status, "normalization_unavailable");
    }
    #[test]
    fn boundary_exclusion_has_no_fabricated_zero_score_in_unified_output() {
        let row = evaluate(
            vec!["ref".into(); 7],
            Local {
                signal: 1,
                adjacent: 1,
            },
            vec![Local {
                signal: 100,
                adjacent: 1,
            }],
            &[Normalization { nc: 100, tr: 100 }],
            Mode::Rrna,
            true,
        );
        assert_eq!(row.status, "boundary_excluded");
        assert!(row.total_score.is_none());
        assert!(!row.candidate());
        assert!(
            unified_values(&row, &[Normalization { nc: 100, tr: 100 }], Mode::Rrna)
                .contains(&"not_evaluated".into())
        );
    }
}
