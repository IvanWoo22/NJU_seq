# Scoring formats and rules

**English** | [简体中文](scoring.zh-CN.md)

`nju-seq score` jointly evaluates one NC against all TR samples. mRNA and rRNA share an output format but differ in input counts, normalization scope, and candidate coordinates.

## Options

Build from the repository root:

```bash
cargo build --release --offline
target/release/nju-seq --version
```

Command interface:

```text
nju-seq score --mode mrna|rrna --nc NC.tsv --tr TR1.tsv [--tr TR2.tsv ...]
  --output scores.tsv [--summary scores.json]
  [--format unified|legacy] [--rows candidates|all]
  [--reference-id 18s]
```

`--mode` defaults to `mrna`; `--format` defaults to `unified`. rRNA requires a `--reference-id` without whitespace, such as `18s`, `28s`, or `5-8s`. The mode controls both input parsing and scoring rules. Four-column rRNA counts cannot be passed directly to mRNA mode.

Sample names come from the parent directories of the input files. Put the NC and each TR in distinct sample directories, such as `NJU6096/rrna_18s.tsv`. The unified format rejects duplicate column names. The output and optional summary must not already exist; invalid inputs, options, or failed calculations produce no result files.

| Option | mRNA | rRNA |
|---|---|---|
| Input | Nine-column genomic read-end counts without a header; end counts must be positive integers | Four columns without a header: `Position/Base/StartCount/EndCount`; positions are consecutive from 1 to L, and counts may be zero |
| Default rows | `candidates`: joint candidates only | `all`: candidate coordinates 1…L−1 |
| `--rows all` | Intersection of sites with positive counts in every TR; not all genomic positions | Positions 1…L−1, including rows below thresholds or unavailable for evaluation |
| `--rows candidates` | Only rows with `Candidate=true` | Only rows with `Candidate=true` |

rRNA inputs must have the same length and base at every position. Empty references, missing or duplicate positions, and negative counts are errors. A reference of length 1 produces only a header; references of length 40 or less have no interior positions outside the boundary mask.

Provide `--summary` to distinguish unavailable normalization, an empty intersection, and successful evaluation with no candidates when the score table is empty. The summary also includes version, mode, row selection, rules, and each TR's NC/TR totals and normalization factor. `joint_evaluation=available` means normalization is available; it does not guarantee an interior position can be scored.

`nju-seq run` saves read and optional UMI score summaries automatically. See [rRNA counting and scoring](rrna.md) for the separate rRNA commands.

## Formulas

Let `signal` be the raw read-end count assigned to a candidate, `adjacent` the neighboring raw count used in the local ratio, and `n` the number of TR samples. Both modes use:

```text
score_j       = TR_signal_j / max(TR_adjacent_j, 1)
                − NC_score_numerator / max(NC_adjacent, 1)
corrected_j   = ceil(TR_signal_j × NC_TOTAL_j / TR_TOTAL_j)
TotalScore    = sum(score_j)
AveScore      = TotalScore / n
TotalCorrectedCount = sum(corrected_j)
TotalFoldChange     = TotalCorrectedCount / max(NC_signal, 1)
AveFC               = TotalFoldChange / n

Candidate = evaluable and outside the boundary mask
            and TotalScore >= 30 × n
            and TotalCorrectedCount >= 3 × n × NC_signal
            and TotalCorrectedCount > 0
```

`corrected_j` is calculated by multiplying first, dividing second, and rounding up last, using the unrounded normalization factor. The joint score threshold applies to the sum across all TR samples; an individual TR score may be below 30.

**A zero NC count at one site is different from a zero NC normalization total.** When the site count is zero but normalization is available, the fold-change criterion still uses the raw NC=0 and requires a positive corrected sum. In this case, `FC_criteria=meet` need not be equivalent to `AveFC>=3`. The displayed average fold change alone cannot reconstruct the candidate decision.

| Rule | `--mode mrna` | `--mode rrna` |
|---|---|---|
| NC_TOTAL | Sum NC counts at the sites detected in each TR, separately for that TR | Sum NC counts across the entire current rRNA reference |
| TR_TOTAL | All read-end counts in that TR input | Sum TR counts across the entire current rRNA reference |
| Sites eligible for scoring | Intersection of positive-count sites across all TR samples | Full reference; a TR may have zero counts at an individual site |
| NC score numerator when NC signal is zero | 1 | 0 |
| Adjacent count | Candidate genomic coordinate −1 on the positive strand, +1 on the negative strand | Reference position k for a candidate at k; the numerator uses position k+1 |
| Boundary mask | None added | Exclude candidate k≤20 or k>L−20; the final position has no k+1 and is omitted |

At exon junctions, mRNA scoring still uses neighboring genomic coordinates. These counts need not correspond to adjacent positions in the transcript.

## Unavailable scores

Every TR must have `NC_TOTAL>0` and `TR_TOTAL>0`. If either total is unavailable for any TR, the entire joint evaluation is marked `normalization_unavailable`. The program does not drop that TR, reduce n, or allow strong signals in other TR samples to compensate for a library with no counts.

| `EvaluationStatus` | Meaning | Scores, corrected counts, total/average FC | `Candidate` / threshold fields |
|---|---|---|---|
| `evaluated` | Normalization is available and the position passes boundary checks | Numeric values, which may be zero | true/false; meet/unmeet |
| `normalization_unavailable` | At least one required TR-specific total is zero | `NA` | false; not_evaluated |
| `boundary_excluded` | Position lies within the rRNA boundary mask | `NA` | false; not_evaluated |

Unavailable normalization takes precedence over boundary exclusion. Raw signal, adjacent counts, and totals are retained to explain the status. Normalization factors are reported individually and are `NA` only when their own totals are invalid. In plots, distinguish unavailable evaluations from valid zero scores.

## Unified output fields

Both modes use the same header. With one NC and three TR samples, the table has 46 columns covering coordinates, annotations, counts, scores, normalization quantities, and evaluation status.

| Field | Definition |
|---|---|
| `Chromosome / Position / Direction` | Candidate position, 1-based; genomic coordinates for mRNA, coordinates within the explicit reference ID and strand `+` for rRNA |
| `CoordinateSystem` | `genome_1based` or `rrna_reference_1based`; rRNA positions are not hg38 coordinates |
| `LeftBase / Base / RightBase / GeneID` | Input annotations for mRNA; neighboring reference bases and the reference ID as GeneID for rRNA |
| `<NC>_count(NC)`, `<TR>_count` | Raw **signal** counts assigned to the displayed candidate position |
| `<sample>_adjacent_count` | Raw count at the denominator position of the local ratio |
| `<NC>_score_numerator` | Numerator actually used in the NC local ratio, including the mode-specific substitution of 1 for zero |
| `<sample>_score_denominator` | Denominator actually used: `max(adjacent_count,1)` |
| `<TR>_count_corr / <TR>_score` | Normalized count for a TR and the difference between its local ratio and the NC ratio |
| `<TR>_NC_TOTAL / TR_TOTAL / normalization_factor` | Totals used for that TR and their ratio `NC_TOTAL/TR_TOTAL` |
| `TotalFoldChange / AveFC` | Sum of corrected TR fold changes / their arithmetic mean |
| `TotalScore / AveScore` | Sum of TR scores / their arithmetic mean |
| `TotalCorrectedCount` | Sum of normalized counts across all TR samples |
| `ScoreMode / EvaluationStatus / Candidate` | Scoring mode, evaluation status, and final joint candidate decision |
| `FC_criteria / Score_criteria` | Outcome of each joint threshold; not_evaluated when evaluation is unavailable |

The candidate is the base immediately preceding the physical read end in the direction of transcription. In the `unified` format, an rRNA row at k reports **read-end counts from position k+1** in `*_count` and counts from position k in `*_adjacent_count`. These are the numerator and denominator counts of the local ratio, respectively.

### Legacy format

Use `--format legacy` when downstream scripts require the header and column order of Perl `score.pl`. mRNA outputs candidates only and does not support `--rows all` in this format. rRNA uses a `#Position` column; count values describe physical read ends at the displayed position k, with numeric values formatted to two decimal places. `FoldChange` corresponds to `TotalFoldChange` in the unified format.

Both formats use the same normalization requirements. Legacy rRNA tables represent unavailable evaluations as 0/unmeet. Use `unified` or a `--summary` file to distinguish unavailable normalization from boundary exclusion. Check the required column names and coordinate definitions before selecting plotting or annotation scripts.
