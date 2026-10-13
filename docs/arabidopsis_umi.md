# Arabidopsis Col_NC31 / Col_TR31: analysis with 8 nt UMIs

**English** | [简体中文](arabidopsis_umi.zh-CN.md)

Col_NC31 is the NC control and Col_TR31 is the TR treatment sample. The libraries contain 8 nt UMIs added **after MgR treatment and before PCR**, read from the 3′ end of R1. The 5′ end of R2 contains the same tag in the reverse-complement orientation. This guide describes the Rust `nju-seq` workflow, UMI counting rules, and analysis results for these samples.

[Setup](#software-and-inputs) · [References](#arabidopsis-references) · [Small subset](#run-a-small-subset) · [UMI processing](#umi-processing) · [Results](#sample-results) · [rRNA](#optional-rrna-analysis)

## Software and inputs

Work from the repository root. Create the `nju-seq` environment as described in the [main README](../README.md#1-installation), then install UMI-tools:

```bash
conda activate nju-seq
conda install --freeze-installed --override-channels -c conda-forge -c bioconda umi_tools=1.1.6
cargo build --release --offline
umi_tools --version
```

UMI-tools extracts the tags; Rust performs molecule grouping and counting. See [Dependencies](dependencies.md) for the required software.

Place the original paired FASTQ files or symlinks at the following paths, or edit the [sample manifest](../config/samples.arabidopsis_umi.tsv):

```text
inputs/arabidopsis/Col_NC31/R1.fq.gz
inputs/arabidopsis/Col_NC31/R2.fq.gz
inputs/arabidopsis/Col_TR31/R1.fq.gz
inputs/arabidopsis/Col_TR31/R2.fq.gz
```

Both manifest rows use `umi=r1_3p8` and `alignment=-`. Relative paths are resolved from the manifest's directory. For FASTQ files that have already been trimmed and had UMIs extracted, use `--from prepared` and retain the extracted `_UMI` suffix in both read IDs. Passing these files to `--from raw` would trim and extract them again.

## Arabidopsis references

Prepare the following matching references. Complete references and indices are not included in the Git repository.

| Path or index prefix | Purpose |
|---|---|
| `references/bacterial/bac_rna` | Bacterial rRNA index |
| `references/arabidopsis/ath_rrna_total` | Arabidopsis rRNA index containing 4,131 sequences |
| `references/arabidopsis/ath_protein_coding` | mRNA index containing 48,321 transcripts |
| `references/arabidopsis/exon.info` | Five-column exon table matching the mRNA index |

Each index prefix has six `.bt2` files; copy all six. `ath.fa`, `ath_main_transcript.txt`, `transcript_region.info`, and `exon_distance.info` support downstream sequence and region annotation and are not read by the main counting command.

Download links for Ensembl Plants 46 / TAIR10 are in [Reference sources](reference_sources.md). The commands below use the indices and matching exon table above. `prepare_human_reference.py` accepts human GENCODE data; use prepared indices for Arabidopsis.

## Run a small subset

Start with the first 1,000 read pairs from each sample. The subsetting script preserves the UMI setting and checks paired read IDs within the selected range.

```bash
python3 scripts/subset_fastq.py \
  --samples config/samples.arabidopsis_umi.tsv \
  --pairs 1000 --skip 0 --out output/arabidopsis_umi_subset

target/release/nju-seq plan \
  --samples output/arabidopsis_umi_subset/samples.tsv --from raw \
  --bac-index references/bacterial/bac_rna \
  --rrna-index references/arabidopsis/ath_rrna_total \
  --mrna-index references/arabidopsis/ath_protein_coding \
  --exons references/arabidopsis/exon.info \
  --out output/arabidopsis_umi_example --threads 2 --jobs 1

target/release/nju-seq run \
  --samples output/arabidopsis_umi_subset/samples.tsv --from raw \
  --bac-index references/bacterial/bac_rna \
  --rrna-index references/arabidopsis/ath_rrna_total \
  --mrna-index references/arabidopsis/ath_protein_coding \
  --exons references/arabidopsis/exon.info \
  --out output/arabidopsis_umi_example --threads 2 --jobs 1
```

`plan` checks paths and lists commands without running the analysis. The subsetting and analysis output directories must not already exist. A small subset lets you check the installation, paired-read handling, and output format; it does not support estimates of full-library candidate counts. Insufficient coverage can produce `normalization_unavailable`.

For the full dataset, keep the index arguments, change `--samples` to `config/samples.arabidopsis_umi.tsv`, and set `--out` to a new directory such as `output/arabidopsis_umi_full`. `--threads` sets threads per alignment task; `--jobs` sets concurrent samples. Intermediate FASTQ files are removed after success by default. Add `--storage all` to retain them and `--qc yes` if you also want a QC report in that mode. See [Helper scripts](scripts.md) for chunked alignment and resource limits.

## UMI processing

`umi=r1_3p8` enables the following steps, performed separately for each library:

1. Cutadapt uses R1 adapter `AGATCGGAAGAGCACA` and R2 adapter `GATCGTCGGACTGTAG`, with `-U 8` to remove the tag at the start of R2. The remaining options are `-O 6 -m 18:10 -e 0.1 --discard-untrimmed`. The minimum R1 length includes its 8 nt UMI, which has not yet been extracted.
2. UMI-tools extracts the last 8 nt of R1 with `--extract-method=string --bc-pattern=NNNNNNNN --3prime` and appends the tag to both read IDs. Extraction alone does not group molecules. The R2 tag is not used for joint error correction.
3. Reads mapping to bacterial rRNA and then Arabidopsis rRNA are removed before mRNA alignment. Both rRNA steps use `-k 5`; mRNA uses `-a`. All retain `--reorder --seed 0`.
4. After mRNA alignment selection, the workflow keeps a read-count branch and groups molecules by an exact match on **chromosome, strand, the complete genomic exon path of R1, and UMI**. Different endpoints or splice paths remain separate; the same fragment mapped to equivalent transcripts can be grouped.
5. Read and UMI end counts are generated separately, merged at genomic coordinates, and jointly scored between NC and TR. Joint UMI scoring is produced only when all samples in the manifest have UMIs.

Each UMI family is represented by the read with the fewest CIGAR errors, with ties resolved by read ID. UMIs containing `N` are retained individually to avoid grouping unknown tags. Rust uses exact matching without directional error correction. Grouping uses the full key above, rather than the UMI string or end coordinate alone.

| Output | Meaning |
|---|---|
| `<sample>/mrna.almostuniquematch.filter.tmp` | Read alignments passing mRNA selection |
| `<sample>/mrna.umi.filter.tmp` | One representative alignment per UMI family |
| `<sample>/mrna.umi.families.tsv.gz`, `mrna.umi.summary.json` | Family membership and grouping statistics, both in the sample directory |
| `<sample>/reads/mrna.filter.tsv` | Read counts |
| `<sample>/umi/mrna.filter.tsv` | UMI counts |
| `mrna.reads.scored.tsv`, `mrna.umi.scored.tsv` | Scores for the two counting methods; candidates are identified by the `Candidate` column |
| `mrna.reads.score_summary.json`, `mrna.umi.score_summary.json` | Normalization status and candidate counts |
| `qc.html`, `SUCCESS` | Quality report and successful completion marker |

## Sample results

The following full-dataset results for Col_NC31 / Col_TR31 were calculated with the Perl/Python workflow using exact UMI grouping.

| Metric | Col_NC31 | Col_TR31 |
|---|---:|---:|
| Raw read pairs | 89,597,181 | 38,795,779 |
| Read pairs after adapter trimming | 31,663,953 | 12,546,074 |
| Read pairs after UMI extraction | 31,663,953 | 12,546,074 |
| Reads passing mRNA selection | 3,239,437 | 575,214 |
| Count records after exact UMI grouping | 874,050 | 395,251 |
| Of these, records with an `N` in the UMI, retained individually | 545 | 9,284 |

The read branch yielded 134 candidates and the exact UMI branch yielded 110: 107 were shared, 27 occurred only in the read branch, and 3 only in the UMI branch. These differences reflect the counting method and do not by themselves demonstrate improved accuracy. There is one NC and one TR sample; the scoring thresholds do not constitute a statistical test across replicates or provide FDR control. See [Scoring](scoring.md).

## Optional rRNA analysis

The main workflow removes rRNA for mRNA analysis and does not produce rRNA UMI molecule counts. For rRNA sites, follow [rRNA counting and scoring](rrna.md), replacing the index and FASTQ paths with the Arabidopsis inputs.

Count against the same reference sequences used for alignment. Export them with `bowtie2-inspect references/arabidopsis/ath_rrna_total`, save the `18s`, `25s`, `5-8s`, and `5s` reference IDs as separate single-sequence FASTA files, and pass each to `scripts/count_rrna.py` with a repeated `--reference FILE` option. Preserve the original FASTA IDs so that alignments and counting coordinates agree.
