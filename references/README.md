# Building reference files

**English** | [简体中文](README.zh-CN.md)

Run the commands below from the repository root after [installing the environment](../README.md#1-installation). They download public source files, select the required sequences, convert annotations, check transcript sequences against genomic exons, and build Bowtie2 indices. Complete references are generated locally; Git contains the scripts, source checksums, and small sequence/interval definitions.

## Files created by the commands

| Directory | Generated contents | Command |
|---|---|---|
| `downloads/` | Compressed public sources, checked against pinned SHA-256 values | All three commands below |
| `bacterial/` | `bac_rna.fa`, six `bac_rna.*.bt2` files, and `audit.json` | `build_references.py bacterial` |
| `human_source/` | `hsa_protein_coding.fa`, `hsa_exon.info`, GRCh38 genome, and `sources.json` | `build_references.py human` |
| `human_rebuilt/` | `mrna`, `rrna`, and `bac` index prefixes; `exon.info`; FASTA files and source information | `build_references.py human` |
| `arabidopsis_source/` | Coding transcript FASTA, exon table, TAIR10 genome, and expanded rRNA FASTA | `build_references.py arabidopsis` |
| `arabidopsis/` | `ath_protein_coding`, `ath_rrna_total`, and `bac` index prefixes; `exon.info`; FASTA files and source information | `build_references.py arabidopsis` |

Paths in this table are relative to `references/`. The output and source directories must be new. If they already exist, select new paths with `--out` and `--source-out` and use those paths in subsequent analysis commands. The `bac` prefix is omitted when building [without bacterial rRNA](#building-without-bacterial-rrna).

## 1. Bacterial rRNA

```bash
python3 scripts/build_references.py bacterial \
  --download-dir references/downloads \
  --out references/bacterial --threads 4
```

This downloads the **SILVA 138.2 SSURef NR99 and LSURef NR99 truncated FASTA** files and selects entries whose taxonomy contains one of three complete family names: `Enterobacteriaceae`, `Lactobacillaceae`, or `Streptococcaceae`. It converts RNA `U` to DNA `T`, retains source order (SSU followed by LSU), and rejects duplicate sequence IDs. The result contains **34,982 SSU + 21,135 LSU = 56,117 sequences**. This is a defined bacterial rRNA subset, not the whole SILVA database.

The script verifies the selected FASTA checksum before running `bowtie2-build`. The index prefix used by the analysis is `references/bacterial/bac_rna`. Human and Arabidopsis workflows require this index only when bacterial rRNA removal is enabled (`--rrna-filter both` or `bac`).

## 2. Human: GENCODE v34 / GRCh38

```bash
python3 scripts/build_references.py human \
  --download-dir references/downloads \
  --source-out references/human_source \
  --bac-index references/bacterial/bac_rna \
  --out references/human_rebuilt --threads 4 --copy-inputs
```

The script downloads the v34 transcript FASTA, comprehensive GFF3 annotation, and GRCh38 primary-assembly genome. It selects FASTA records whose eighth `|`-separated field is `protein_coding`, then selects their GFF3 `exon` features. GFF3 columns 1, 4, 5, 7, and transcript/gene attributes become the five-column exon table. Noncoding transcripts and representative-transcript filtering are not part of this selection.

The resulting mRNA reference contains **84,068 transcripts and 762,910 exon rows**. Every transcript must match its genomic exons, including negative-strand reverse complements and `_PAR_Y` identifiers. The three host rRNA inputs are [18s.fa](../data/hsa_rrna/18s.fa), [28s.fa](../data/hsa_rrna/28s.fa), and [5-8s.fa](../data/hsa_rrna/5-8s.fa).

Use `--index references/human_rebuilt --exons references/human_rebuilt/exon.info` with `nju-seq`; the [main README](../README.md#6-human-example-hek293t) supplies the complete human example.

## 3. Arabidopsis: Ensembl Plants 46 / TAIR10

```bash
python3 scripts/build_references.py arabidopsis \
  --download-dir references/downloads \
  --source-out references/arabidopsis_source \
  --bac-index references/bacterial/bac_rna \
  --out references/arabidopsis --threads 4 --copy-inputs
```

The mRNA reference uses the release-46 TAIR10 genome, cDNA FASTA, and GFF3. The script selects `transcript_biotype:protein_coding` records, retains matching exon features, and checks all **48,321 transcripts** against the genome. It preserves the `Parent=transcript:AT…` annotation syntax used by the Rust workflow.

The **4,131-sequence rRNA reference** is built separately from:

- Four rRNA sequences in [data/ath_rrna](../data/ath_rrna), with the explicit substitutions and insertions in [arabidopsis_rrna.json](../config/arabidopsis_rrna.json). Input and output sequence checksums are checked.
- **4,127 genomic intervals** from the public **Col-PEK assembly GWHBFWQ00000000.1**, listed in [arabidopsis_rrna.intervals.bed](../config/arabidopsis_rrna.intervals.bed). These intervals include rRNA copies and fragments and are used for rRNA removal.

The BED file uses **0-based, half-open** coordinates and genomic forward-strand sequences. Its fourth column supplies stable FASTA IDs. Chr4 IDs retain the 100-base leading-gap offset of the Col-PEK naming scheme; use BED columns 2–3 for extraction coordinates instead of parsing the IDs. These are rRNA reference labels, not TAIR10 mRNA coordinates. The complete reconstructed rRNA FASTA must match its pinned SHA-256 checksum before index construction.

The output prefixes are `ath_protein_coding` and `ath_rrna_total`. Continue with the [Arabidopsis UMI workflow](../docs/arabidopsis_umi.md#run-a-small-subset). The counting command uses the generated `exon.info`; the genome in the package is named `genome.fa`.

## Building without bacterial rRNA

To remove host rRNA only, skip the bacterial build and pass `--without-bacterial` instead of `--bac-index` to either host builder. For example:

```bash
python3 scripts/build_references.py human \
  --download-dir references/downloads \
  --source-out references/human_host_source \
  --without-bacterial \
  --out references/human_host --threads 4 --copy-inputs

target/release/nju-seq run \
  --samples tests/data/human/samples.tsv --from raw --rrna-filter host \
  --index references/human_host --exons references/human_host/exon.info \
  --out output/human_host_example --threads 2 --jobs 1
```

This downloads only the human sources and builds `mrna` and `rrna`, plus `exon.info`. It neither downloads SILVA nor requires an existing bacterial index. The same build option applies to `arabidopsis`; its host prefixes remain `ath_protein_coding` and `ath_rrna_total`. In either workflow, select `--rrna-filter host` during analysis to skip bacterial alignment and removal. See [all filtering choices](../README.md#rrna-filtering), including `none` for direct mRNA alignment.

## Downloads, storage, and checks

Exact URLs and SHA-256 values are in [reference_sources.json](../config/reference_sources.json); the [reference source list](../docs/reference_sources.md) also includes mouse GENCODE M25 and optional references. Downloaded data retain their original filenames. To download without building, for example:

```bash
python3 scripts/build_references.py human \
  --download-dir references/downloads --download-only
```

Run this separately for `bacterial` and `arabidopsis` if needed. Completed downloads are reused only after checksum verification. Add `--offline` to a build command to require a complete local cache. To use files downloaded elsewhere, put them under the exact filenames listed in the JSON manifest; the same checksum checks apply.

`--copy-inputs` copies the genome and any included bacterial index into each host package, allowing that package to move independently. Without it, relative symlinks are created, so preserve the relative locations of the source and reference directories. The human build needs several GB of disk space for the decompressed genome and indices; copying inputs requires additional space. Complete references and generated indices are excluded by `.gitignore`.

Each package includes source/sequence checksums, reference counts, and Bowtie2 build logs. A failed build leaves a `.partial` directory for inspection. Correct the problem and choose new output/source directories; successful packages are never overwritten. Download failures leave a `.partial` file and a subsequent attempt downloads that file again.

## Using separately prepared human inputs

The lower-level script accepts matching source files directly:

```bash
python3 scripts/prepare_human_reference.py \
  --mrna references/human_source/hsa_protein_coding.fa \
  --exons references/human_source/hsa_exon.info \
  --genome references/human_source/GRCh38.primary_assembly.genome.fa \
  --bac-index references/bacterial/bac_rna \
  --out references/human_custom --threads 4 \
  --label 'GENCODE v34 / GRCh38' --copy-inputs
```

Transcript FASTA headers use the GENCODE `|`-separated format. The first two fields are versioned ENST and ENSG identifiers; the eighth is `protein_coding`. The exon table has five tab-separated columns without a header:

```text
chromosome    start    end    strand    attributes
chr1          101      200    +         Parent=ENST000001.1;gene_id=ENSG000001.1
```

Spaces above are for display only. Exon coordinates are **1-based and inclusive**. Nine-column GFF3 cannot be passed directly to `--exons`; `build_references.py` performs this conversion. The genome must be uncompressed; transcript and rRNA FASTA inputs may be gzipped. If supplied, the bacterial index must contain all six `.bt2` or `.bt2l` files in one format. Omit `--bac-index` from this lower-level command when running with `--rrna-filter host` or `none`. Repeat `--rrna FILE` to use custom human rRNA sequences.
