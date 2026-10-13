# Reference files

**English** | [简体中文](README.zh-CN.md)

`data/` contains small references bundled with the repository, `references/` holds complete references and indices, and `inputs/` holds sequencing FASTQ files or symlinks. Prepare complete references and original FASTQ files separately.

Download links for human GENCODE v34, mouse M25, and Arabidopsis Ensembl Plants 46 are in [Reference sources](../docs/reference_sources.md). For Arabidopsis indices and libraries with UMIs, see the [Arabidopsis workflow](../docs/arabidopsis_umi.md).

## Directory layout

| Directory | Contents |
|---|---|
| `human/` | Prepared human mRNA, rRNA, and bacterial rRNA indices, plus `exon.info` |
| `human_source/` | Human transcript FASTA, five-column exon table, GRCh38 genome, and other annotations |
| `bacterial/` | Six files comprising the bacterial rRNA Bowtie2 index |
| `arabidopsis/` | Arabidopsis references and indices |

The mRNA reference for the human example contains 84,068 protein-coding transcripts and 762,910 exon rows. Its rRNA references are the three sequences in `data/hsa_rrna/18s.fa`, `28s.fa`, and `5-8s.fa`.

The transcript FASTA, exon table, and genome must match. Basic annotations, representative transcripts, region tables, and codon tables in `human_source/` support downstream annotation. The main workflow uses the full transcript set and its matching exon table.

`data/` also contains rRNA references for other species, plus `MCMV.fa` and `MHV.fa`. For another species, prepare matching transcripts and annotations and use an annotation parser that supports their format.

## Building and copying references

Run `scripts/prepare_human_reference.py` as shown in the [main README](../README.md#2-human-references). It creates a new directory containing the `mrna`, `rrna`, and `bac` index prefixes, `exon.info`, FASTA files, and `audit.json`.

- With `--copy-inputs`, the genome and bacterial index are copied, so the resulting directory can be moved independently.
- Without this option, relative symlinks are used; preserve the relative locations of the reference and source directories.
- To use other host rRNA sequences, repeat `--rrna FILE` for each FASTA input.

The script checks each transcript against its sequence reconstructed from genomic exons. On failure, it retains a `.partial` directory for inspection; correct the inputs and use a new output directory. The prepared `human/` directory also uses relative symlinks: copy it together with `human_source/` and `bacterial/`, or dereference symlinks when copying it alone, for example with `cp -RL references/human /destination/human`.

## Input format for human reference preparation

Transcript FASTA headers follow the GENCODE `|`-separated format. The first two fields are versioned ENST and ENSG identifiers; the eighth is `protein_coding`. `_PAR_Y` identifiers are supported.

The exon table has five columns without a header, separated by actual tabs:

```text
chromosome    start    end    strand    attributes
chr1          101      200    +         Parent=ENST000001.1;gene_id=ENSG000001.1
```

Spaces above are for display only. Coordinates are 1-based and inclusive. The script selects transcripts present in the input FASTA and preserves gene associations, negative strands, and PAR-Y information. Convert nine-column GFF3 to this format before passing it to `--exons`.

The genome must be an uncompressed FASTA; transcript and rRNA FASTA files may be gzipped. The bacterial index must contain all six `.bt2` or `.bt2l` files, using one format consistently.
