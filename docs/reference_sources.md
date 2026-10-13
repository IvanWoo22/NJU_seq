# Reference sources: human, mouse, and Arabidopsis

**English** | [简体中文](reference_sources.zh-CN.md)

This page lists reference releases, download links, and filenames for human, mouse, and Arabidopsis. The abbreviations `hsa`, `mmu`, and `ath` identify the three species. Store complete downloaded references under `references/`; small references bundled with the repository are in `data/`. [Building reference files](../references/README.md) provides executable commands from these downloads to complete indices. The fixed URLs and SHA-256 values used by the builder are in [reference_sources.json](../config/reference_sources.json).

| Species | Release | Genome assembly | Official source |
|---|---|---|---|
| Human (*Homo sapiens*) | GENCODE v34 | GRCh38.p13 | [GENCODE release 34](https://www.gencodegenes.org/human/release_34.html) |
| Mouse (*Mus musculus*) | GENCODE M25 | GRCm38.p6 | [GENCODE release M25](https://www.gencodegenes.org/mouse/release_M25.html) |
| Arabidopsis (*Arabidopsis thaliana*) | Ensembl Plants release 46 | TAIR10 | [Release 46 directory](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/) |

## Human: GENCODE v34

| Contents | File and download link | Local shorthand |
|---|---|---|
| GFF3 annotation | [gencode.v34.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gff3.gz) | `hsa.gff3.gz` |
| All transcript sequences | [gencode.v34.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz) | `hsa_transcript.fa.gz` |

`build_references.py human` selects `protein_coding` transcripts, converts their GFF3 exons, and validates their sequences before building indices. It also downloads the matching genome FASTA: [GRCh38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz) is available in the same release directory.

## Mouse: GENCODE M25

| Contents | File and download link | Local shorthand |
|---|---|---|
| GFF3 annotation | [gencode.vM25.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz) | `mmu.gff3.gz` |
| All transcript sequences | [gencode.vM25.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz) | `mmu_transcript.fa.gz` |

To check transcript sequences against the genome, use [GRCm38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz) from the same directory. M25 uses GRCm38; its annotations and coordinates cannot be mixed with GRCm39.

The Rust workflow supports human ENST/ENSG and Arabidopsis AT annotation formats. Mouse ENSMUST/ENSMUSG identifiers require an extension to the annotation parser before use.

## Arabidopsis: Ensembl Plants release 46 / TAIR10

| Contents | File and download link | Local shorthand |
|---|---|---|
| Genome FASTA | [Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz) | `ath_dna.fa.gz` |
| GFF3 annotation | [Arabidopsis_thaliana.TAIR10.46.gff3.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz) | `ath.gff3.gz` |
| cDNA FASTA | [Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz) | `ath_transcript.fa.gz` |
| ncRNA FASTA | [Arabidopsis_thaliana.TAIR10.ncrna.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz) | `ath_ncrna.fa.gz` |

Use the genome, transcript, and annotation files from the same release-46 directory so that sequences and coordinates agree. `build_references.py arabidopsis` prepares these mRNA inputs automatically. The expanded rRNA removal index also uses [Col-PEK GWHBFWQ00000000.1](https://ngdc.cncb.ac.cn/gwh/Assembly/83527/show): [genome FASTA](https://download.cncb.ac.cn/gwh/Plants/Arabidopsis_thaliana_Col-PEK_GWHBFWQ00000000.1/GWHBFWQ00000000.1.genome.fasta.gz). Its extraction intervals and the four canonical rRNA sequence edits are supplied in [Building reference files](../references/README.md#3-arabidopsis-ensembl-plants-46--tair10).

## rRNA and other references

| Reference | Location or source | Usage |
|---|---|---|
| Human rRNA | [data/hsa_rrna](../data/hsa_rrna): 18s, 28s, 5-8s | The three sequences used by the human example |
| Mouse rRNA | [data/mmu_rrna](../data/mmu_rrna): 18s, 28s, 5-8s | Mouse rRNA reference sequences |
| Arabidopsis rRNA | [data/ath_rrna](../data/ath_rrna): 18s, 25s, 5-8s, 5s | Arabidopsis rRNA sequences; export the FASTA for counting from the actual alignment index, as described in the [UMI workflow](arabidopsis_umi.md#optional-rrna-analysis) |
| Bacterial rRNA | SILVA 138.2 [SSURef NR99 truncated FASTA](https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz) and [LSURef NR99 truncated FASTA](https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_LSURef_NR99_tax_silva_trunc.fasta.gz) | `build_references.py bacterial` selects Enterobacteriaceae, Lactobacillaceae, and Streptococcaceae and builds `references/bacterial/bac_rna` |
| miRNA (optional) | [miRBase downloads](https://www.mirbase.org/download/): [mature.fa](https://www.mirbase.org/download/mature.fa) | Save the release and file checksum; compress the FASTA if saving it as `mirna.fa.gz` |

Use matching transcript, GFF3, and genome releases and retain SHA-256 checksums for downloaded files. See [Reference files](../references/README.md) for directory layout and exon-table format.
