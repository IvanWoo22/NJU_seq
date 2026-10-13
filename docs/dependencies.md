# Dependencies

**English** | [简体中文](dependencies.zh-CN.md)

Create a Conda environment with [environment.yml](../environment.yml), using packages from [conda-forge](https://conda-forge.org/) and [Bioconda](https://bioconda.github.io/recipes/bowtie2/README.html). The commands below target Linux; Windows users can run them in a Linux environment.

On Linux x86_64, you can install pinned package builds, including indirect dependencies:

```bash
conda create -n nju-seq --file config/conda-linux-64.lock
conda activate nju-seq
```

Use either this method or the environment creation command in the [main README](../README.md#1-installation). The [lock file](../config/conda-linux-64.lock) fixes package builds and checksums. Use `environment.yml` on other platforms.

## Core analysis

Building `nju-seq` requires Rust / Cargo 1.88 or later. No third-party Rust crates need to be downloaded.

Starting from raw FASTQ requires Cutadapt, Bowtie2, and gzip. Bowtie2's launch scripts also need Python and Perl. Reference preparation additionally uses `bowtie2-build` and `samtools`. The Python helper scripts use the standard library and require Python ≥3.11.

The Conda environment includes these tools and a linker. For a manual installation, provide a C linker and check that the commands below use the intended environment:

```bash
cargo --version
rustc --version
python3 --version
cutadapt --version
bowtie2 --version
samtools --version
perl -e 'print "$^V\n"'
gzip --version
```

## UMI libraries

For `umi=r1_3p8` libraries, install UMI-tools in the active `nju-seq` environment:

```bash
conda install --freeze-installed --override-channels -c conda-forge -c bioconda umi_tools=1.1.6
umi_tools --version
```

UMI-tools extracts tags; Rust groups molecules by exact matching. See the [Arabidopsis UMI workflow](arabidopsis_umi.md) for commands.

## Perl analysis scripts

Install the modules required by the transcript-region annotation, Perl counting, and related scripts:

```bash
conda install --freeze-installed --override-channels -c conda-forge -c bioconda \
  perl-perlio-gzip=0.20 perl-scalar-util-numeric=0.40 \
  perl-app-cpanminus make
cpanm --mirror https://cpan.metacpan.org --mirror-only AlignDB::IntSpan@1.1.1

perl -MPerlIO::gzip -MAlignDB::IntSpan -e 'print "Perl modules OK\n"'
```

`AlignDB::IntSpan` handles intervals and annotations; `PerlIO::gzip` reads compressed files. The standalone rRNA counter, `readend_count.pl`, uses only Perl core modules.

Install [Scalar::Util::Numeric](https://bioconda.github.io/recipes/perl-scalar-util-numeric/README.html), a dependency of `AlignDB::IntSpan`, as a precompiled Conda package to avoid local Perl extension compiler configuration issues.

Tools are resolved through `PATH` by default. For `scripts/count_rrna.py`, use `--perl /path/to/perl` or the `NJU_PERL` environment variable to select an interpreter. Supply an actual path and install any modules in that Perl environment.

## Optional plotting

Open `qc.html` directly in a browser. For R plotting scripts, install the packages used by the script you need:

| Script | R packages |
|---|---|
| `presentation/seq_logo.R` | ggplot2, ggseqlogo |
| `presentation/point_venn.R` | VennDiagram |
| `presentation/point_distribution.R` | ggplot2, gridExtra, RColorBrewer, optparse |
| `presentation/signature_count.R` | ggplot2, reshape2, dplyr |
| `presentation/vertical_barplot.R` | ggplot2, dplyr, forcats, gridExtra |
| `presentation/codon_distance*.R` | ggplot2, gridExtra; splines is included with R |
| `presentation/splice_pie.R`, `tool/draw_table.R` | ggplot2, ggpubr |
| `quality_control/picture_draw.R` | gridExtra, ggpubr |
| `presentation/exon_distance.R`, `gene_ontology.R` | ggplot2 |

`grid` is included with R. Input tables, coordinate definitions, and annotation versions used for plotting and region analysis must match the analysis results.
