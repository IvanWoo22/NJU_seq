# 软件依赖

按 [environment.yml](../environment.yml) 建立 Conda 环境，软件包来自 [conda-forge](https://conda-forge.org/) 和 [Bioconda](https://bioconda.github.io/recipes/bowtie2/README.html)。以下命令适用于 Linux；Windows 用户可在 Linux 环境中运行。

Linux x86_64 可使用固定包版本安装，包括间接依赖：

```bash
conda create -n nju-seq --file config/conda-linux-64.lock
conda activate nju-seq
```

这与 README 的环境创建方式二选一。[版本清单](../config/conda-linux-64.lock) 固定了包构建和校验值；其他平台应使用 `environment.yml`。

## 核心分析

编译 `nju-seq` 需要 Rust / Cargo 1.88 或更高版本，无须下载第三方 Rust crate。

运行原始 FASTQ 需要 Cutadapt、Bowtie2 和 gzip；Bowtie2 的启动脚本还需要 Python、Perl。准备参考另需 `bowtie2-build` 和 `samtools`。Python 辅助脚本使用标准库，要求 Python ≥3.11。

Conda 环境已包含上述工具和链接器。自行安装时还需提供 C 链接器，并确保以下命令来自预期环境：

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

## UMI 文库

分析 `umi=r1_3p8` 文库时，在已激活的 `nju-seq` 环境中安装 UMI-tools：

```bash
conda install --freeze-installed --override-channels -c conda-forge -c bioconda umi_tools=1.1.6
umi_tools --version
```

UMI-tools 负责提取标签，Rust 程序按 exact 规则合并分子。操作步骤见 [拟南芥 UMI 分析](arabidopsis_umi.md)。

## Perl 分析脚本

使用转录本区域注释、Perl 计数等脚本时，安装对应模块：

```bash
conda install --freeze-installed --override-channels -c conda-forge -c bioconda \
  perl-perlio-gzip=0.20 perl-scalar-util-numeric=0.40 \
  perl-app-cpanminus make
cpanm --mirror https://cpan.metacpan.org --mirror-only AlignDB::IntSpan@1.1.1

perl -MPerlIO::gzip -MAlignDB::IntSpan -e 'print "Perl modules OK\n"'
```

`AlignDB::IntSpan` 用于区间和注释处理，`PerlIO::gzip` 用于读取压缩文件。独立 rRNA 计数器 `readend_count.pl` 只需要 Perl 核心模块。

`AlignDB::IntSpan` 的依赖 [Scalar::Util::Numeric](https://bioconda.github.io/recipes/perl-scalar-util-numeric/README.html) 使用 Conda 提供的预编译包，避免本地 Perl 扩展编译器配置影响安装。

默认使用 `PATH` 中的工具。`scripts/count_rrna.py` 可用 `--perl /path/to/perl` 或环境变量 `NJU_PERL` 指定解释器；请使用实际路径，并确保模块安装在对应的 Perl 环境中。

## 可选绘图

`qc.html` 可直接用浏览器打开。运行 R 绘图脚本时，按所用脚本安装包：

| 脚本 | R 包 |
|---|---|
| `presentation/seq_logo.R` | ggplot2、ggseqlogo |
| `presentation/point_venn.R` | VennDiagram |
| `presentation/point_distribution.R` | ggplot2、gridExtra、RColorBrewer、optparse |
| `presentation/signature_count.R` | ggplot2、reshape2、dplyr |
| `presentation/vertical_barplot.R` | ggplot2、dplyr、forcats、gridExtra |
| `presentation/codon_distance*.R` | ggplot2、gridExtra；splines 随 R 提供 |
| `presentation/splice_pie.R`、`tool/draw_table.R` | ggplot2、ggpubr |
| `quality_control/picture_draw.R` | gridExtra、ggpubr |
| `presentation/exon_distance.R`、`gene_ontology.R` | ggplot2 |

`grid` 随 R 提供。绘图与区域注释脚本所用的输入表、坐标和注释版本须与分析结果一致。
