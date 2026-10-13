# 软件依赖

主 README 使用 [environment.yml](../environment.yml) 建立环境，版本选择以本项目的测试结果为准。Conda 包来自 [conda-forge](https://conda-forge.org/) 和 [Bioconda](https://bioconda.github.io/recipes/bowtie2/README.html)。目前验证的平台是 Linux x86_64；Windows 请在 Linux 环境中运行。

Linux x86_64 上也可按已测试的完整包版本安装，包括间接依赖：

```bash
conda create -n nju-seq --file config/conda-linux-64.lock
conda activate nju-seq
```

这与 README 的环境创建方式二选一。[版本清单](../config/conda-linux-64.lock) 固定了包构建和校验值；其他平台应使用 `environment.yml`。

## 核心分析

`nju-seq` 由 Rust 标准库实现。源码使用 Rust 2024 edition 和 [Rust 1.88 引入的 let chains](https://blog.rust-lang.org/2025/06/26/Rust-1.88.0/)，因此最低编译版本为 1.88。

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

## UMI 文库与完整测试

在已激活的 `nju-seq` 环境中补装：

```bash
conda install --freeze-installed --override-channels -c conda-forge -c bioconda \
  umi_tools=1.1.6 perl-perlio-gzip=0.20 perl-scalar-util-numeric=0.40 \
  perl-app-cpanminus make
cpanm --mirror https://cpan.metacpan.org --mirror-only AlignDB::IntSpan@1.1.1

umi_tools --version
perl -MPerlIO::gzip -MAlignDB::IntSpan -e 'print "Perl modules OK\n"'
```

UMI-tools 用于 `r1_3p8` 标签提取及 Python 对照测试。Rust 中的 UMI 合并采用 exact/unique 规则。`AlignDB::IntSpan` 和 `PerlIO::gzip` 用于旧 Perl 计数、注释脚本及对照测试；独立 rRNA 的 `readend_count.pl` 只使用 Perl 核心模块。

`AlignDB::IntSpan` 的依赖 [Scalar::Util::Numeric](https://bioconda.github.io/recipes/perl-scalar-util-numeric/README.html) 使用 Conda 提供的预编译包，避免本地 Perl 扩展编译器配置影响安装。

默认使用 `PATH` 中的工具。多环境并存时，可用 `NJU_PERL=/path/to/perl` 指定 Perl 对照解释器，`NJU_NATIVE=/path/to/nju-seq` 指定测试二进制；这些路径需替换为实际位置。不要将另一套 Perl 的模块目录混入当前环境。

## 可选绘图

Rust 生成的 `qc.html` 不需要 R。运行旧绘图脚本时，按所用脚本安装 R 包：

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

`grid` 随 R 提供。绘图与区域注释还需相应的输入表，不能只安装软件就直接套用任意评分文件。
