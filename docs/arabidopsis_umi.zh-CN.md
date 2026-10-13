# 拟南芥 Col_NC31 / Col_TR31：8 nt UMI 分析

[English](arabidopsis_umi.md) | **简体中文**

Col_NC31 为 NC 对照，Col_TR31 为 TR 处理样本。两份文库在 **MgR 处理后、PCR 前**加入 8 nt UMI，由 R1 的 3′ 端读取；R2 的 5′ 端为同一标签的反向互补方向。本文介绍 Rust 版 `nju-seq` 的运行方法、UMI 计数规则及两份样本的分析结果。

[软件与输入](#软件与输入) · [参考](#使用配套的拟南芥参考) · [小样本](#先跑小样本) · [UMI 处理](#umi-在流程中的位置) · [结果](#样本分析结果) · [rRNA](#可选-rrna-分析)

## 软件与输入

在仓库根目录操作。先按 [主 README](../README.zh-CN.md#1-安装与构建) 建立 `nju-seq` 环境，再补装 UMI-tools：

```bash
conda activate nju-seq
conda install --freeze-installed --override-channels -c conda-forge -c bioconda umi_tools=1.1.6
cargo build --release --offline
umi_tools --version
```

UMI 标签由 UMI-tools 提取，分子合并和计数由 Rust 完成。软件依赖见 [依赖说明](dependencies.zh-CN.md)。

将原始双端 FASTQ 放在以下位置，也可建立链接或修改 [样本清单](../config/samples.arabidopsis_umi.tsv) 中的路径：

```text
inputs/arabidopsis/Col_NC31/R1.fq.gz
inputs/arabidopsis/Col_NC31/R2.fq.gz
inputs/arabidopsis/Col_TR31/R1.fq.gz
inputs/arabidopsis/Col_TR31/R2.fq.gz
```

两行清单的 `umi` 均为 `r1_3p8`，`alignment` 为 `-`。相对路径以清单所在目录为基准。不要将已去接头、已提取 UMI 的 FASTQ 作为 `--from raw` 输入；这类文件应使用 `--from prepared`，两端 read ID 都须保留提取后的 `_UMI` 后缀。

## 使用配套的拟南芥参考

从公开下载文件生成配套参考：

```bash
python3 scripts/build_references.py bacterial \
  --download-dir references/downloads \
  --out references/bacterial --threads 4

python3 scripts/build_references.py arabidopsis \
  --download-dir references/downloads \
  --source-out references/arabidopsis_source \
  --bac-index references/bacterial/bac_rna \
  --out references/arabidopsis --threads 4 --copy-inputs
```

如果已经生成 `references/bacterial/`，跳过第一个命令。输出目录和源文件目录须为新目录，需要时用 `--out` 和 `--source-out` 指定其他路径。

若不筛除细菌 rRNA，跳过第一个命令，将拟南芥构建命令中的 `--bac-index references/bacterial/bac_rna` 换成 `--without-bacterial`。此时不需要下载细菌参考或准备细菌索引。

| 路径或索引前缀 | 用途 |
|---|---|
| `references/bacterial/bac_rna` | SILVA 138.2 细菌 rRNA 子集，56,117 条序列；仅启用细菌筛除时需要 |
| `references/arabidopsis/ath_rrna_total` | 拟南芥扩展 rRNA 索引，4,131 条序列 |
| `references/arabidopsis/ath_protein_coding` | Ensembl Plants 46 / TAIR10 mRNA 索引，48,321 条转录本 |
| `references/arabidopsis/exon.info` | 配套的五列外显子表 |

mRNA 构建使用 release 46 的转录本、GFF3 和 TAIR10 基因组。扩展 rRNA 由四条指定 rRNA 序列及 4,127 个 Col-PEK 基因组区段组成。脚本会下载 Col-PEK 基因组，按随仓库提供的区段清单提取序列，并校验完整 rRNA FASTA 后建立索引。mRNA 计数仍使用 TAIR10 坐标。

每个索引前缀对应六个 `.bt2` 文件。`--copy-inputs` 使宿主参考目录可独立搬迁，其中的基因组文件名为 `genome.fa`。来源 URL、区段清单、序列修订及离线用法见[参考文件生成](../references/README.zh-CN.md#3-拟南芥ensembl-plants-46--tair10)和[参考来源清单](reference_sources.zh-CN.md)。

## 先跑小样本

每个样本先取前 1,000 对 reads。脚本保留清单中的 UMI 类型，并检查截取部分的双端 ID：

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

`plan` 只检查路径并列出命令，不执行分析。抽样和运行的输出目录均须尚不存在。小样本用于检查安装、配对处理与输出格式，不能用来估计整库候选数；覆盖不足时可能出现 `normalization_unavailable`。

以上命令默认依次筛除细菌和宿主 rRNA（`--rrna-filter both`）。仅筛除拟南芥 rRNA 时，在 `plan` 和 `run` 中都加上 `--rrna-filter host`，并省略 `--bac-index references/bacterial/bac_rna`。这样会跳过细菌比对与筛除，保留 UMI 提取和计数。全部四种取值见 [rRNA 筛除设置](../README.zh-CN.md#rrna-筛除设置)。

全量分析使用相同索引参数，将 `--samples` 改为 `config/samples.arabidopsis_umi.tsv`，`--out` 改为新的 `output/arabidopsis_umi_full`。`--threads` 是每个比对任务的线程数，`--jobs` 是并行样本数。默认成功后删除中间 FASTQ；需要保留时加 `--storage all`，同时需要 QC 报告时再加 `--qc yes`。分批比对与资源限制见 [辅助脚本](scripts.zh-CN.md)。

## UMI 在流程中的位置

`umi=r1_3p8` 会启用以下处理，两份文库分别运行：

1. Cutadapt 使用 R1 接头 `AGATCGGAAGAGCACA`、R2 接头 `GATCGTCGGACTGTAG`，以 `-U 8` 裁去 R2 起始标签；其余参数为 `-O 6 -m 18:10 -e 0.1 --discard-untrimmed`。R1 的最低长度包括尚未提取的 8 nt UMI。
2. UMI-tools 以 `--extract-method=string --bc-pattern=NNNNNNNN --3prime` 提取 R1 末端 8 nt，将标签追加到两端 read ID。提取标签本身不合并分子，R2 标签不参与联合纠错。
3. 默认依次排除细菌 rRNA 和拟南芥 rRNA，再比对 mRNA。使用 `--rrna-filter host` 时，仅筛除拟南芥 rRNA。启用的 rRNA 步骤使用 `-k 5`，mRNA 使用 `-a`；均保留 `--reorder --seed 0`。
4. 完成 mRNA 定位筛选后，保留 read 计数分支，并按 **染色体、链、R1 完整基因组外显子路径和 UMI** 做 exact 合并。不同起止点或剪接路径分别保留；等价转录本中的同一片段可合并。
5. 分别生成 read 与 UMI 末端计数，合并到基因组坐标，再进行 NC/TR 联合评分。只有清单中全部样本都声明 UMI 时，才生成联合 UMI 评分。

每个 UMI 家族选取 CIGAR 错误数最少的 read，平局时按 read ID 选择。含 `N` 的 UMI 逐条保留，避免把未知标签误合并。当前 Rust 版使用 exact 规则，不进行 directional 纠错；也不能只按 UMI 字符串或末端坐标全局去重。

| 输出 | 含义 |
|---|---|
| `<sample>/mrna.almostuniquematch.filter.tmp` | 通过 mRNA 筛选的 read 定位 |
| `<sample>/mrna.umi.filter.tmp` | 每个 UMI 家族的一条代表定位 |
| `<sample>/mrna.umi.families.tsv.gz`、`mrna.umi.summary.json` | 家族组成与合并统计，均位于样本目录内 |
| `<sample>/reads/mrna.filter.tsv` | read 计数表 |
| `<sample>/umi/mrna.filter.tsv` | UMI 计数表 |
| `mrna.reads.scored.tsv`、`mrna.umi.scored.tsv` | 两种计数下的评分结果，按 `Candidate` 列识别候选 |
| `mrna.reads.score_summary.json`、`mrna.umi.score_summary.json` | 归一化状态及候选数 |
| `qc.html`、`SUCCESS` | 质量报告、全部步骤成功标记 |

## 样本分析结果

以下为 Col_NC31 / Col_TR31 使用 Perl/Python 流程、按 exact UMI 规则计算的全量结果。

| 指标 | Col_NC31 | Col_TR31 |
|---|---:|---:|
| 原始 read pairs | 89,597,181 | 38,795,779 |
| 去接头后 read pairs | 31,663,953 | 12,546,074 |
| UMI 提取后 read pairs | 31,663,953 | 12,546,074 |
| 通过 mRNA 筛选的 reads | 3,239,437 | 575,214 |
| exact UMI 合并后的计数记录 | 874,050 | 395,251 |
| 其中含 `N` UMI、单独保留的记录 | 545 | 9,284 |

read 分支得到 134 个候选，exact UMI 分支得到 110 个；共有 107 个，read 分支独有 27 个，UMI 分支独有 3 个。候选变化反映计数方式的影响，不能直接解释为准确率提升。样本各只有一份 NC、TR；评分阈值不等同于重复间统计检验或 FDR，详见 [评分说明](scoring.zh-CN.md)。

## 可选 rRNA 分析

主流程排除 rRNA 是为了分析 mRNA，不会生成 rRNA 的 UMI 分子计数。需要 rRNA 位点结果时，按 [rRNA 独立计数与评分](rrna.zh-CN.md) 操作，并将索引和 FASTQ 换为本批拟南芥输入。

rRNA 计数使用与比对索引一致的参考序列。用 `bowtie2-inspect references/arabidopsis/ath_rrna_total` 导出序列，按 `18s`、`25s`、`5-8s`、`5s` 的参考 ID 分别保存为单条 FASTA，再向 `scripts/count_rrna.py` 重复传入 `--reference FILE`。FASTA 中保留原参考 ID，以便对应比对结果和计数坐标。
