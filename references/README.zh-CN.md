# 参考文件生成

[English](README.md) | **简体中文**

按[安装说明](../README.zh-CN.md#1-安装与构建)准备环境后，在仓库根目录运行以下命令。脚本从公开来源下载文件、筛选序列、转换注释、校验转录本与基因组外显子序列，并建立 Bowtie2 索引。完整参考在本地生成；Git 提供脚本、来源校验值及小型序列修订和区段清单。

## 命令生成的文件

| 目录 | 生成内容 | 对应命令 |
|---|---|---|
| `downloads/` | 公开来源的压缩文件，按固定 SHA-256 校验 | 以下三个命令 |
| `bacterial/` | `bac_rna.fa`、六个 `bac_rna.*.bt2` 文件及 `audit.json` | `build_references.py bacterial` |
| `human_source/` | `hsa_protein_coding.fa`、`hsa_exon.info`、GRCh38 基因组及 `sources.json` | `build_references.py human` |
| `human_rebuilt/` | `mrna`、`rrna`、`bac` 索引前缀，`exon.info`、FASTA 及来源信息 | `build_references.py human` |
| `arabidopsis_source/` | 蛋白编码转录本 FASTA、外显子表、TAIR10 基因组及扩展 rRNA FASTA | `build_references.py arabidopsis` |
| `arabidopsis/` | `ath_protein_coding`、`ath_rrna_total`、`bac` 索引前缀，`exon.info`、FASTA 及来源信息 | `build_references.py arabidopsis` |

表中路径均相对于 `references/`。输出目录和源文件目录须为新目录；如果已经存在，用 `--out` 和 `--source-out` 指定其他路径，并在后续分析命令中使用对应路径。

## 1. 细菌 rRNA

```bash
python3 scripts/build_references.py bacterial \
  --download-dir references/downloads \
  --out references/bacterial --threads 4
```

下载 **SILVA 138.2 SSURef NR99 和 LSURef NR99 的 truncated FASTA**，按分类注释中的完整科名筛选 `Enterobacteriaceae`、`Lactobacillaceae` 和 `Streptococcaceae`。脚本将 RNA 的 `U` 转换为 DNA 的 `T`，保留原文件顺序（先 SSU、后 LSU），并检查序列 ID 是否重复。输出为 **34,982 条 SSU + 21,135 条 LSU，共 56,117 条序列**，属于指定细菌科的 rRNA 子集，并非整个 SILVA 数据库。

脚本校验筛选后 FASTA 的校验值，再调用 `bowtie2-build`。分析使用的索引前缀为 `references/bacterial/bac_rna`，人源和拟南芥流程均使用这一索引。

## 2. 人：GENCODE v34 / GRCh38

```bash
python3 scripts/build_references.py human \
  --download-dir references/downloads \
  --source-out references/human_source \
  --bac-index references/bacterial/bac_rna \
  --out references/human_rebuilt --threads 4 --copy-inputs
```

脚本下载 v34 转录本 FASTA、完整 GFF3 注释和 GRCh38 primary assembly 基因组。先筛选 FASTA 标题中以 `|` 分隔的第八项为 `protein_coding` 的记录，再筛选对应的 GFF3 `exon` 特征。GFF3 第 1、4、5、7 列和转录本/基因属性组成五列外显子表。此处筛选蛋白编码转录本，不进行代表转录本筛选。

生成的 mRNA 参考含 **84,068 条转录本和 762,910 行 exon**。每条转录本须与基因组外显子重建的序列一致，包括负链的反向互补和 `_PAR_Y` 标识。宿主 rRNA 使用仓库自带的 [18s.fa](../data/hsa_rrna/18s.fa)、[28s.fa](../data/hsa_rrna/28s.fa) 和 [5-8s.fa](../data/hsa_rrna/5-8s.fa)。

在 `nju-seq` 中使用 `--index references/human_rebuilt --exons references/human_rebuilt/exon.info`，完整人源示例见[主 README](../README.zh-CN.md#4-比对计数与评分)。

## 3. 拟南芥：Ensembl Plants 46 / TAIR10

```bash
python3 scripts/build_references.py arabidopsis \
  --download-dir references/downloads \
  --source-out references/arabidopsis_source \
  --bac-index references/bacterial/bac_rna \
  --out references/arabidopsis --threads 4 --copy-inputs
```

mRNA 参考使用 release 46 的 TAIR10 基因组、cDNA FASTA 和 GFF3。脚本筛选 `transcript_biotype:protein_coding` 记录、保留配套的外显子特征，并逐条核对 **48,321 条转录本**与基因组序列。外显子表保留 Rust 流程使用的 `Parent=transcript:AT…` 属性格式。

**4,131 条序列的 rRNA 参考**单独由以下内容生成：

- [data/ath_rrna](../data/ath_rrna) 中的四条 rRNA，按 [arabidopsis_rrna.json](../config/arabidopsis_rrna.json) 中列明的替换和插入修订，校验修订前后的序列 SHA-256。
- 公开 **Col-PEK 基因组 GWHBFWQ00000000.1** 中的 **4,127 个区段**，坐标见 [arabidopsis_rrna.intervals.bed](../config/arabidopsis_rrna.intervals.bed)。这些区段包含 rRNA 拷贝和片段，用于去除 rRNA reads。

BED 使用 **0-based、左闭右开**坐标，提取基因组正向序列；第四列为固定的 FASTA ID。Chr4 的 ID 保留 Col-PEK 命名中 100 bp 起始缺口的偏移，提取时应使用 BED 第 2–3 列坐标，不能从 ID 解析坐标。这些名称是 rRNA 参考标签，不是 TAIR10 的 mRNA 坐标。重建的完整 rRNA FASTA 须通过固定 SHA-256 校验，才会建立索引。

生成的索引前缀为 `ath_protein_coding` 和 `ath_rrna_total`。随后按[拟南芥 UMI 流程](../docs/arabidopsis_umi.zh-CN.md#先跑小样本)分析样本。计数命令使用生成的 `exon.info`，参考目录中的基因组文件名为 `genome.fa`。

## 下载、存储与校验

具体 URL 和 SHA-256 见 [reference_sources.json](../config/reference_sources.json)；[参考来源清单](../docs/reference_sources.zh-CN.md)还列有小鼠 GENCODE M25 和可选参考。下载文件保留原始文件名。若只下载、不建索引，可运行：

```bash
python3 scripts/build_references.py human \
  --download-dir references/downloads --download-only
```

需要时分别对 `bacterial` 和 `arabidopsis` 执行同样操作。已下载的文件通过校验后会直接复用。给构建命令添加 `--offline` 可要求所有源文件已在本地。使用其他方式下载的文件时，按 JSON 清单中的原始文件名放入下载目录，脚本仍会检查校验值。

`--copy-inputs` 将基因组和细菌索引复制进宿主参考目录，生成的参考目录可单独搬迁。不加此选项则创建相对链接，需要保留源文件目录和参考目录的相对位置。人源构建需数 GB 空间存放解压的基因组和索引；复制输入会额外占用空间。完整参考和生成的索引已由 `.gitignore` 排除。

参考目录包含来源及序列校验值、参考条目数和 Bowtie2 构建日志。构建失败时保留 `.partial` 目录供检查，修正问题后请使用新的输出目录和源文件目录；脚本不会覆盖已经生成的参考。下载失败会留下 `.partial` 文件，下次执行时重新下载该文件。

## 使用已准备的人源输入

底层脚本可以直接接受配套的源文件：

```bash
python3 scripts/prepare_human_reference.py \
  --mrna references/human_source/hsa_protein_coding.fa \
  --exons references/human_source/hsa_exon.info \
  --genome references/human_source/GRCh38.primary_assembly.genome.fa \
  --bac-index references/bacterial/bac_rna \
  --out references/human_custom --threads 4 \
  --label 'GENCODE v34 / GRCh38' --copy-inputs
```

转录本 FASTA 使用 GENCODE 风格的 `|` 分隔标题，前两项为带版本号的 ENST 和 ENSG，第八项为 `protein_coding`。外显子表无表头，用实际制表符分隔五列：

```text
chromosome    start    end    strand    attributes
chr1          101      200    +         Parent=ENST000001.1;gene_id=ENSG000001.1
```

上面的空格仅用于展示。外显子坐标为 **1-based 闭区间**。九列 GFF3 不能直接传给 `--exons`，`build_references.py` 会完成转换。基因组须解压，转录本和 rRNA FASTA 可用 gzip。细菌索引须包含同一种格式的六个 `.bt2` 或 `.bt2l` 文件。需要自定义人源 rRNA 时，重复传入 `--rrna FILE`。
