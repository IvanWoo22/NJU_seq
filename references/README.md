# 人源参考与原 data

`data/` 保存仓库自带的小型参考，`references/` 存放本地完整参考与索引，`inputs/` 存放测序 FASTQ 或其链接。完整参考和 FASTQ 不纳入 Git；新克隆需要另外复制这些文件。

## 当前参考

| 目录 | 内容 |
|---|---|
| `human/` | 已准备的人源 mRNA、rRNA、细菌 rRNA 索引和 `exon.info` |
| `human_source/` | 原人源转录本 FASTA、五列 exon 表、GRCh38 基因组及其他注释 |
| `bacterial/` | 原细菌 rRNA 的六片 Bowtie2 索引 |
| `arabidopsis/` | 保留的拟南芥参考 |

人源 mRNA 参考含 84,068 条 protein_coding 转录本，对应 762,910 行 exon。人 rRNA 使用原 `data/hsa_rrna/18s.fa`、`28s.fa`、`5-8s.fa`；它们只覆盖这三条参考。

原人源注释的确切发布号及细菌索引建库来源尚未确认。复现已有结果应使用保留的原文件，不能仅凭 `hg38` 文件名另选 GENCODE 版本。`human_source/` 中的 basic、代表转录本、区域和 codon 表供后续分析参考，不能替代当前全转录本与 exon 的配套输入。

`data/` 中其他物种的 rRNA 及 `MCMV.fa`、`MHV.fa` 均保留。换物种时需要同时准备匹配的转录本和注释。

## 重建与搬迁

按 [主 README](../README.md) 的命令运行 `scripts/prepare_human_reference.py`。输出为一个新目录，包含 `mrna`、`rrna`、`bac` 三个索引前缀、`exon.info`、FASTA 和 `audit.json`。

- 加 `--copy-inputs`：复制基因组和细菌索引，生成的目录可单独搬迁。
- 不加此选项：使用相对链接，应保持参考目录与源文件目录的相对位置。
- 使用其他宿主 rRNA：重复传入 `--rrna FILE`，每个 FASTA 分别读取。

脚本逐条校验转录本与基因组 exon 重建序列。失败时保留 `.partial` 目录供检查，修正输入后改用新输出目录。既有 `human/` 也使用相对链接，应连同 `human_source/` 和 `bacterial/` 一起复制；仅复制 `human/` 时需要解引用链接，例如 `cp -RL references/human /destination/human`。

## 输入格式

转录本 FASTA 使用 GENCODE 风格的 `|` 分隔标题：前两项为带版本号的 ENST、ENSG，第八项为 `protein_coding`，支持 `_PAR_Y`。

exon 表无表头，用实际制表符分隔五列：

```text
chromosome    start    end    strand    attributes
chr1          101      200    +         Parent=ENST000001.1;gene_id=ENSG000001.1
```

上面的空格仅用于展示。坐标为 1-based 闭区间；脚本仅选取输入 FASTA 中出现的转录本，并保留基因对应关系、负链和 PAR-Y 信息。九列 GFF3 需要先转换，不能直接作为 `--exons`。

基因组须为未压缩 FASTA，转录本和 rRNA FASTA 可用 gzip。细菌索引须是同一种格式的六个 `.bt2` 或 `.bt2l` 文件。

旧目录对应关系见 [迁移说明](../docs/migration.md)。已有 `human/audit.json` 中的旧路径是来源记录，运行时不读取这些路径。
