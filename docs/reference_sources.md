# 人、小鼠、拟南芥参考来源

本页按原 README 整理人、小鼠和拟南芥的参考版本、下载地址及文件名。文件简称沿用 `hsa`、`mmu`、`ath`，便于对应分析脚本。下载的完整参考放在 `references/` 下，仓库自带的小型参考位于 `data/`。

| 物种 | 发布版本 | 基因组版本 | 官方入口 |
|---|---|---|---|
| 人（Homo sapiens） | GENCODE v34 | GRCh38.p13 | [GENCODE release 34](https://www.gencodegenes.org/human/release_34.html) |
| 小鼠（Mus musculus） | GENCODE M25 | GRCm38.p6 | [GENCODE release M25](https://www.gencodegenes.org/mouse/release_M25.html) |
| 拟南芥（Arabidopsis thaliana） | Ensembl Plants release 46 | TAIR10 | [release-46 文件目录](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/) |

## 人：GENCODE v34

| 内容 | 文件与下载地址 | 文件简称 |
|---|---|---|
| GFF3 注释 | [gencode.v34.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gff3.gz) | `hsa.gff3.gz` |
| 全部转录本 FASTA | [gencode.v34.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz) | `hsa_transcript.fa.gz` |

从全部转录本中筛选 `protein_coding` 后建立 mRNA 索引。人源参考准备脚本还需匹配的基因组 FASTA，可从同一发布目录获取 [GRCh38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz)。

## 小鼠：GENCODE M25

| 内容 | 文件与下载地址 | 文件简称 |
|---|---|---|
| GFF3 注释 | [gencode.vM25.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz) | `mmu.gff3.gz` |
| 全部转录本 FASTA | [gencode.vM25.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz) | `mmu_transcript.fa.gz` |

如需从基因组核对转录本，可补充同目录的 [GRCm38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz)。M25 对应 GRCm38，不能混用 GRCm39 注释或坐标。

Rust 流程目前支持人源 ENST/ENSG 和拟南芥 AT 注释格式；小鼠 ENSMUST/ENSMUSG 需扩展注释解析后使用。

## 拟南芥：Ensembl Plants release 46 / TAIR10

| 内容 | 文件与下载地址 | 文件简称 |
|---|---|---|
| 基因组 FASTA | [Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz) | `ath_dna.fa.gz` |
| GFF3 注释 | [Arabidopsis_thaliana.TAIR10.46.gff3.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz) | `ath.gff3.gz` |
| cDNA FASTA | [Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz) | `ath_transcript.fa.gz` |
| ncRNA FASTA | [Arabidopsis_thaliana.TAIR10.ncrna.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz) | `ath_ncrna.fa.gz` |

使用同一 release-46 目录下的基因组、转录本和注释文件，保持序列与坐标对应。

## 原 data 与补充参考

| 参考 | 文件位置或来源 | 用法 |
|---|---|---|
| 人 rRNA | [data/hsa_rrna](../data/hsa_rrna)：18s、28s、5-8s | 当前人源示例使用这三条序列 |
| 小鼠 rRNA | [data/mmu_rrna](../data/mmu_rrna)：18s、28s、5-8s | 小鼠 rRNA 参考序列 |
| 拟南芥 rRNA | [data/ath_rrna](../data/ath_rrna)：18s、25s、5-8s、5s | 拟南芥 rRNA 参考序列；计数所用 FASTA 从实际索引导出，见 [UMI 流程](arabidopsis_umi.md#可选-rrna-分析) |
| 细菌 rRNA | `references/bacterial/bac_rna.*.bt2` | 比对前准备齐全的六片索引 |
| miRNA（可选） | [miRBase CURRENT/mature.fa.gz](ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz)，文件简称 `mirna.fa.gz` | `CURRENT` 随发布更新，使用时保存实际版本和文件校验值 |

转录本、GFF3 和基因组需使用配套版本，下载后保存文件的 SHA-256 校验值。目录布局与 exon 格式见 [参考说明](../references/README.md)。
