# 人、小鼠、拟南芥参考来源

以下沿用重构前 README（提交 `bd1c685`）的物种、发布版本和文件，不改为最新版本。旧文档以 `hsa`、`mmu`、`ath` 命名下载文件；表中保留这些名称，便于对应原脚本。下载的大文件应放在本地 `references/` 下，仓库自带的 `data/` 继续保留。

2026-10-13 已检查下列八个原始文件及两个补充基因组文件的 HTTPS 地址，均返回 HTTP 200；此次仅检查响应头，未下载完整参考。

| 物种 | 原 README 版本 | 基因组版本 | 官方入口 |
|---|---|---|---|
| 人（Homo sapiens） | GENCODE v34 | GRCh38.p13 | [GENCODE release 34](https://www.gencodegenes.org/human/release_34.html) |
| 小鼠（Mus musculus） | GENCODE M25 | GRCm38.p6 | [GENCODE release M25](https://www.gencodegenes.org/mouse/release_M25.html) |
| 拟南芥（Arabidopsis thaliana） | Ensembl Plants release 46 | TAIR10 | [release-46 文件目录](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/) |

## 人：GENCODE v34

| 内容 | 原文件与下载地址 | 旧本地名称 |
|---|---|---|
| GFF3 注释 | [gencode.v34.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gff3.gz) | `hsa.gff3.gz` |
| 全部转录本 FASTA | [gencode.v34.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz) | `hsa_transcript.fa.gz` |

原流程从全部转录本中筛选 `protein_coding` 后建立 mRNA 索引。当前人源参考准备脚本还需匹配的基因组 FASTA，可从同一发布目录获取 [GRCh38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz)；这是对旧下载清单的补充。

## 小鼠：GENCODE M25

| 内容 | 原文件与下载地址 | 旧本地名称 |
|---|---|---|
| GFF3 注释 | [gencode.vM25.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz) | `mmu.gff3.gz` |
| 全部转录本 FASTA | [gencode.vM25.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz) | `mmu_transcript.fa.gz` |

如需从基因组核对转录本，可补充同目录的 [GRCm38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz)。M25 对应 GRCm38，不能混用 GRCm39 注释或坐标。

本表提供原小鼠参考来源。当前 Rust 流程自动解析的注释格式为人源 ENST/ENSG 与拟南芥 AT 转录本；小鼠 ENSMUST/ENSMUSG 尚未完成适配和验证，不能仅替换索引就直接套用人源示例。

## 拟南芥：Ensembl Plants release 46 / TAIR10

| 内容 | 原文件与下载地址 | 旧本地名称 |
|---|---|---|
| 基因组 FASTA | [Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz) | `ath_dna.fa.gz` |
| GFF3 注释 | [Arabidopsis_thaliana.TAIR10.46.gff3.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz) | `ath.gff3.gz` |
| cDNA FASTA | [Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz) | `ath_transcript.fa.gz` |
| ncRNA FASTA | [Arabidopsis_thaliana.TAIR10.ncrna.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-46/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz) | `ath_ncrna.fa.gz` |

旧 README 使用 `ftp://ftp.ensemblgenomes.org/pub/plants/release-46/`；上表使用 EMBL-EBI 的 HTTPS 地址，发布目录和文件名保持不变。TAIR10 是基因组组装名称，不能代替注释发布号；不同 Ensembl Plants 版本不能据此互换。

## 原 data 与补充参考

| 参考 | 保留位置或来源 | 用法 |
|---|---|---|
| 人 rRNA | [data/hsa_rrna](../data/hsa_rrna)：18s、28s、5-8s | 当前人源示例使用这三条序列 |
| 小鼠 rRNA | [data/mmu_rrna](../data/mmu_rrna)：18s、28s、5-8s | 原参考文件继续保留 |
| 拟南芥 rRNA | [data/ath_rrna](../data/ath_rrna)：18s、25s、5-8s、5s | 与保留的 `ath_rrna_total` 索引并不完全相同，见 [UMI 流程](arabidopsis_umi.md#可选-rrna-分析) |
| 细菌 rRNA | `references/bacterial/bac_rna.*.bt2` | 复用原六片索引；旧 README 未给出可确认的建库来源 |
| miRNA（可选） | 原地址：[miRBase CURRENT/mature.fa.gz](ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz)，旧名 `mirna.fa.gz` | 原清单未固定发布号；如使用，应另记实际版本和文件校验值 |

这些来源是旧 README 的下载清单，不是对当前保留文件发布号的认定。原人源和拟南芥索引仍缺少确切的下载版本记录；要重现已有计数，应使用配套保留文件。另行下载时，将转录本、GFF3、基因组及其 SHA-256 一并保存，避免同名文件来自不同版本。当前目录布局与 exon 格式见 [参考说明](../references/README.md)。
