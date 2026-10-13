# 人源示例数据

[English](README.md) | **简体中文**

NJU6096（NC）与 NJU6097、NJU6098、NJU6099（TR），各取原始双端 FASTQ 的前 **1,000 对 reads**，均为 `umi=none`。实验信息见 [sample_list.csv](../../../sample_list.csv)：HEK293T rRNA-free RNA，NC 未经 MgR 处理，TR 经 MgR 处理。

[样本清单](samples.tsv) 使用相对路径，可随数据文件夹一起复制。[数据来源与校验值](audit.json) 包含所取范围及 FASTQ 的 SHA-256。这份连续子集用于演示分析步骤，不代表整库分布，不能用于估计整库候选数。

按 [主 README](../../../README.zh-CN.md) 的 `--from raw` 命令，使用配套人源参考运行。各步计数及 `normalization_unavailable` 状态的解释见 [示例结果](../../../docs/results.zh-CN.md)。

已有完整 FASTQ 时，可在仓库根目录抽取同一范围：

```bash
python3 scripts/subset_fastq.py --samples config/samples.human.tsv \
  --pairs 1000 --skip 0 --out output/human_subset
```
