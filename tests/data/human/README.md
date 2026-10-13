# 人源小份测试数据

NJU6096（NC）与 NJU6097、NJU6098、NJU6099（TR），各原始双端 FASTQ 的前 **1,000 对 reads**，均 `umi=none`。实验来源沿用根目录 [sample_list.csv](../../../sample_list.csv)：HEK293T rRNA-free RNA，NC without MgR、TR with MgR。

[样本清单](samples.tsv) 为相对路径，可随文件夹搬迁。[audit.json](audit.json) 记录完整输入的实际来源、大小/mtime、所取范围及输出 SHA-256。抽取时仅读取所需片段，未检查完整源 gzip 的 CRC。这是流程测试用连续子集，不能代表整库分布或预估候选数。

用根目录 README 的 raw 命令和完整人源参考运行。1,000 对的测试可能缺少 NC/TR 共同覆盖，使用保留参考得到 0 个联合候选与 `normalization_unavailable`，不视为程序失败；合成测试另外检查非空候选和预期计数。

再生成同一范围：

```bash
python3 scripts/subset_fastq.py --samples config/samples.human.tsv \
  --pairs 1000 --skip 0 --out output/recreated_human_test
```
