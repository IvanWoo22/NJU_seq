# 示例结果与输出说明

按 [主 README](../README.md) 完成人源分析后，可用本页解释计数、评分状态和输出文件。

## 人源小样本

[示例数据](../tests/data/human/README.md) 包含 NJU6096（NC）及 NJU6097–6099（TR），每个样本 1,000 对 reads，无 UMI。使用示例中的配套人源参考和软件版本，各阶段计数如下：

| 样本 | 排除的细菌 rRNA pairs | 排除的人 rRNA pairs | 通过筛选的 mRNA reads |
|---|---:|---:|---:|
| NJU6096 | 6 | 14 | 18 |
| NJU6097 | 81 | 41 | 13 |
| NJU6098 | 73 | 13 | 5 |
| NJU6099 | 78 | 29 | 2 |

每个通过 mRNA 筛选的 read 对应一条定位，read 分支按 R1 末端计数。rRNA 排除数按双端 pairs 计。

NJU6098、NJU6099 对应的 NC 归一化总量为 0，因此联合评分状态为 `normalization_unavailable`，候选数为 0。这表示小样本共同覆盖不足，不能据此判断完整文库是否存在 Nm 位点。

## 查看输出

| 文件 | 读取方式 |
|---|---|
| `SUCCESS` | 全部分析步骤成功结束 |
| `FAILED`、`logs/`、`stages.tsv` | 查看失败步骤、错误信息及退出状态 |
| `qc.html` | 浏览器离线查看质量报告；默认 `--storage compact` 生成 |
| `<sample>/reads/mrna.filter.tsv` | 九列基因组末端计数，无表头 |
| `mrna.reads.scored.tsv` | 默认只包含通过联合阈值的候选 |
| `mrna.reads.score_summary.json` | 查看归一化状态、参与评分的样本和候选数 |

评分表只有表头时，应同时查看摘要。归一化不可用、共同位点为空和完成评分但没有候选是不同情况，详见 [评分规则](scoring.md)。

## 含 UMI 的结果

声明 `umi=r1_3p8` 的样本还会生成 `<sample>/umi/mrna.filter.tsv`、UMI 家族表及合并统计。只有 NC 和全部 TR 均声明 UMI 时，才生成 `mrna.umi.scored.tsv` 与相应摘要。

read 计数以通过筛选的 R1 为单位；UMI 计数以每个样本内 exact 合并后的分子为单位，含 `N` 的标签逐条保留。拟南芥 Col_NC31 / Col_TR31 的样本结果及解释见 [UMI 分析流程](arabidopsis_umi.md#样本分析结果)。
