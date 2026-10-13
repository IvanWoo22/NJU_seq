# 评分格式与规则

`nju-seq score` 将一个 NC 与全部 TR 联合计算。mRNA 与 rRNA 共用输出格式，但输入计数、归一化范围和候选坐标不同。

## 运行选项

在仓库根目录构建：

```bash
cargo build --release --offline
target/release/nju-seq --version
```

通用接口：

```text
nju-seq score --mode mrna|rrna --nc NC.tsv --tr TR1.tsv [--tr TR2.tsv ...]
  --output scores.tsv [--summary scores.json]
  [--format unified|legacy] [--rows candidates|all]
  [--reference-id 18s]
```

`--mode` 默认为 `mrna`，`--format` 默认为 `unified`。rRNA 必须指定不含空白的 `--reference-id`，例如 `18s`、`28s`、`5-8s`。模式包含输入解析及评分规则，不能将 rRNA 四列输入直接交给 mRNA 模式比较。

样本名取自输入文件的父目录。因此 NC 和每个 TR 应分别位于名称不同的样本目录，如 `NJU6096/rrna_18s.tsv`。统一格式拒绝重复列名。输出及可选摘要均不得已存在；输入、选项或计算失败时不写出结果。

| 选项 | mRNA | rRNA |
|---|---|---|
| 输入 | 无表头九列基因组末端计数，末端计数为正整数 | 无表头四列 `Position/Base/StartCount/EndCount`；位置从 1 连续至 L，计数可为 0 |
| 默认输出范围 | `candidates`，仅联合候选 | `all`，候选位置 1…L−1 |
| `--rows all` | 所有 TR 都有正计数的交集；不是全基因组所有位点 | 位置 1…L−1，包含未通过阈值和无法评价的行 |
| `--rows candidates` | 仅 `Candidate=true` 的行 | 仅 `Candidate=true` 的行 |

rRNA 输入必须具有一致的长度和逐位点碱基；空参考、断号、重复位置和负计数均报错。长度为 1 时仅输出表头；长度不超过 40 时没有通过边界屏蔽的内部位置。

建议总是提供 `--summary`：候选表为空时，摘要仍能区分归一化不可用、交集为空和正常评价后没有候选。摘要还记录版本、模式、输出范围、规则、每个 TR 的 NC/TR 总量和归一化因子。`joint_evaluation=available` 表示归一化可用，不代表一定存在可评分的内部位点。

`nju-seq run` 自动保存 read 及可选 UMI 的评分摘要。rRNA 计数与评分用法见 [rRNA 分支](rrna.md)。

## 计算公式

令 `signal` 为归属于候选的原始末端计数，`adjacent` 为局部比值的相邻原始计数，`n` 为 TR 数量。两种模式共同使用：

```text
score_j       = TR_signal_j / max(TR_adjacent_j, 1)
                − NC_score_numerator / max(NC_adjacent, 1)
corrected_j   = ceil(TR_signal_j × NC_TOTAL_j / TR_TOTAL_j)
TotalScore    = sum(score_j)
AveScore      = TotalScore / n
TotalCorrectedCount = sum(corrected_j)
TotalFoldChange     = TotalCorrectedCount / max(NC_signal, 1)
AveFC               = TotalFoldChange / n

Candidate = 可评价且不在屏蔽区
            且 TotalScore >= 30 × n
            且 TotalCorrectedCount >= 3 × n × NC_signal
            且 TotalCorrectedCount > 0
```

`corrected_j` 按先乘、后除、最后向上取整的顺序计算，使用未经舍入的归一化因子。联合阈值按所有 TR 的得分总和判定，单个 TR 的 score 可以低于 30。

**NC 单个位点为 0 不等于 NC 归一化总体为 0。** 当位点为 0 而总体有效时，FC 门槛仍按原始 NC=0 并要求校正和为正；此时 `FC_criteria=meet` 不必等价于 `AveFC>=3`。不能用显示的平均 FC 独自重建候选判定。

| 规则 | `--mode mrna` | `--mode rrna` |
|---|---|---|
| NC_TOTAL | 分别汇总各 TR 检出位点上的 NC 计数 | 当前 rRNA 参考全部位置的 NC 计数 |
| TR_TOTAL | 该 TR 输入中的全部末端计数 | 当前 rRNA 参考全部位置的 TR 计数 |
| 可评分支持集合 | 所有 TR 正计数位点的交集 | 完整参考；允许某个 TR 在单个位点计数为 0 |
| NC 信号为 0 时的 score 分子 | 1 | 0 |
| 相邻计数 | 正链候选基因组坐标 −1，负链 +1 | 候选参考位置 k；分子为位置 k+1 |
| 首尾屏蔽 | 不额外屏蔽 | 候选 k≤20 或 k>L−20 排除；最后一位没有 k+1，不输出 |

mRNA 在外显子连接处仍使用基因组相邻计数，不能将其解释为转录本上相邻位置的计数。

## 无法评分的情况

每个 TR 均须满足 `NC_TOTAL>0` 且 `TR_TOTAL>0`。任一不可用时，整个联合评分标记为 `normalization_unavailable`，不自动删除该 TR 或减少 n，也不允许其余 TR 的强信号补偿一个完全无计数的文库。

| `EvaluationStatus` | 含义 | 分数、校正计数、总/平均 FC | `Candidate` / 门槛字段 |
|---|---|---|---|
| `evaluated` | 归一化可用，且通过位置边界检查 | 实际数值，可为 0 | true/false；meet/unmeet |
| `normalization_unavailable` | 至少一个 TR 所需总体为 0 | `NA` | false；not_evaluated |
| `boundary_excluded` | rRNA 首尾屏蔽区 | `NA` | false；not_evaluated |

归一化不可用优先于边界状态。输出保留原始 signal/adjacent 和总体计数，便于解释原因；各归一化因子单独显示，只有对应总体无效时才为 `NA`。绘图时应按状态区分不可评价位置与有效零分。

## 统一输出字段含义

两种模式共用表头；一个 NC、三个 TR 时均为 46 列，包含坐标、注释、计数、评分、归一化量及评价状态。

| 字段 | 统一定义 |
|---|---|
| `Chromosome / Position / Direction` | 候选位置，1-based；mRNA 为基因组坐标，rRNA 为显式 reference-id 内坐标、方向 `+` |
| `CoordinateSystem` | `genome_1based` 或 `rrna_reference_1based`，防止将 rRNA 坐标当成 hg38 坐标 |
| `LeftBase / Base / RightBase / GeneID` | mRNA 使用输入注释；rRNA 为参考相邻碱基，GeneID 为 reference-id |
| `<NC>_count(NC)`、`<TR>_count` | 用于评分的原始 **signal** 计数；均归属于显示的候选位置 |
| `<sample>_adjacent_count` | 局部比值分母位置的原始计数 |
| `<NC>_score_numerator` | 实际用于 NC 局部比值的分子，包含该模式对零值的补 1 处理 |
| `<sample>_score_denominator` | 实际分母，即 `max(adjacent_count,1)` |
| `<TR>_count_corr / <TR>_score` | 单个 TR 的归一化计数及与 NC 的局部比值之差 |
| `<TR>_NC_TOTAL / TR_TOTAL / normalization_factor` | 该 TR 所用总体及 `NC_TOTAL/TR_TOTAL` |
| `TotalFoldChange / AveFC` | TR 校正倍数之和 / 其算术平均值 |
| `TotalScore / AveScore` | TR score 之和 / 其算术平均值 |
| `TotalCorrectedCount` | 所有 TR 归一化计数之和 |
| `ScoreMode / EvaluationStatus / Candidate` | 规则模式、评价状态、最终联合候选判断 |
| `FC_criteria / Score_criteria` | 两个联合阈值各自的结果；无法评价时为 not_evaluated |

候选位置是沿转录方向位于物理 read 末端前一位的碱基。在 `unified` 格式中，rRNA 行 k 的 `*_count` 为 **位置 k+1 的末端计数**，`*_adjacent_count` 为位置 k 的末端计数。这两个字段分别对应局部比值的分子和分母。

### legacy 格式

下游脚本需要 Perl `score.pl` 的表头和列序时，可指定 `--format legacy`。mRNA 仅输出候选表，不支持 `--rows all`。rRNA 使用 `#Position` 列，count 表示该行位置 k 的物理末端计数，数值保留两位小数；`FoldChange` 对应 `unified` 格式的 `TotalFoldChange`。

两种格式使用相同的归一化条件。`legacy` 的 rRNA 表将无法评价的位置表示为 0/unmeet；需要区分归一化不可用和边界屏蔽时，使用 `unified` 格式或 `--summary` 摘要。选择绘图与注释脚本时，应核对其所需的列名和坐标定义。
