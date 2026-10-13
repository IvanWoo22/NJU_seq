# 合成数据测试

这里的序列和 ENST/ENSG ID **均为生成的测试数据**，不是真实人基因。固定随机种子为 6096，生成器为 [make_fixture.py](../../make_fixture.py)。

包含一个 500 nt 转录本、匹配的微型基因组/五列 exon 表、两条 240 nt 的 rRNA/细菌参考，以及 NC + 3 TR 的无 UMI 原始双端 reads。每样本含 2 对 bac、3 对宿主 rRNA；mRNA 的已知端点计数见 [truth.json](truth.json)。基因组候选位置 219 应通过联合阈值。

在仓库根目录执行 `python3 tests/smoke.py`。测试临时建索引、重建并搬迁参考目录，分别运行 compact/all 两种存储方式，核对计数、候选、rRNA 删除、QC 与清理结果；结束自动删除生成物。使用主 README 的软件环境即可运行。

这些合成序列用于检查程序的计数和评分，不能用于生物学推断。正负链、PAR-Y、零覆盖、UMI 和失败传播另由 Rust 与 Python 回归测试覆盖。
