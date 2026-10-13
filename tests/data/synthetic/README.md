# 合成示例数据

本目录的序列和 ENST/ENSG ID 均为人工生成的示例，不对应真实人基因。生成器为 [make_fixture.py](../../make_fixture.py)，随机种子为 6096。

数据包含一个 500 nt 转录本、配套微型基因组和五列 exon 表、两条 240 nt 的 rRNA/细菌参考，以及一个 NC 和三个 TR 的无 UMI 双端 reads。每个样本含 2 对细菌 rRNA、3 对宿主 rRNA；mRNA 的预设末端计数见 [truth.json](truth.json)，候选位置为基因组坐标 219。

安装 [主 README](../../../README.md) 所列环境后，在仓库根目录运行：

```bash
python3 tests/smoke.py
```

命令使用临时目录建立索引、完成比对及评分，结束后删除生成的文件。该示例无须完整参考，可用于熟悉安装和运行方式；合成序列的结果不能用于生物学推断。
