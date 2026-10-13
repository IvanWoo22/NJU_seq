# 测试与预期结果

先按 [主 README](../README.md) 安装环境、构建程序。合成数据测试只需要基本环境；其余 Python 测试还需要 [UMI-tools 和 Perl 模块](dependencies.md#umi-文库与完整测试)。所有命令均在仓库根目录运行。

```bash
cargo test --offline
cargo fmt --check
cargo clippy --offline --all-targets -- -D warnings
python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/umi_parity.py
python3 tests/smoke.py
```

测试使用临时目录，结束后自动清理。`NJU_NATIVE` 可指定另一个已构建的 `nju-seq`，`NJU_PERL` 可指定装有上述模块的 Perl。

## 人源小样本

使用保留的原人源参考，四个样本各输入 1,000 对 reads，预期结果如下：

| 样本 | 删除的细菌 rRNA pairs | 删除的人 rRNA pairs | 通过筛选的 mRNA reads |
|---|---:|---:|---:|
| NJU6096 | 6 | 14 | 18 |
| NJU6097 | 81 | 41 | 13 |
| NJU6098 | 73 | 13 | 5 |
| NJU6099 | 78 | 29 | 2 |

流程应生成 `SUCCESS` 和 `qc.html`。NJU6098、NJU6099 对应的 NC 归一化总量为 0，联合评分状态为 `normalization_unavailable`，候选为 0。小样本缺乏共同覆盖，不能据此判断完整数据是否存在 Nm 位点。合成数据测试另外检查候选位置 219 和已知末端计数。

## 最近一次检查：2026-10-13

在同一 Linux x86_64 主机上，新建 Conda 环境，并将源码和所需参考复制到一个含空格的新目录。清空继承的环境变量，重新编译程序，按 README 重建参考并运行真实小样本；未使用原 `target/`、旧日期目录或原分析环境的软件。

| 检查 | 结果 |
|---|---|
| Rust 1.88 构建、59 项单元测试、fmt、Clippy | 通过 |
| Python 测试 | 94 项通过 |
| exact UMI 对照 | 5,000 reads → 1,976 个分子；Rust/Python 定位表、family 表及统计一致 |
| 合成数据完整流程 | 计数、非空候选、QC、两种存储方式均符合预期 |
| 参考搬迁 | 相对链接连同源目录搬迁、复制包脱离源目录，两种方式均通过 |
| 完整人源参考重建 | 84,068 条转录本、762,910 行 exon；基因组重建序列不一致数为 0 |
| 四个人源小样本 | 全部成功，计数与上表一致 |
| 从原 FASTQ 重新截取小样本 | 八个 gzip 文件逐字节相同 |
| 原 `data/` 和 `sample_list.csv` | 29 个文件校验值不变 |

核心环境的完整版本见 [Conda 版本清单](../config/conda-linux-64.lock)。完整回归测试另安装 UMI-tools 1.1.6、PerlIO::gzip 0.20、AlignDB::IntSpan 1.1.1 等可选依赖；这些 Perl 包使用 Perl 5.32.1。简明结果与命令保存在 [reproduction.json](validation/reproduction.json)，各样本步骤结果见 [human_1k.json](validation/human_1k.json)。

已检查当前程序、辅助脚本和测试中的旧路径与脚本调用。旧机器路径仅保留在历史脚本、迁移表和来源记录中。完整测序数据未重跑，R 绘图未包含在此次测试中；新机器仍需另行提供完整参考。
