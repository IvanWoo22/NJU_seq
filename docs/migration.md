# 旧脚本与当前用法

## 目录

```text
Cargo.toml / Cargo.lock / src/ / assets/   Rust 主程序
scripts/                                 准备数据与运行辅助脚本
config/samples.human.tsv                  原始人源数据清单
inputs/human/                            原 FASTQ 链接（本地）
data/                                    原始小型参考（保留）
references/                              配套完整参考（本地）
tests/data/human/                         4 × 1,000 对真实 reads
tests/data/synthetic/                     合成数据、预期计数与微型参考
tests/                                   自动清理输出的回归测试
output/ / target/                        生成物，Git 忽略
docs/history/                            清理说明；详细历史归档仅在本地保留
```

## 旧参考目录

| 旧位置 | 当前位置 |
|---|---|
| `test_1001/references/hg38/` | `references/human/` |
| `test_1001/index/` | `references/human_source/` |
| `test_0924/index/bac_rna.*.bt2` | `references/bacterial/bac_rna.*.bt2` |
| `test_0924/index/` 其余参考 | `references/arabidopsis/` |

## 替换关系

| 旧操作 | 推荐入口 | 说明 |
|---|---|---|
| `test_rust_0926` crate | 仓库根 `cargo build --release --offline`；`target/release/nju-seq` | 版本 0.5.0，crate 名 `nju_seq` |
| 日期目录 `run_hg38.sh`、HPCC 固定路径 shell、旧 Nextflow | `nju-seq plan/run --samples ... --index ... --exons ... --out ...` | 本地执行，外层调度由调用者管理 |
| 固定 NJU6096–6099 的抽样代码 | `scripts/subset_fastq.py --samples ... --pairs ... --out ...` | 任意 NC/TR 清单；显式 UMI 设计 |
| 固定工作区的 `prepare_hg38.py references` | `scripts/prepare_human_reference.py` | 支持已验证的人源 GENCODE 风格输入；原 data/hsa_rrna 为默认 rRNA |
| 任意行数切割 SAM、GNU parallel 管道 | `nju-seq align-select`；`scripts/align_chunked.py` | 保留完整 read 组和 mRNA `-a` |
| SAM awk、质量/多重比对过滤、定位 dedup | Rust `filter` 或组合命令 `align-select` | 人源 ENST/PAR-Y 与拟南芥 AT 已验证 |
| 旧 `almostunique.sh` | Rust `almostunique`；主 run 已集成 | Rust 对所有 read 核对 R1；旧 shell 的单定位绕过核对政策只保留作历史对照 |
| 双端删除 read ID 的旧 shell | Rust `delete-pairs` | 同时核对两端，完成后写出结果 |
| mRNA `count.pl → sort → merge.pl → score.pl` | Rust `count → merge → score`；主 run 已集成 | 坐标/评分政策保持，默认评分表为 unified |
| 固定三种 rRNA 的 parallel 计数 | `scripts/count_rrna.py` | 包装既有 Perl 计数器；可重复指定单条 FASTA；之后 Rust score |
| R 单一 QC 汇总 | Rust `qc` / `qc-report` | 生成离线 HTML 和 TSV |
| 旧固定路径资源守护 | `scripts/guard_run.py` | Linux `/proc` / 进程组；默认上限适合小测试，完整数据需显式配置 |

`mrna_analysis/`、`rrna_analysis/`、`presentation/`、`quality_control/`、`tool/`、`lib/` 中的 Perl/R 脚本用于独立核对和后续分析。`log/` 中的旧机器路径及 Nextflow 配置仅供历史查阅。Rust 流程保持原比对、计数和评分规则；多重比对预筛、短序列剔除和错配救援未启用为默认选项。

`--index DIR` 默认前缀从 `bac_rna / ath_rrna_total / ath_protein_coding` 改为 `bac / rrna / mrna`。已有自定义/拟南芥包使用三个显式 `--*-index` 参数。

## 分批 mRNA 比对

先准备已经过接头处理、完成可选 UMI 提取及 rRNA 排除的双端 FASTQ。以下以一个样本为例，参数指向自己的输入：

```bash
mkdir -p output/batch_sample
python3 scripts/align_chunked.py \
  --native "$PWD/target/release/nju-seq" \
  --index references/human_rebuilt/mrna --exons references/human_rebuilt/exon.info \
  --r1 /path/to/sample.after_rrna.R1.fq.gz \
  --r2 /path/to/sample.after_rrna.R2.fq.gz \
  --output output/batch_sample/accepted.tsv \
  --log output/batch_sample/bowtie2.log \
  --stats output/batch_sample/selection.json --threads 4 --chunk-pairs 20000
```

`chunk-pairs` 按完整配对记录切分，每批调用同一个 Rust `align-select`，不是对 SAM 任意切行。没有将 `-a` 换成有限 `-k`；大量多重匹配仍可能消耗较多资源。多个样本完成后写一份 `--from unique` 清单（R1 指回实际输入、alignment 指向 accepted TSV），交给 Rust 统一计数评分。

断点通过 `scripts/chunk_checkpoint.py --help` 的 export 接口导出已提交批次，再传给 `align_chunked.py --resume-checkpoint DIR`。恢复时核对输入、参考、Rust/Bowtie2 文件和已完成批次的内容；只接受由当前脚本生成的断点。导出前停止生产者，不复用未提交尾部。该功能仅针对 mRNA 分批步骤，不能恢复整个 raw 工作流。

## 小测试资源保护

```bash
python3 scripts/guard_run.py --output output/human_guarded \
  --max-gib 4 --max-file-gib 1 --min-free-gib 20 --timeout 1800 -- \
  target/release/nju-seq run --samples tests/data/human/samples.tsv --from raw \
  --index references/human_rebuilt --exons references/human_rebuilt/exon.info \
  --out output/human_guarded --threads 2 --jobs 1
```

`--output` 与 `nju-seq --out` 必须指向同一新目录，父目录先建好。限制触发时终止进程组并写资源记录；输出目录/内存限制按间隔采样，单文件大小另有硬限制。正式全量任务应自行设置适用的上限。输出不支持覆盖。
