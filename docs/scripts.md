# 辅助脚本

以下脚本用于准备数据、分批比对和控制资源占用。所有命令在仓库根目录运行，文件路径按实际样本替换。

| 用途 | 脚本 | 主要输入与输出 |
|---|---|---|
| 准备人源参考 | `scripts/prepare_human_reference.py` | 转录本、exon 表、基因组及细菌索引 → 配套参考目录 |
| 抽取双端小样本 | `scripts/subset_fastq.py` | 六列样本清单 → 配对 FASTQ 和新的样本清单 |
| 分批 mRNA 比对 | `scripts/align_chunked.py` | 去接头、去 rRNA 后的 FASTQ → 筛选后的五列定位表 |
| 导出比对断点 | `scripts/chunk_checkpoint.py` | 已完成批次 → 可供继续运行的断点目录 |
| rRNA 计数 | `scripts/count_rrna.py` | 五列定位表、单条 rRNA FASTA → 四列计数表 |
| 限制资源占用 | `scripts/guard_run.py` | 分析命令及磁盘、内存、时间上限 → 分析结果及资源统计 |

参考准备和抽样的完整命令见 [主 README](../README.md)，rRNA 用法见 [rRNA 计数与评分](rrna.md)。各脚本支持 `--help`。

## 分步处理

`nju-seq run` 可完成从 FASTQ 到评分的分析，也可根据已有输入调用单独的命令：

| 操作 | 命令 |
|---|---|
| 比对并筛选定位 | `nju-seq align-select` |
| 筛选已有 SAM 或五列比对表 | `nju-seq filter` |
| 结合 R1 序列筛选每个 read 的定位 | `nju-seq almostunique` |
| 按 ID 同时删除双端 reads | `nju-seq delete-pairs` |
| exact UMI 合并 | `nju-seq umi` |
| 末端计数、基因组坐标合并、联合评分 | `nju-seq count`、`merge`、`score` |
| FASTQ / 定位表质量统计、HTML 报告 | `nju-seq qc`、`qc-report` |

参数用法见 `target/release/nju-seq --help`。`--index DIR` 使用 `bac`、`rrna`、`mrna` 三个索引前缀；其他命名方式分别指定 `--bac-index`、`--rrna-index`、`--mrna-index`。

## 分批 mRNA 比对

输入为已经去接头、完成可选 UMI 提取并排除 rRNA 的双端 FASTQ。以下以一个样本为例：

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

`--chunk-pairs` 指定每批的完整 read pairs 数；`--threads` 指定比对线程数。每批仍使用 mRNA `-a` 全定位比对，多重匹配较多的样本可能需要较多内存。`--output`、`--log`、`--stats` 指向三个不同的新文件。

各样本完成后，准备 `--from unique` 的六列样本清单：`r1` 指向该样本比对时使用的 R1，`r2` 可填 `-`，`alignment` 指向 `accepted.tsv`，`umi` 填文库实际类型。随后用 `nju-seq run --from unique` 完成计数与评分。

## 从中断的批次继续

先停止正在写入的比对进程。若至少有一个批次完成，可从 `accepted.tsv.chunks/` 导出断点：

```bash
python3 scripts/chunk_checkpoint.py \
  --chunks output/batch_sample/accepted.tsv.chunks \
  --accepted output/batch_sample/accepted.tsv.chunks/accepted.partial \
  --log output/batch_sample/accepted.tsv.chunks/bowtie2.partial.log \
  --native "$PWD/target/release/nju-seq" \
  --index references/human_rebuilt/mrna --exons references/human_rebuilt/exon.info \
  --r1 /path/to/sample.after_rrna.R1.fq.gz \
  --r2 /path/to/sample.after_rrna.R2.fq.gz \
  --output output/sample_checkpoint
```

保持 FASTQ、索引、exon 表及 Rust/Bowtie2 程序一致，使用新输出目录继续：

```bash
mkdir -p output/batch_sample_resumed
python3 scripts/align_chunked.py \
  --native "$PWD/target/release/nju-seq" \
  --index references/human_rebuilt/mrna --exons references/human_rebuilt/exon.info \
  --r1 /path/to/sample.after_rrna.R1.fq.gz \
  --r2 /path/to/sample.after_rrna.R2.fq.gz \
  --output output/batch_sample_resumed/accepted.tsv \
  --log output/batch_sample_resumed/bowtie2.log \
  --stats output/batch_sample_resumed/selection.json --threads 4 --chunk-pairs 20000 \
  --resume-checkpoint output/sample_checkpoint
```

断点只包含已完成批次；未完成的批次会重新计算。此功能适用于 `align_chunked.py` 的 mRNA 比对，`nju-seq run` 不支持整条流程断点恢复。

## Linux 资源限制

以下以人源小样本为例，限制输出目录大小、单文件大小、可用磁盘空间和运行时间：

```bash
mkdir -p output
python3 scripts/guard_run.py --output output/human_guarded \
  --max-gib 4 --max-file-gib 1 --min-free-gib 20 --timeout 1800 -- \
  target/release/nju-seq run --samples tests/data/human/samples.tsv --from raw \
  --index references/human_rebuilt --exons references/human_rebuilt/exon.info \
  --out output/human_guarded --threads 2 --jobs 1
```

| 参数 | 含义 |
|---|---|
| `--max-gib` | 输出目录的大小上限，GiB |
| `--max-file-gib` | 每个输出文件的大小硬限制，GiB |
| `--min-free-gib` | 输出所在磁盘的最低可用空间，GiB |
| `--timeout` | 运行时间上限，秒 |
| `--max-rss-gib` | 进程组常驻内存上限，GiB；0 表示不限制 |
| `--poll` | 资源检查间隔，秒 |

`--output` 与 `nju-seq --out` 须指向同一新目录，父目录需已存在。达到限制时，脚本终止进程组，并在输出目录旁写入 `<目录名>.resources.json`。输出目录和内存按间隔采样；运行全量数据时，按样本规模调整这些上限。
