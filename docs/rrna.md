# rRNA 独立计数与评分

主 `run` 的人 rRNA 步骤输出删除 ID 供 mRNA 分析；删除 ID 不是 rRNA 计数表。独立 rRNA 分支需保留筛选后的五列比对行，计数与评分均在与比对完全相同的参考坐标上执行。

以下从一个样本已去接头、完成可选 UMI 提取的 FASTQ 开始；若要与主流程的人 rRNA 步骤精确对照，还须先应用相同的细菌 rRNA 排除。根据实际输入替换路径，输出目录先创建。Linux/macOS shell 使用 `pipefail` 检查 Bowtie2 与 Rust 两端状态：

```bash
set -euo pipefail
mkdir -p output/rrna_work
bowtie2 --end-to-end --reorder --seed 0 -t --no-unal --no-mixed --no-discordant \
  -D 20 -R 3 -N 0 -L 10 -i S,1,0.75 --xeq -p 4 -k 5 --np 0 \
  -x references/human_rebuilt/rrna \
  -1 /path/to/NJU6096.after_bac.R1.fq.gz \
  -2 /path/to/NJU6096.after_bac.R2.fq.gz \
  2> output/rrna_work/NJU6096.bowtie2.log \
  | target/release/nju-seq filter --input /dev/stdin --format sam --mode rrna \
      --output output/rrna_work/NJU6096.selected.tsv

# 默认逐条读取原 data/hsa_rrna/{18s,28s,5-8s}.fa。
python3 scripts/count_rrna.py --input output/rrna_work/NJU6096.selected.tsv \
  --out output/rrna_counts/NJU6096
```

先确认整个管道退出码为零，再计数：上游失败时，筛选程序可能已看到合法的部分流，不能仅据 TSV 存在认定成功。其他参考用重复的 `--reference FILE` 显式指定；每个文件只能有一条参考，其 ID 和序列须匹配实际索引。计数单位是筛选后的定位记录，多个有效 rRNA 定位沿用原政策；不是 UMI 分子数。

按相同处理步骤分别生成 NC 与所有 TR 的计数，再联合评分：

```bash
target/release/nju-seq score --mode rrna --reference-id 18s \
  --nc output/rrna_counts/NJU6096/18s.tsv \
  --tr output/rrna_counts/NJU6097/18s.tsv \
  --tr output/rrna_counts/NJU6098/18s.tsv \
  --tr output/rrna_counts/NJU6099/18s.tsv \
  --output output/rrna_counts/18s.scored.tsv \
  --summary output/rrna_counts/18s.score_summary.json
```

28S、5.8S 分别使用 `28s`、`5-8s` 及对应计数表。样本名称取自计数文件父目录，必须唯一。默认输出所有候选坐标行及评价状态；只需通过阈值的行可加 `--rows candidates`。详见 [评分说明](scoring.md)。
