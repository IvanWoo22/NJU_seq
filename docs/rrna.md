# rRNA 独立计数与评分

`nju-seq run` 排除 rRNA 后进行 mRNA 分析。需要 rRNA 位点结果时，单独保留筛选后的五列比对表，并在比对所用的参考序列上计数和评分。

以下使用已去接头、完成可选 UMI 提取并排除细菌 rRNA 的双端 FASTQ。将示例路径替换为实际输入；`pipefail` 使比对或筛选失败时停止后续计数：

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

# 默认读取 data/hsa_rrna/{18s,28s,5-8s}.fa。
python3 scripts/count_rrna.py --input output/rrna_work/NJU6096.selected.tsv \
  --out output/rrna_counts/NJU6096
```

使用其他参考时，重复传入 `--reference FILE`。每个 FASTA 文件只能包含一条序列，其 ID 和序列须与比对索引一致。每条通过筛选的定位单独计数，因此同一 read 的多个有效 rRNA 定位可分别贡献计数；这些计数以定位记录为单位，不进行 UMI 分子合并。

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
