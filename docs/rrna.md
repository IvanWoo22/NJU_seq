# rRNA counting and scoring

**English** | [简体中文](rrna.zh-CN.md)

`nju-seq run` removes rRNA before mRNA analysis. To analyze rRNA sites separately, retain the selected five-column alignment table and count and score against the reference sequences used for alignment.

The following commands use paired FASTQ after adapter trimming, any UMI extraction, and bacterial rRNA removal. Replace the example paths with your inputs. `pipefail` stops subsequent counting if alignment or filtering fails.

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

# Uses data/hsa_rrna/{18s,28s,5-8s}.fa by default.
python3 scripts/count_rrna.py --input output/rrna_work/NJU6096.selected.tsv \
  --out output/rrna_counts/NJU6096
```

For other references, repeat `--reference FILE` for each FASTA. Each file must contain one sequence, with an ID and sequence matching the alignment index. Every accepted alignment contributes separately to the counts, so multiple valid rRNA alignments of one read can each contribute. These are alignment-record counts, without UMI molecule grouping.

Generate counts for the NC and all TR samples using the same steps, then score them jointly:

```bash
target/release/nju-seq score --mode rrna --reference-id 18s \
  --nc output/rrna_counts/NJU6096/18s.tsv \
  --tr output/rrna_counts/NJU6097/18s.tsv \
  --tr output/rrna_counts/NJU6098/18s.tsv \
  --tr output/rrna_counts/NJU6099/18s.tsv \
  --output output/rrna_counts/18s.scored.tsv \
  --summary output/rrna_counts/18s.score_summary.json
```

For 28S and 5.8S, use `28s` and `5-8s` with the corresponding count tables. Sample names come from the parent directories of the count files and must be unique. By default, output includes all candidate-coordinate rows and their evaluation status. Add `--rows candidates` to retain only rows passing the thresholds. See [Scoring](scoring.md).
