# Human example data

**English** | [简体中文](README.zh-CN.md)

These files contain the first **1,000 read pairs** from each of NJU6096 (NC) and NJU6097, NJU6098, and NJU6099 (TR), all with `umi=none`. See [sample_list.csv](../../../sample_list.csv) for experimental details: HEK293T rRNA-free RNA, with no MgR treatment for NC and MgR treatment for TR.

The [sample manifest](samples.tsv) uses relative paths and can be copied with the data directory. [Provenance and checksums](audit.json) record the selected range and FASTQ SHA-256 values. This consecutive subset illustrates the workflow; it does not represent the full-library distribution or support estimates of full-library candidate counts.

Run the `--from raw` commands in the [main README](../../../README.md) with the matching human references. See [Example results](../../../docs/results.md) for stage counts and the meaning of `normalization_unavailable`.

If you have the complete FASTQ files, extract the same range from the repository root:

```bash
python3 scripts/subset_fastq.py --samples config/samples.human.tsv \
  --pairs 1000 --skip 0 --out output/human_subset
```
