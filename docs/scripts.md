# Helper scripts

**English** | [简体中文](scripts.zh-CN.md)

These scripts prepare inputs, align reads in batches, and limit resource usage. Run all commands from the repository root and replace example paths with your own inputs.

| Task | Script | Main inputs and outputs |
|---|---|---|
| Download and build references | `scripts/build_references.py` | Pinned public downloads → bacterial, human, or Arabidopsis references; see [construction commands](../references/README.md) |
| Prepare human references | `scripts/prepare_human_reference.py` | Transcripts, exon table, genome, and optional bacterial index → matching reference directory |
| Subset paired reads | `scripts/subset_fastq.py` | Six-column sample manifest → paired FASTQ and a new manifest |
| Align mRNA in batches | `scripts/align_chunked.py` | FASTQ after trimming and rRNA removal → selected five-column alignment table |
| Export an alignment checkpoint | `scripts/chunk_checkpoint.py` | Completed batches → checkpoint directory for resuming alignment |
| Count rRNA ends | `scripts/count_rrna.py` | Five-column alignment table and single-sequence rRNA FASTA → four-column count table |
| Limit resource usage | `scripts/guard_run.py` | Analysis command and disk, memory, and time limits → results and resource statistics |

Complete reference preparation and subsetting commands are in the [main README](../README.md); rRNA commands are in [rRNA counting and scoring](rrna.md). Each script supports `--help`.

## Running individual steps

`nju-seq run` performs the workflow from FASTQ to scores. You can also run individual commands on existing inputs:

| Task | Command |
|---|---|
| Align reads and select alignments | `nju-seq align-select` |
| Filter an existing SAM or five-column alignment table | `nju-seq filter` |
| Select each read's alignment using its R1 sequence | `nju-seq almostunique` |
| Remove both mates by read ID | `nju-seq delete-pairs` |
| Group UMIs by exact matching | `nju-seq umi` |
| Count ends, merge at genomic coordinates, and jointly score samples | `nju-seq count`, `merge`, `score` |
| Calculate FASTQ / alignment QC metrics and generate HTML reports | `nju-seq qc`, `qc-report` |

See `target/release/nju-seq --help` for options. `--index DIR` uses the prefixes `bac`, `rrna`, and `mrna`; specify `--bac-index`, `--rrna-index`, and `--mrna-index` separately for other names.

For `raw` and `prepared` inputs, `--rrna-filter host` skips bacterial alignment and removal and requires only the host rRNA and mRNA indices. `both` (default), `bac`, and `none` select both, bacterial only, or neither rRNA filter. See [rRNA filtering](../README.md#rrna-filtering). To prepare references without a bacterial index, use `build_references.py human|arabidopsis --without-bacterial`, or omit `--bac-index` from `prepare_human_reference.py`.

## Chunked mRNA alignment

Input FASTQ must already have been trimmed, had any UMIs extracted, and had rRNA reads removed. The following example processes one sample:

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

`--chunk-pairs` sets the number of complete read pairs per batch; `--threads` sets alignment threads. Each batch still uses mRNA `-a` alignment, so samples with many multiple matches may use substantial memory. `--output`, `--log`, and `--stats` must name three distinct files that do not already exist.

After processing each sample, prepare a six-column manifest for `--from unique`: set `r1` to the R1 used for alignment, `r2` to `-` if omitted, `alignment` to `accepted.tsv`, and `umi` to the library's UMI type. Then run `nju-seq run --from unique` for counting and scoring.

## Resuming interrupted batches

First stop any alignment process that is still writing files. If at least one batch has completed, export a checkpoint from `accepted.tsv.chunks/`:

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

Keep the FASTQ, index, exon table, and Rust/Bowtie2 executables unchanged, and resume with new output paths:

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

A checkpoint contains completed batches only; incomplete batches are processed again. This feature resumes mRNA alignment performed by `align_chunked.py`. `nju-seq run` does not resume an entire interrupted workflow.

## Linux resource limits

The following human-subset example limits output directory size, individual file size, free disk space, and elapsed time:

```bash
mkdir -p output
python3 scripts/guard_run.py --output output/human_guarded \
  --max-gib 4 --max-file-gib 1 --min-free-gib 20 --timeout 1800 -- \
  target/release/nju-seq run --samples tests/data/human/samples.tsv --from raw \
  --index references/human_rebuilt --exons references/human_rebuilt/exon.info \
  --out output/human_guarded --threads 2 --jobs 1
```

| Option | Meaning |
|---|---|
| `--max-gib` | Maximum size of the output directory, in GiB |
| `--max-file-gib` | Hard size limit for each output file, in GiB |
| `--min-free-gib` | Minimum free space on the output filesystem, in GiB |
| `--timeout` | Maximum elapsed time, in seconds |
| `--max-rss-gib` | Maximum resident memory of the process group, in GiB; 0 disables this limit |
| `--poll` | Resource-check interval, in seconds |

`--output` and `nju-seq --out` must point to the same new directory, whose parent must already exist. On reaching a limit, the script terminates the process group and writes `<directory-name>.resources.json` beside the output directory. Directory size and memory are sampled at intervals. Adjust the limits for the size of full datasets.
