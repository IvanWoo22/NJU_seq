#!/usr/bin/env python3
"""Bound mRNA alignment batches while preserving all -a alignments and read order.

Uses the unchanged native align-select in each batch. No SAM is materialized.
Outputs are published only after all batches and paired gzip inputs succeed.
"""
import argparse
import gzip
import json
from pathlib import Path
import os
import shutil
import subprocess
import time
from chunk_checkpoint import (COUNTERS, load as load_checkpoint, prerequisites,
                              prerequisite_digest, verify_prerequisites)
from provenance import save


def record(handle):
    first = handle.readline()
    if not first:
        return None
    row = [first, handle.readline(), handle.readline(), handle.readline()]
    if not all(row) or not row[0].startswith(b'@') or not row[2].startswith(b'+'):
        raise ValueError('Truncated or malformed FASTQ')
    if len(row[1].rstrip(b'\r\n')) != len(row[3].rstrip(b'\r\n')):
        raise ValueError('FASTQ sequence/quality length mismatch')
    return row


def identity(row):
    name = row[0].split()[0][1:]
    return name[:-2] if name.endswith((b'/1', b'/2')) else name


def run(args):
    args.native = args.native.resolve()
    targets = [args.output, args.log, args.stats]
    if len({p.absolute() for p in targets}) != 3 or any(p.exists() or p.is_symlink() for p in targets):
        raise ValueError('Require three distinct new output paths')
    work = args.output.with_name(args.output.name + '.chunks')
    work.mkdir()
    accepted = work / 'accepted.partial'
    logpath = work / 'bowtie2.partial.log'
    counts = dict.fromkeys(COUNTERS, 0)
    batches = pairs = 0
    previous_id = None
    started = time.monotonic()
    checkpoint = getattr(args, 'resume_checkpoint', None)
    reused_pairs = reused_batches = 0
    if checkpoint:
        meta = load_checkpoint(checkpoint, args.native, args.index, args.exons, args.r1, args.r2)
        batches = reused_batches = meta['completed_batches']
        pairs = reused_pairs = meta['completed_pairs']
        counts = meta['counts'].copy()
        snapshot = meta['prerequisites']
    else:
        snapshot = prerequisites(args.native, args.index, args.exons, args.r1, args.r2)
    digest = prerequisite_digest(snapshot)
    save(work / 'prerequisites.json', snapshot)
    def verify(*, hashes=False):
        verify_prerequisites(snapshot, args.native, args.index, args.exons, args.r1, args.r2, hashes=hashes)
    verify(hashes=True)
    def progress(status, **extra):
        data = dict(status=status, completed_batches=batches, completed_pairs=pairs,
                    accepted_reads=counts['selected_reads'], elapsed_seconds=time.monotonic()-started,
                    threads=args.threads, reused_batches=reused_batches, reused_pairs=reused_pairs,
                    prerequisites_sha256=digest, **extra)
        temporary = work / 'progress.partial.json'
        temporary.write_text(json.dumps(data, indent=2) + '\n')
        temporary.replace(work / 'progress.json')
    def source(path):
        return gzip.open(path, 'rb') if str(path).endswith('.gz') else path.open('rb')
    with source(args.r1) as one, source(args.r2) as two, accepted.open('xb') as dst, \
            logpath.open('xb') as log, (work / 'batches.jsonl').open('x') as audit:
        if checkpoint:
            with (checkpoint/'accepted.tsv').open('rb') as src:
                shutil.copyfileobj(src, dst)
            with (checkpoint/'bowtie2.log').open('rb') as src:
                shutil.copyfileobj(src, log)
            with (checkpoint/'batches.jsonl').open() as src:
                shutil.copyfileobj(src, audit)
            audit.flush(); dst.flush(); log.flush()
            progress('validating_checkpoint_input', checked_prefix_pairs=0)
            for i in range(pairs):
                x, y = record(one), record(two)
                if x is None or y is None or not identity(x) or identity(x) != identity(y):
                    raise ValueError(f'Checkpoint FASTQ mismatch at pair {i + 1}')
                if identity(x) == previous_id:
                    raise ValueError('Adjacent duplicate fragment ID in checkpoint prefix')
                previous_id = identity(x)
                if (i + 1) % 1000000 == 0:
                    progress('validating_checkpoint_input', checked_prefix_pairs=i + 1)
            progress('running')
        while True:
            verify()
            folder = work / f'{batches:06d}'
            folder.mkdir()
            r1, r2 = folder / 'R1.fq', folder / 'R2.fq'
            size = 0
            with r1.open('xb') as a, r2.open('xb') as b:
                while size < args.chunk_pairs:
                    x, y = record(one), record(two)
                    if x is None and y is None:
                        break
                    if x is None or y is None or not identity(x) or identity(x) != identity(y):
                        raise ValueError(f'FASTQ mate mismatch at pair {pairs + size + 1}')
                    if identity(x) == previous_id:
                        raise ValueError('Adjacent duplicate fragment ID; cannot preserve a unique read group')
                    previous_id = identity(x)
                    a.write(b''.join(x)); b.write(b''.join(y))
                    size += 1
            if not size:
                shutil.rmtree(folder)
                break
            argv = [str(args.native), 'align-select', '--mode', 'mrna', '--index', str(args.index),
                    '--r1', str(r1), '--r2', str(r2), '--threads', str(args.threads),
                    '--exons', str(args.exons), '--output', str(folder / 'accepted.tsv'),
                    '--log', str(folder / 'bowtie2.log'), '--stats', str(folder / 'selection.json')]
            tick = time.monotonic()
            with (folder / 'command.log').open('xb') as command_log:
                child = subprocess.run(argv, stdout=command_log, stderr=subprocess.STDOUT)
            if child.returncode:
                raise RuntimeError(f'Batch {batches} failed with exit {child.returncode}; see {folder}')
            verify()
            stats = json.loads((folder / 'selection.json').read_text())
            if stats['sam_materialized'] or stats['screen_kind'] != 'none':
                raise ValueError('Unexpected alignment policy')
            for key in COUNTERS:
                counts[key] += stats[key]
            with (folder / 'accepted.tsv').open('rb') as src:
                shutil.copyfileobj(src, dst)
            log.write(f'\n=== batch {batches}, pairs {pairs + 1}-{pairs + size} ===\n'.encode())
            with (folder / 'bowtie2.log').open('rb') as src:
                shutil.copyfileobj(src, log)
            audit.write(json.dumps(dict(batch=batches, first_pair=pairs + 1, pairs=size, argv=argv,
                                        elapsed_seconds=time.monotonic() - tick, stats=stats,
                                        prerequisites_sha256=digest)) + '\n')
            audit.flush(); dst.flush(); log.flush()
            pairs += size; batches += 1
            progress('running')
            shutil.rmtree(folder)
    verify(hashes=True)
    counts.update(sam_materialized=False, screen_multimap=None, screen_kind='none',
                  chunk_pairs=args.chunk_pairs, completed_batches=batches, input_pairs=pairs)
    if checkpoint:
        counts.update(reused_pairs=reused_pairs, reused_batches=reused_batches)
    statpath = work / 'selection.partial.json'
    statpath.write_text(json.dumps(counts) + '\n')
    published = []
    try:
        for source, target in zip([accepted, logpath, statpath], targets):
            os.link(source, target)
            published.append(target)
    except BaseException:
        for target in published:
            target.unlink()
        raise
    for p in [accepted, logpath, statpath]:
        p.unlink()
    progress('success')
    print(json.dumps(counts), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('native', 'index', 'r1', 'r2', 'exons', 'output', 'log', 'stats'):
        parser.add_argument('--' + option, required=True, type=Path)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--chunk-pairs', type=int, default=20000)
    parser.add_argument('--resume-checkpoint', type=Path)
    args = parser.parse_args()
    if args.threads < 1 or args.chunk_pairs < 1:
        parser.error('threads and chunk-pairs must be positive')
    run(args)
