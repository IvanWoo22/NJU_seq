#!/usr/bin/env python3
"""Create a bounded, deterministic paired FASTQ subset from a six-column manifest.

Read only the requested prefix. This is a workflow fixture, not random sampling
or a validation of the unconsumed source gzip stream.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re


def read_record(handle):
    first = handle.readline()
    if not first:
        return None
    row = [first, handle.readline(), handle.readline(), handle.readline()]
    if not all(row) or not first.startswith(b'@') or not row[2].startswith(b'+'):
        raise ValueError('Malformed or truncated FASTQ record')
    seq, qual = (row[i].rstrip(b'\r\n') for i in (1, 3))
    if not seq or len(seq) != len(qual):
        raise ValueError('FASTQ sequence/quality length mismatch or empty sequence')
    return row


def fragment_id(row):
    name = row[0].split()[0][1:]
    if name.endswith((b'/1', b'/2')):
        name = name[:-2]
    if not name:
        raise ValueError('Empty FASTQ fragment ID')
    return name


def read_manifest(path):
    path = path.resolve()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['sample', 'role', 'r1', 'r2', 'alignment', 'umi']:
            raise ValueError('Expected sample, role, r1, r2, alignment, umi TSV columns')
        rows = list(reader)
    names = set()
    for row in rows:
        name = row['sample']
        if (None in row or None in row.values() or not re.fullmatch(r'[A-Za-z0-9_.-]+', name)
                or name in ('.', '..') or name in names):
            raise ValueError(f'Invalid or repeated sample: {name!r}')
        names.add(name)
        if row['role'] not in ('NC', 'TR') or row['umi'] not in ('none', 'r1_3p8'):
            raise ValueError(f'Invalid role or UMI design: {name}')
        if row['alignment'] != '-':
            raise ValueError('Subsets require a FASTQ manifest with alignment=-')
        for mate in ('r1', 'r2'):
            source = (path.parent / row[mate]).resolve()
            if row[mate] == '-' or not source.is_file():
                raise ValueError(f'Missing {name} {mate}: {source}')
            row[mate] = str(source)
        if row['r1'] == row['r2']:
            raise ValueError(f'R1 and R2 must be distinct inputs: {name}')
    if sum(r['role'] == 'NC' for r in rows) != 1 or not any(r['role'] == 'TR' for r in rows):
        raise ValueError('Require exactly one NC and at least one TR')
    return rows


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def source_info(path):
    p = Path(path)
    s = p.stat()
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def subset(manifest, out, pairs=1000, skip=0):
    if pairs < 1 or skip < 0:
        raise ValueError('pairs must be positive; skip must be nonnegative')
    rows = read_manifest(manifest)
    out = out.absolute()
    stage = out.with_name(out.name + '.partial')
    if any(p.exists() or p.is_symlink() for p in (out, stage)):
        raise ValueError(f'Refusing existing output or staging directory: {out}')
    out.parent.mkdir(parents=True, exist_ok=True)
    stage.mkdir()
    report = dict(method='consecutive prefix after skip; not random', pairs_per_sample=pairs,
                  skipped_pairs=skip, source_gzip_scope='consumed prefix only; full CRC not checked',
                  samples=[])
    with (stage / 'samples.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            folder = stage / row['sample']
            folder.mkdir()
            sources = [source_info(row[m]) for m in ('r1', 'r2')]
            seen = set()
            def opener(p):
                return gzip.open(p, 'rb') if p.endswith('.gz') else open(p, 'rb')
            with opener(row['r1']) as one, opener(row['r2']) as two, \
                    (folder / 'R1.fq.gz').open('xb') as raw1, (folder / 'R2.fq.gz').open('xb') as raw2, \
                    gzip.GzipFile(filename='', fileobj=raw1, mode='wb', mtime=0) as dst1, \
                    gzip.GzipFile(filename='', fileobj=raw2, mode='wb', mtime=0) as dst2:
                for i in range(skip + pairs):
                    a, b = read_record(one), read_record(two)
                    if a is None or b is None:
                        raise ValueError(f"{row['sample']}: EOF before requested pair {i + 1}")
                    name = fragment_id(a)
                    if name != fragment_id(b):
                        raise ValueError(f"{row['sample']}: mismatched mates at pair {i + 1}")
                    if i < skip:
                        continue
                    if name in seen:
                        raise ValueError(f"{row['sample']}: duplicate fragment ID in subset")
                    seen.add(name)
                    dst1.write(b''.join(a))
                    dst2.write(b''.join(b))
            if sources != [source_info(row[m]) for m in ('r1', 'r2')]:
                raise ValueError('Source changed during sampling')
            hashes = {f'R{m}.fq.gz': digest(folder / f'R{m}.fq.gz') for m in (1, 2)}
            report['samples'].append(dict(sample=row['sample'], role=row['role'], umi=row['umi'],
                                          pairs=len(seen), sources=sources, sha256=hashes))
            writer.writerow(dict(row, r1=f"{row['sample']}/R1.fq.gz", r2=f"{row['sample']}/R2.fq.gz"))
    (stage / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    stage.rename(out)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--pairs', type=int, default=1000)
    parser.add_argument('--skip', type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(subset(args.samples, args.out, args.pairs, args.skip), indent=2))
