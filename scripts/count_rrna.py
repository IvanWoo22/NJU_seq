#!/usr/bin/env python3
"""Count selected rRNA alignments with the maintained Perl algorithm and local FASTA.

Input is the five-column output of nju-seq filter --mode rrna. This step does not
perform adapter trimming, rRNA removal, or UMI molecule deduplication.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

from prepare_human_reference import fasta, fingerprint

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    references = []
    names = set()
    for path in args.reference:
        entries = list(fasta(path))
        if len(entries) != 1:
            raise ValueError(f'Expected single-record FASTA: {path}')
        name = entries[0][0].split()[0]
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', name) or name in ('.', '..') or name in names:
            raise ValueError(f'Unsafe or duplicate reference ID: {name}')
        names.add(name)
        references.append((name, path, entries[0][1]))
    # The preserved Perl reader accepts plain FASTA/TSV; normalize references below.
    if not args.input.is_file() or args.input.suffix == '.gz':
        raise ValueError('Input must be a plain five-column alignment TSV')
    source = fingerprint(args.input, digest=False)
    out = args.out.absolute()
    stage = out.with_name(out.name + '.partial')
    if any(p.exists() or p.is_symlink() for p in (out, stage)):
        raise ValueError(f'Refusing existing output: {out}')
    out.parent.mkdir(parents=True, exist_ok=True)
    stage.mkdir()
    (stage / 'references').mkdir()
    counts = {}
    for name, path, sequence in references:
        normalized = stage / 'references' / f'{name}.fa'
        normalized.write_text(f'>{name}\n{sequence}\n')
        argv = [args.perl, str(ROOT / 'rrna_analysis/readend_count.pl'),
                str(normalized), str(args.input.resolve()), name]
        with (stage / f'{name}.tsv').open('x') as output, (stage / f'{name}.log').open('x') as log:
            subprocess.run(argv, stdout=output, stderr=log, check=True)
        rows = [line.split('\t') for line in (stage / f'{name}.tsv').read_text().splitlines()]
        if len(rows) != len(sequence) or any(int(r[0]) != i + 1 or r[1] != sequence[i] for i, r in enumerate(rows)):
            raise ValueError(f'Invalid reference coordinates: {name}')
        starts, ends = (sum(int(r[i]) for r in rows) for i in (2, 3))
        if starts != ends:
            raise ValueError(f'Count conservation failed: {name}')
        counts[name] = dict(alignment_records=ends, length=len(sequence), reference=fingerprint(path))
    if source != fingerprint(args.input, digest=False):
        raise ValueError('Alignment input changed during counting')
    (stage / 'audit.json').write_text(json.dumps(dict(input=source, counts=counts,
                                                      unit='selected alignment records, not UMI molecules'), indent=2) + '\n')
    (stage / 'SUCCESS').write_text('rRNA counts completed\n')
    stage.rename(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--reference', type=Path, action='append', help='Repeat for each single-reference FASTA')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--perl', default=os.environ.get('NJU_PERL', 'perl'))
    args = parser.parse_args()
    args.reference = args.reference or [ROOT / 'data/hsa_rrna' / f'{n}.fa' for n in ('18s', '28s', '5-8s')]
    run(args)
