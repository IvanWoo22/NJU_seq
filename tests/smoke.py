#!/usr/bin/env python3
"""Check reference relocation and raw FASTQ analysis against known counts."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests/data/synthetic'
BIN = Path(os.environ.get('NJU_NATIVE', ROOT / 'target/release/nju-seq')).resolve()


def command(argv, log):
    result = subprocess.run(list(map(str, argv)), stdout=log, stderr=log)
    if result.returncode:
        log.seek(0)
        raise RuntimeError(f'Command failed: {argv}\n{log.read()}')


def main():
    truth = json.loads((FIXTURE / 'truth.json').read_text())
    with tempfile.TemporaryDirectory(prefix='nju-smoke-') as name:
        work = Path(name)
        with (work / 'commands.log').open('w+') as log:
            bundle = work / 'original bundle'
            source = bundle / 'sources'
            source.mkdir(parents=True)
            for name in ('mrna.fa', 'exon.info', 'genome.fa', 'bac.fa', 'rrna.fa'):
                shutil.copyfile(FIXTURE / name, source / name)
            command(['bowtie2-build', '--threads', '1', source / 'bac.fa', source / 'bac'], log)
            prepare = [sys.executable, ROOT / 'scripts/prepare_human_reference.py',
                       '--mrna', source / 'mrna.fa', '--exons', source / 'exon.info',
                       '--genome', source / 'genome.fa', '--bac-index', source / 'bac',
                       '--rrna', source / 'rrna.fa', '--threads', '1',
                       '--label', 'synthetic regression fixture']
            command([*prepare, '--out', bundle / 'reference'], log)
            command([*prepare, '--out', work / 'standalone', '--copy-inputs'], log)
            moved = work / 'moved bundle'
            bundle.rename(moved)
            standalone = work / 'copied reference'
            (work / 'standalone').rename(standalone)
            for link in (moved / 'reference').iterdir():
                if link.is_symlink():
                    assert not link.readlink().is_absolute(), link
                    assert link.is_file() and link.resolve().is_relative_to(moved), link
            assert not any(p.is_symlink() for p in standalone.iterdir())
            for storage in ('compact', 'all'):
                reference = moved / 'reference' if storage == 'compact' else standalone
                if storage == 'all':
                    shutil.rmtree(moved)  # The copied reference must work without any sources.
                out = work / storage
                command([BIN, 'run', '--samples', FIXTURE / 'samples.tsv', '--from', 'raw',
                         '--index', reference, '--exons', reference / 'exon.info', '--out', out,
                         '--threads', '1', '--jobs', '2', '--storage', storage,
                         '--qc', 'yes' if storage == 'compact' else 'no'], log)
                assert (out / 'SUCCESS').is_file()
                for sample, expected in truth['samples'].items():
                    folder = out / sample
                    rows = [r.split('\t') for r in (folder / 'reads/mrna.filter.tsv').read_text().splitlines()]
                    observed = {r[1]: int(r[8]) for r in rows}
                    assert observed == expected['genome_end_counts'], (sample, observed)
                    accepted = folder / 'mrna.almostuniquematch.filter.tmp'
                    assert len(accepted.read_text().splitlines()) == expected['accepted']
                    if storage == 'compact':
                        for mode, field in [('bac', 'bac_removed'), ('rrna', 'host_removed')]:
                            assert len((folder / f'{mode}.removed.ids').read_text().splitlines()) == expected[field]
                with (out / 'mrna.reads.scored.tsv').open() as handle:
                    candidates = list(csv.DictReader(handle, delimiter='\t'))
                assert [int(r['Position']) for r in candidates] == truth['candidate_positions'], candidates
                assert all(r['Candidate'] == 'true' for r in candidates)
            assert (work / 'compact/qc.html').is_file()
            assert not list((work / 'compact').glob('*/R*.fq.gz'))
            assert (work / 'compact/mrna.reads.scored.tsv').read_bytes() == (work / 'all/mrna.reads.scored.tsv').read_bytes()
    print('PASS: linked/copied reference relocation, raw pipeline, rRNA removal, counts, scores, QC, compact/all parity')


if __name__ == '__main__':
    main()
