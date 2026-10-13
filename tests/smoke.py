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


def check_counts(out, truth):
    for sample, expected in truth['samples'].items():
        folder = out / sample
        rows = [r.split('\t') for r in (folder / 'reads/mrna.filter.tsv').read_text().splitlines()]
        assert {r[1]: int(r[8]) for r in rows} == expected['genome_end_counts'], (out, sample)
        assert len((folder / 'mrna.almostuniquematch.filter.tmp').read_text().splitlines()) == expected['accepted']


def check_optional_filters(work, reference, truth, log):
    prepared = work / 'prepared.tsv'
    prepared.write_text('sample\trole\tr1\tr2\talignment\tumi\n' + ''.join(
        f'{s}\t{"NC" if s == "NC" else "TR"}\t{work}/all/{s}/R1.cut.fq.gz\t'
        f'{work}/all/{s}/R2.cut.fq.gz\t-\tnone\n' for s in truth['samples']))
    for mode, enabled in [('both', ('bac', 'rrna')), ('host', ('rrna',)),
                          ('bac', ('bac',)), ('none', ())]:
        # Only the indices selected by the user exist in this package.
        indices = work / f'indices-{mode}'
        indices.mkdir()
        for prefix in (*enabled, 'mrna'):
            for source in reference.glob(f'{prefix}.*.bt2'):
                (indices / source.name).symlink_to(source)
        for entry, storage, qc in [('raw', 'compact', 'yes'), ('raw', 'all', 'no'),
                                   ('prepared', 'compact', 'no'), ('prepared', 'all', 'yes')]:
            out = work / f'{mode}-{entry}-{storage}'
            command([BIN, 'run', '--samples', FIXTURE / 'samples.tsv' if entry == 'raw' else prepared,
                     '--from', entry, '--index', indices, '--rrna-filter', mode,
                     '--exons', reference / 'exon.info', '--out', out, '--threads', '1', '--jobs', '2',
                     '--storage', storage, '--qc', qc], log)
            assert (out / 'SUCCESS').is_file(), out
            check_counts(out, truth)
            plan = json.loads((out / 'plan.json').read_text())
            assert plan['rrna_filter'] == mode
            assert (out / 'qc.html').exists() == (qc == 'yes')
            for sample, expected in truth['samples'].items():
                folder = out / sample
                remaining = expected['raw_pairs']
                for prefix, field, stage in [('bac', 'bac_removed', 'after_bac'),
                                             ('rrna', 'host_removed', 'after_host')]:
                    if prefix in enabled:
                        remaining -= expected[field]
                        selected = folder / f'{prefix}.{"selected.tsv" if entry == "raw" and qc == "no" else "removed.ids"}'
                        assert len(selected.read_text().splitlines()) == expected[field], (out, sample)
                        if qc == 'yes':
                            with (folder / f'qc_{stage}/summary.tsv').open() as src:
                                rows = list(csv.DictReader(src, delimiter='\t'))
                            assert int(next(r['reads'] for r in rows if r['mate'] == 'R1')) == remaining
                    else:
                        assert not list(folder.glob(f'{prefix}*')), (out, sample, prefix)
                        assert not (folder / f'qc_{stage}').exists()
                if qc == 'yes':
                    # Skipped rRNA reads must reach mRNA alignment even though they do not map there.
                    assert f'{remaining} reads; of these:' in (folder / 'mrna.bowtie2.log').read_text()
                if storage == 'all':
                    for prefix, suffix in [('bac', 'filter'), ('rrna', 'filter.mrna')]:
                        assert (folder / f'R1.{suffix}.fq.gz').exists() == (prefix in enabled)

    if not shutil.which('umi_tools'):
        print('SKIP: UMI extraction checks (install UMI-tools to include them)')
        return
    # Exercise actual UMI extraction and molecule counting with no bacterial index.
    tagged = work / 'tagged'
    tagged.mkdir()
    manifest = ['sample\trole\tr1\tr2\talignment\tumi']
    for sample in truth['samples']:
        for mate in (1, 2):
            lines = (FIXTURE / f'{sample}.R{mate}.fq').read_text().splitlines()
            with (tagged / f'{sample}.R{mate}.fq').open('w') as dst:
                for i in range(0, len(lines), 4):
                    sequence = (lines[i + 1].removesuffix('AGATCGGAAGAGCACA') + 'ACGTCTGAAGATCGGAAGAGCACA'
                                if mate == 1 else 'TCAGACGT' + lines[i + 1])
                    dst.write(f'{lines[i].removesuffix(f"/{mate}")}\n{sequence}\n+\n{"I" * len(sequence)}\n')
        manifest.append(f'{sample}\t{"NC" if sample == "NC" else "TR"}\t'
                        f'{sample}.R1.fq\t{sample}.R2.fq\t-\tr1_3p8')
    (tagged / 'samples.tsv').write_text('\n'.join(manifest) + '\n')
    for storage in ('compact', 'all'):
        out = work / f'umi-{storage}'
        command([BIN, 'run', '--samples', tagged / 'samples.tsv', '--from', 'raw',
                 '--rrna-filter', 'host', '--rrna-index', reference / 'rrna',
                 '--mrna-index', reference / 'mrna', '--exons', reference / 'exon.info',
                 '--out', out, '--threads', '1', '--jobs', '2', '--storage', storage], log)
        check_counts(out, truth)
        for sample in truth['samples']:
            folder = out / sample
            assert not list(folder.glob('bac*'))
            assert len((folder / 'mrna.umi.filter.tmp').read_text().splitlines()) == 3
        assert (out / 'mrna.umi.scored.tsv').is_file()
        assert (out / 'SUCCESS').is_file()
    for branch in ('reads', 'umi'):
        assert (work / f'umi-compact/mrna.{branch}.scored.tsv').read_bytes() == (work / f'umi-all/mrna.{branch}.scored.tsv').read_bytes()
    print('PASS: UMI extraction and counting without a bacterial index')


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
            check_optional_filters(work, standalone, truth, log)
    print('PASS: reference relocation; raw/prepared pipelines; all four rRNA filter modes; '
          'read counts; scores; QC; compact/all parity')


if __name__ == '__main__':
    main()
