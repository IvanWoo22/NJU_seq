"""Exercise portable preparation against malformed inputs and relocated fixtures."""
import gzip
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from subset_fastq import subset


class Subsets(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='nju-subset-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = self.root / 'samples.tsv'
        self.manifest.write_text('sample\trole\tr1\tr2\talignment\tumi\n'
                                 'a\tNC\tR1.fq\tR2.fq\t-\tnone\n'
                                 'b\tTR\tR1.fq\tR2.fq\t-\tr1_3p8\n')
        for mate in (1, 2):
            (self.root / f'R{mate}.fq').write_text(''.join(
                f'@read{i}/{mate}\nACGT\n+\nIIII\n' for i in range(4)))

    def test_relocatable_paired_manifest_and_deterministic_gzip(self):
        first, second = self.root / 'first', self.root / 'second'
        report = subset(self.manifest, first, pairs=2, skip=1)
        subset(self.manifest, second, pairs=2, skip=1)
        self.assertEqual(report['samples'][0]['pairs'], 2)
        self.assertEqual(report['samples'][1]['umi'], 'r1_3p8')
        for sample in ('a', 'b'):
            for mate in (1, 2):
                file = f'{sample}/R{mate}.fq.gz'
                self.assertEqual((first / file).read_bytes(), (second / file).read_bytes())
                with gzip.open(first / file, 'rt') as handle:
                    self.assertEqual(handle.readline().strip(), f'@read1/{mate}')
        first.rename(self.root / 'moved')
        subset(self.root / 'moved/samples.tsv', self.root / 'third', pairs=1)
        with self.assertRaisesRegex(ValueError, 'Refusing existing'):
            subset(self.manifest, second, pairs=1)

    def test_failure_does_not_publish_subset(self):
        variants = ['@wrong/2\nACGT\n+\nIIII\n', '@read0/2\nACGT\n+\n',
                    '@read0/2\nACGT\n+\nIII\n', '']
        for i, text in enumerate(variants):
            (self.root / 'R2.fq').write_text(text)
            out = self.root / f'bad{i}'
            with self.assertRaises(ValueError):
                subset(self.manifest, out, pairs=1)
            self.assertFalse(out.exists())

    def test_invalid_design_paths_and_sample_names(self):
        original = self.manifest.read_text()
        for text in [original.replace('r1_3p8', 'guessed'), original.replace('a\tNC', '../a\tNC'),
                     original.replace('a\tNC', 'b\tNC'), original.replace('NC', 'TR')]:
            self.manifest.write_text(text)
            with self.assertRaises(ValueError):
                subset(self.manifest, self.root / 'bad', pairs=1)
        self.assertFalse((self.root / 'bad').exists())

    def test_source_gzip_error_and_short_input(self):
        (self.root / 'bad.gz').write_bytes(b'not gzip')
        self.manifest.write_text(self.manifest.read_text().replace('R1.fq', 'bad.gz'))
        with self.assertRaises(gzip.BadGzipFile):
            subset(self.manifest, self.root / 'bad', pairs=1)
        self.manifest.write_text(self.manifest.read_text().replace('bad.gz', 'R1.fq'))
        with self.assertRaisesRegex(ValueError, 'EOF'):
            subset(self.manifest, self.root / 'short', pairs=10)


class Plans(unittest.TestCase):
    def test_umi_parity_uses_selected_binary(self):
        binary = Path(os.environ.get('NJU_NATIVE', ROOT / 'target/release/nju-seq')).resolve()
        with tempfile.TemporaryDirectory(prefix='nju-parity-') as name:
            folder = Path(name)
            marker = folder / 'called'
            wrapper = folder / 'selected binary'
            wrapper.write_text(f'#!/bin/sh\ntouch {shlex.quote(str(marker))}\n'
                               f'exec {shlex.quote(str(binary))} "$@"\n')
            wrapper.chmod(0o755)
            subprocess.run([sys.executable, str(ROOT / 'tests/umi_parity.py'), '--reads', '100'],
                           env=dict(os.environ, NJU_NATIVE=str(wrapper)),
                           capture_output=True, check=True)
            self.assertTrue(marker.is_file())

    def test_custom_adapters_and_mixed_index_rejected(self):
        binary = Path(os.environ.get('NJU_NATIVE', ROOT / 'target/release/nju-seq'))
        fixture = ROOT / 'tests/data/synthetic'
        with tempfile.TemporaryDirectory(prefix='nju-plan-') as name:
            folder = Path(name)
            for prefix in ('bac', 'rrna', 'mrna'):
                for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
                    (folder / f'{prefix}.{part}.bt2').touch()
            args = [str(binary), 'plan', '--samples', str(fixture / 'samples.tsv'), '--from', 'raw',
                    '--index', str(folder), '--exons', str(fixture / 'exon.info'),
                    '--out', str(folder / 'run'), '--adapter-r1', 'ACGTCAGT', '--adapter-r2', 'GTCAGTCA']
            result = subprocess.run(args, capture_output=True, text=True, check=True)
            plan = json.loads(result.stdout)
            cut = next(s for s in plan['samples'][0]['steps'] if s['stage'] == '01_cutadapt')['argv']
            self.assertEqual(cut[cut.index('-a') + 1], 'ACGTCAGT')
            self.assertEqual(cut[cut.index('-A') + 1], 'GTCAGTCA')
            self.assertFalse((folder / 'run').exists())
            (folder / 'mrna.1.bt2').rename(folder / 'mrna.1.bt2l')
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('six components of one format', result.stderr)
            self.assertFalse((folder / 'run').exists())


class RrnaCounts(unittest.TestCase):
    def test_original_human_references_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix='nju-rrna-') as name:
            folder = Path(name)
            alignment = folder / 'selected.tsv'
            alignment.write_text('q1\t18s\t51\t10=\tACGTACGTAC\n'
                                 'q2\t28s\t101\t3=1I6=\tACGTACGTAC\n')
            out = folder / 'NC'
            argv = [sys.executable, str(ROOT / 'scripts/count_rrna.py'),
                    '--input', str(alignment), '--out', str(out)]
            subprocess.run(argv, check=True, capture_output=True)
            audit = json.loads((out / 'audit.json').read_text())
            self.assertEqual({k: v['alignment_records'] for k, v in audit['counts'].items()},
                             {'18s': 1, '28s': 1, '5-8s': 0})
            rows = [r.split('\t') for r in (out / '28s.tsv').read_text().splitlines()]
            self.assertEqual(rows[100][2], '1')
            self.assertEqual(rows[108][3], '1')  # Insertion does not consume reference bases.
            result = subprocess.run(argv, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertTrue((out / 'SUCCESS').is_file())


if __name__ == '__main__':
    unittest.main()
