#!/usr/bin/env python3
"""Small CLI acceptance tests for hg38 identities and bounded-run failures."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from subset_fastq import read_record, fragment_id
from prepare_human_reference import fasta

BIN = Path(os.environ.get('NJU_NATIVE', ROOT/'target/release/nju-seq'))
GUARD = ROOT/'scripts/guard_run.py'


class Cases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='nju-human-cases-')
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.root = Path(cls.tmp.name)
        ref = cls.root/'references'; ref.mkdir()
        for name in ('bac', 'rrna', 'mrna'):
            for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
                (ref/f'{name}.{part}.bt2').touch()

    def setUp(self):
        self.out = self.root/self._testMethodName
        self.out.mkdir()
        self.exons = self.write('exon.info',
            'chrX\t101\t110\t+\tParent=ENST00001.1;gene_id=ENSG00009.2\n'
            'chrX\t201\t210\t+\tParent=ENST00001.1;gene_id=ENSG00009.2\n'
            'chrY\t101\t110\t+\tParent=ENST00001.1_PAR_Y;gene_id=ENSG00009.2\n'
            'chrY\t201\t210\t+\tParent=ENST00001.1_PAR_Y;gene_id=ENSG00009.2\n'
            'chr2\t301\t310\t-\tParent=ENST00002.1;gene_id=ENSG00008.3\n'
            'chr2\t401\t410\t-\tParent=ENST00002.1;gene_id=ENSG00008.3\n')

    def write(self,name,text):
        p=self.out/name
        p.write_text(text)
        return p

    def command(self,*args,success=True):
        r=subprocess.run(list(map(str,args)),capture_output=True,text=True,timeout=30)
        with (self.out/'commands.jsonl').open('a') as f:
            f.write(json.dumps(dict(argv=list(map(str,args)),exit=r.returncode,stderr=r.stderr))+'\n')
        self.assertEqual(r.returncode==0,success,r.stderr)
        return r

    def test_reference_fasta_validation(self):
        source = self.write('valid.fa', '>28s description\nacgt\nNNry\n\n>5-8s\nACG\n')
        self.assertEqual(list(fasta(source)), [('28s description', 'ACGTNNRY'), ('5-8s', 'ACG')])
        malformed = [
            ('>28s\nACGT>5-8s\nACG\n', 'header must start a new line'),
            ('>5-8s first\nACG\n>5-8s second\nACG\n', 'duplicate FASTA ID'),
            ('>\nACG\n', 'empty FASTA header'),
            ('>28s\n>5-8s\nACG\n', 'empty FASTA record'),
            ('>28s\n', 'empty FASTA record'),
            ('ACGT\n>28s\nACG\n', 'sequence before FASTA header'),
            ('>28s\nAC1G\n', 'invalid DNA FASTA sequence character'),
            ('\n', 'empty FASTA file'),
        ]
        for n, (text, error) in enumerate(malformed):
            with self.subTest(error=error, case=n):
                with self.assertRaisesRegex(ValueError, error):
                    list(fasta(self.write(f'invalid{n}.fa', text)))

    def test_minimum_insert_length_with_and_without_umi(self):
        # Physical boundary reads, including a too-short R2, must exercise
        # actual Cutadapt and UMI extraction, not only compare planned argv.
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        fragment = 'TCCTAACGTTAGGCCATG'
        tag = 'ACGTCTGA'
        sizes = [(9, 10), (10, 9), (10, 10), (11, 11), (17, 17), (18, 18)]
        expected = {f'n{a}x{b}': (fragment[:a], rc(fragment[:b]))
                    for a, b in sizes if min(a, b) >= 10}
        libraries = [('NC', 'NC', 'none'), ('TR1', 'TR', 'r1_3p8'), ('TR2', 'TR', 'none')]
        manifest = ['sample\trole\tr1\tr2\talignment\tumi']
        for name, role, design in libraries:
            paths = [self.out/f'{name}.R{mate}.fq.gz' for mate in (1, 2)]
            with gzip.open(paths[0], 'wt') as one, gzip.open(paths[1], 'wt') as two:
                for a, b in sizes:
                    seq1 = fragment[:a] + (tag if design == 'r1_3p8' else '') + 'AGATCGGAAGAGCACA'
                    seq2 = (rc(tag) if design == 'r1_3p8' else '') + rc(fragment[:b]) + 'GATCGTCGGACTGTAG'
                    for dst, seq in [(one, seq1), (two, seq2)]:
                        dst.write(f'@n{a}x{b}\n{seq}\n+\n'+('I'*len(seq))+'\n')
            manifest.append('\t'.join([name, role, *map(str, paths), '-', design]))
        m = self.write('mixed.tsv', '\n'.join(manifest)+'\n')
        ref = self.root/'references'
        for storage in ('all', 'compact'):
            planned = self.out/storage
            result = self.command(BIN, 'plan', '--samples', m, '--from', 'raw', '--exons', self.exons,
                                  '--bac-index', ref/'bac', '--rrna-index', ref/'rrna', '--mrna-index', ref/'mrna',
                                  '--storage', storage, '--out', planned, '--threads', '1')
            for sample in json.loads(result.stdout)['samples']:
                (planned/sample['sample']).mkdir(parents=True)
                for step in sample['steps']:
                    if step['stage'] in ('01_cutadapt', '02_extract'):
                        self.command(step['program'], *step['argv'])
                suffix = 'origin' if sample['umi'] == 'r1_3p8' else 'cut'
                observed = []
                for mate in (1, 2):
                    records = {}
                    with gzip.open(planned/sample['sample']/f'R{mate}.{suffix}.fq.gz', 'rb') as f:
                        while (row := read_record(f)) is not None:
                            name = fragment_id(row).decode()
                            if sample['umi'] == 'r1_3p8':
                                name, observed_tag = name.rsplit('_', 1)
                                self.assertEqual(observed_tag, tag)
                            records[name] = row[1].strip().decode()
                    observed.append(records)
                self.assertEqual(observed[0], {name: pair[0] for name, pair in expected.items()})
                self.assertEqual(observed[1], {name: pair[1] for name, pair in expected.items()})

    def test_human_merge_gene_and_both_strands(self):
        source=self.write('count.tsv',
            'ENST00001.1|ENSG00009.2\t11\tA\tC\tG\t0\t2\n'
            'ENST00001.1_PAR_Y|ENSG00009.2_PAR_Y\t11\tA\tC\tG\t0\t3\n'
            'ENST00002.1|ENSG00008.3\t11\tT\tC\tA\t0\t4\n')
        target=self.out/'merged.tsv'
        self.command(BIN,'merge','--input',source,'--exons',self.exons,'--output',target)
        self.assertEqual(set(target.read_text().splitlines()),{
            'X\t201\t+\tA\tC\tG\tENSG00009\t0\t2',
            'Y\t201\t+\tA\tC\tG\tENSG00009\t0\t3',
            '2\t310\t-\tT\tC\tA\tENSG00008\t0\t4'})

    def test_par_multi_location_is_rejected_as_ambiguous(self):
        fq=self.write('r1.fq','@q\nAAA\n+\nIII\n')
        source=self.write('align.tsv',
            'q\tENST00001.1|ENSG00009.2\t1\t3=\tAAA\n'
            'q\tENST00001.1_PAR_Y|ENSG00009.2_PAR_Y\t1\t3=\tAAA\n')
        filtered=self.out/'filtered.tsv'
        self.command(BIN,'filter','--input',source,'--format','tsv','--mode','mrna',
                     '--exons',self.exons,'--output',filtered)
        self.assertEqual(len(filtered.read_text().splitlines()),2)
        accepted=self.out/'accepted.tsv'
        self.command(BIN,'almostunique','--input',filtered,'--r1',fq,'--output',accepted)
        self.assertEqual(accepted.read_text(),'')

    def test_invalid_annotation_missing_and_out_of_range_fail_atomically(self):
        for n,ref,start in [(0,'ENST00999.1',1),(1,'ENST00001.1',21),(2,'ENST00001.1_bad',1)]:
            source=self.write(f'bad{n}.tsv',f'q\t{ref}\t{start}\t3=\tAAA\n')
            target=self.out/f'bad{n}.out'
            self.command(BIN,'filter','--input',source,'--format','tsv','--mode','mrna',
                         '--exons',self.exons,'--output',target,success=False)
            self.assertFalse(target.exists())
        self.assertFalse(list(self.out.glob('*.partial')))

    def test_empty_TR_joint_score_succeeds_without_false_candidates(self):
        manifest=['sample\trole\tr1\tr2\talignment\tumi']
        for name,role,align in [('NC','NC','q\tENST00001.1|ENSG00009.2\t1\t3=\tAAA\n'),
                              ('TR1','TR',''),('TR2','TR',''),('TR3','TR','')]:
            r1=self.write(name+'.fq','@q\nAAA\n+\nIII\n')
            a=self.write(name+'.tsv',align)
            manifest.append(f'{name}\t{role}\t{r1}\t-\t{a}\tnone')
        m=self.write('samples.tsv','\n'.join(manifest)+'\n')
        out=self.out/'run'
        self.command(BIN,'run','--samples',m,'--from','unique','--exons',self.exons,'--out',out)
        self.assertTrue((out/'SUCCESS').is_file())
        self.assertEqual(len((out/'mrna.reads.scored.tsv').read_text().splitlines()),1)
        self.assertFalse((out/'mrna.umi.scored.tsv').exists())
        self.command(BIN,'run','--samples',m,'--from','unique','--exons',self.exons,'--out',out,success=False)

    def test_invalid_manifest_and_missing_indexes_fail_before_output(self):
        fq=self.write('r.fq','@q\nAAA\n+\nIII\n')
        for n,roles in enumerate([('TR','TR'),('NC','NC')]):
            m=self.write(f'm{n}.tsv','sample\trole\tr1\tr2\talignment\tumi\n'+
                ''.join(f'{s}\t{r}\t{fq}\t-\t{fq}\tnone\n' for s,r in zip(['a','b'],roles)))
            out=self.out/f'run{n}'
            self.command(BIN,'run','--samples',m,'--from','unique','--exons',self.exons,'--out',out,success=False)
            self.assertFalse(out.exists())
        self.command(BIN,'run','--samples',m,'--from','raw','--index',self.out/'missing',
                     '--exons',self.exons,'--out',self.out/'absent',success=False)
        self.assertFalse((self.out/'absent').exists())

    def test_filter_options_require_only_enabled_indices(self):
        fq = self.write('r.fq', '@q\nACGTACGTACGT\n+\nIIIIIIIIIIII\n')
        manifest = self.write('samples.tsv', 'sample\trole\tr1\tr2\talignment\tumi\n' +
            ''.join(f'{name}\t{role}\t{fq}\t{fq}\t-\tnone\n'
                    for name, role in [('NC', 'NC'), ('TR', 'TR')]))
        for mode, required in [('both', ('bac', 'rrna', 'mrna')), ('host', ('rrna', 'mrna')),
                               ('bac', ('bac', 'mrna')), ('none', ('mrna',))]:
            indices = [arg for prefix in required
                       for arg in (f'--{prefix}-index', self.root/'references'/prefix)]
            for entry in ('raw', 'prepared'):
                for storage in ('compact', 'all'):
                    for qc in ('yes', 'no'):
                        with self.subTest(mode=mode, entry=entry, storage=storage, qc=qc):
                            args = [BIN, 'plan', '--samples', manifest, '--from', entry,
                                    '--rrna-filter', mode, '--storage', storage, '--qc', qc,
                                    '--exons', self.exons, '--out', self.out/'absent']
                            self.command(*args, *indices)
                            self.assertFalse((self.out/'absent').exists())
                for missing in required:
                    args = [BIN, 'run', '--samples', manifest, '--from', entry, '--rrna-filter', mode,
                            '--exons', self.exons, '--out', self.out/'absent']
                    for prefix in required:
                        args += [f'--{prefix}-index', self.out/'missing' if prefix == missing
                                 else self.root/'references'/prefix]
                    self.command(*args, success=False)
                    self.assertFalse((self.out/'absent').exists())

    def test_pair_desynchronization_truncation_and_gzip_error(self):
        ids=self.write('ids','')
        good=self.write('good.fq','@q/1\nAAA\n+\nIII\n')
        bad=[self.write('other.fq','@z/2\nAAA\n+\nIII\n'),
             self.write('truncated.fq','@q/2\nAAA\n+\n'),self.write('bad.fq.gz','not gzip')]
        for n,b in enumerate(bad):
            one,two=self.out/f'o{n}.fq.gz',self.out/f't{n}.fq.gz'
            self.command(BIN,'delete-pairs','--names',ids,'--r1',good,'--r2',b,
                         '--out1',one,'--out2',two,success=False)
            self.assertFalse(one.exists());self.assertFalse(two.exists())

    def guard(self,name,code,*limits,success=False):
        out=self.out/name
        self.command(sys.executable,GUARD,'--output',out,'--poll','0.02','--min-free-gib','0',
                     *limits,'--',sys.executable,'-c',code,str(out),success=success)
        return json.loads(out.with_name(out.name+'.resources.json').read_text())

    def test_guard_success_and_no_overwrite(self):
        code='import pathlib,sys;p=pathlib.Path(sys.argv[1]);p.mkdir();(p/"SUCCESS").write_text("ok")'
        r=self.guard('success',code,success=True)
        self.assertEqual(r['status'],'success')
        self.command(sys.executable,GUARD,'--output',self.out/'success','--',sys.executable,'-c','pass',success=False)

    def test_guard_output_budget_and_hard_file_limit(self):
        code='import pathlib,sys,time;p=pathlib.Path(sys.argv[1]);p.mkdir();(p/"large").write_bytes(b"x"*65536);time.sleep(2)'
        r=self.guard('disk',code,'--max-gib',str(1024/1024**3))
        self.assertEqual(r['reason'],'output_size_limit')
        r=self.guard('file',code,'--max-file-gib',str(1024/1024**3))
        self.assertEqual(r['status'],'failed')
        self.assertLessEqual((self.out/'file/large').stat().st_size,1024)

    def test_guard_memory_limits(self):
        code=('import pathlib,sys,time;p=pathlib.Path(sys.argv[1]);p.mkdir();'
              'data=bytearray(128*1024*1024);time.sleep(3)')
        r=self.guard('rss',code,'--max-rss-gib','0.04')
        self.assertEqual(r['reason'],'rss_limit')
        self.assertGreater(r['peak_rss_bytes'],0.04*1024**3)
        r=self.guard('address',code,'--max-address-gib','0.06')
        self.assertEqual(r['status'],'failed')

    def test_guard_timeout_kills_descendants(self):
        child_code = 'import time,pathlib,sys;time.sleep(0.6);pathlib.Path(sys.argv[1]).write_text("leaked")'
        code = ('import pathlib,sys,subprocess,time;p=pathlib.Path(sys.argv[1]);p.mkdir();'
                f'subprocess.Popen([sys.executable,"-c",{child_code!r},str(p/"orphan")]);time.sleep(3)')
        r=self.guard('timeout',code,'--timeout','0.1')
        self.assertEqual(r['reason'],'wall_time_limit')
        # Give a leaked process an opportunity to leave its marker.
        import time
        time.sleep(0.7)
        self.assertFalse((self.out/'timeout/orphan').exists())


if __name__=='__main__':
    unittest.main(verbosity=2)
