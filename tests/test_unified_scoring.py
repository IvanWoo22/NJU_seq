#!/usr/bin/env python3
"""Exercise both native score modes, the independent Perl oracles and output semantics."""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from score_checks import rows,check_mrna,check_rrna

ROOT=Path(__file__).resolve().parents[1]
BIN=Path(os.environ.get('NJU_NATIVE',ROOT/'target/release/nju-seq'))
PERL=os.environ.get('NJU_PERL', 'perl')


class UnifiedScore(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.root=Path(temp.name);self.serial=0

    def counts(self,peaks=(100,100,100),nc_peak=1,length=120,position=60,minus=False,empty_nc=False):
        signals=[[0]+[1]*(length-1) for _ in range(4)]
        if empty_nc:signals[0]=[0]*length
        else:signals[0][position]=nc_peak
        for i,value in enumerate(peaks,1):signals[i][position]=value
        return signals,self.write(signals,minus)

    def write(self,signals,minus=False):
        folder=self.root/f'inputs{self.serial}';self.serial+=1;folder.mkdir();paths={'mrna':[],'rrna':[]}
        for name,values in zip(('NC','TR1','TR2','TR3'),signals):
            d=folder/name;d.mkdir();m=d/'mrna.tsv';r=d/'rrna.tsv'
            m.write_text(''.join(f'1\t{300-k if minus else k}\t{"-" if minus else "+"}\tA\tA\tA\tG1\t0\t{values[k]}\n'
                                 for k in range(1,len(values)) if values[k]>0))
            r.write_text(''.join(f'{i+1}\tA\t0\t{n}\n' for i,n in enumerate(values)))
            paths['mrna'].append(m);paths['rrna'].append(r)
        return paths

    def execute(self,mode,paths,extra=(),success=True):
        self.serial+=1;out=self.root/f'out{self.serial}.tsv';summary=out.with_suffix('.json')
        argv=[str(BIN),'score','--mode',mode,'--nc',str(paths[0])]
        for p in paths[1:]:argv+=['--tr',str(p)]
        if mode=='rrna':argv+=['--reference-id','18s']
        result=subprocess.run([*argv,'--output',str(out),'--summary',str(summary),*extra],capture_output=True,text=True)
        self.assertEqual(result.returncode==0,success,result.stderr)
        if not success:self.assertFalse(out.exists());self.assertFalse(summary.exists())
        return out,summary,result

    def oracle(self,mode,inputs):
        self.serial+=1;out=self.root/f'perl{self.serial}.tsv'
        result=subprocess.run([PERL,str(ROOT/f'{mode}_analysis/score.pl'),*map(str,inputs)],capture_output=True,text=True,check=True)
        out.write_text(result.stdout);return out

    def test_same_physical_signal_has_identical_shared_schema_and_score(self):
        _,paths=self.counts();outputs={}
        for mode in ('mrna','rrna'):
            outputs[mode]=self.execute(mode,paths[mode])[0]
        self.assertEqual(rows(outputs['mrna'])[0],rows(outputs['rrna'])[0])
        a=rows(outputs['mrna'])[1][0];b=next(r for r in rows(outputs['rrna'])[1] if r['Position']=='60')
        for row in (a,b):
            self.assertEqual(row['TotalScore'],'297');self.assertEqual(row['TotalFoldChange'],'165')
            self.assertEqual(row['AveFC'],'55');self.assertEqual(row['TR1_count'],'100')
            self.assertEqual(row['TR1_adjacent_count'],'1');self.assertEqual(row['NC_count(NC)'],'1')
        check_mrna(outputs['mrna'],self.oracle('mrna',paths['mrna']))
        check_rrna(outputs['rrna'],self.oracle('rrna',paths['rrna']),paths['rrna'],'18s')

    def test_negative_strand_candidate_and_adjacent_counts(self):
        _,paths=self.counts(minus=True)
        out,_,_=self.execute('mrna',paths['mrna']);row=rows(out)[1][0]
        self.assertEqual((row['Position'],row['Direction'],row['TR1_adjacent_count']),('240','-','1'))
        self.assertEqual(row['TotalScore'],'297')
        check_mrna(out,self.oracle('mrna',paths['mrna']))

    def test_zero_nc_total_never_calls_in_either_mode_or_legacy(self):
        _,paths=self.counts(empty_nc=True)
        for mode in ('mrna','rrna'):
            with self.subTest(mode=mode):
                out,summary,_=self.execute(mode,paths[mode],['--rows','all'])
                self.assertEqual(json.loads(summary.read_text())['joint_evaluation'],'normalization_unavailable')
                for row in rows(out)[1]:
                    self.assertEqual(row['Candidate'],'false');self.assertEqual(row['EvaluationStatus'],'normalization_unavailable')
                    self.assertEqual(row['TotalFoldChange'],row['TotalScore']);self.assertEqual(row['TotalScore'],'NA')
                if mode=='rrna':check_rrna(out,self.oracle(mode,paths[mode]),paths[mode],'18s')
                legacy,legacy_summary,_=self.execute(mode,paths[mode],['--format','legacy'])
                meaning=json.loads(legacy_summary.read_text())['count_meaning']
                self.assertEqual(meaning.startswith('legacy count'),mode=='rrna')
                _,expected=rows(self.oracle(mode,paths[mode]));_,actual=rows(legacy)
                self.assertEqual(actual,expected)

    def test_zero_overlap_for_one_tr_suppresses_joint_mrna(self):
        signals,_=self.counts(nc_peak=0);signals[3]=[0]*120;signals[3][60]=100;paths=self.write(signals)
        out,summary,_=self.execute('mrna',paths['mrna'])
        self.assertEqual(rows(out)[1],[]);self.assertEqual(rows(self.oracle('mrna',paths['mrna']))[1],[])
        info=json.loads(summary.read_text());self.assertEqual(info['samples'][-1]['status'],'zero_nc_total')
        self.assertEqual(info['joint_evaluation'],'normalization_unavailable')

    def test_empty_tr_invalidates_rrna_joint_even_if_others_are_strong(self):
        signals,_=self.counts(peaks=(200,200,200));signals[3]=[0]*120;paths=self.write(signals)
        out,summary,_=self.execute('rrna',paths['rrna'])
        self.assertEqual(json.loads(summary.read_text())['samples'][-1]['status'],'zero_tr_total')
        check_rrna(out,self.oracle('rrna',paths['rrna']),paths['rrna'],'18s')
        self.assertFalse(any(r['Candidate']=='true' for r in rows(out)[1]))

    def test_zero_site_nc_keeps_explicit_policy_difference(self):
        _,paths=self.counts(peaks=(30,30,30),nc_peak=0)
        for mode,score,numerator,called in [('mrna','87','1','false'),('rrna','90','0','true')]:
            out,_,_=self.execute(mode,paths[mode],['--rows','all'])
            row=next(r for r in rows(out)[1] if r['Position']=='60')
            self.assertEqual((row['TotalScore'],row['NC_score_numerator'],row['Candidate']),(score,numerator,called))
            self.assertEqual(row['EvaluationStatus'],'evaluated')

    def test_rrna_edges_and_short_references(self):
        _,paths=self.counts()
        out,_,_=self.execute('rrna',paths['rrna']);data={int(r['Position']):r for r in rows(out)[1]}
        for p,status in [(20,'boundary_excluded'),(21,'evaluated'),(100,'evaluated'),(101,'boundary_excluded')]:
            self.assertEqual(data[p]['EvaluationStatus'],status)
        for length in (1,2,20,40):
            paths=self.write([[1]*length for _ in range(4)])
            out,_,_=self.execute('rrna',paths['rrna']);data=rows(out)[1]
            self.assertEqual(len(data),length-1)
            self.assertTrue(all(r['EvaluationStatus']=='boundary_excluded' for r in data))

    def test_invalid_rrna_tables_publish_neither_output(self):
        for bad in ('2\tA\t0\t1\n','1\tC\t0\t1\n','1\tA\t0\t-1\n','1\tA\t0\t1\n1\tA\t0\t1\n',''):
            _,paths=self.counts();paths['rrna'][1].write_text(bad)
            self.execute('rrna',paths['rrna'],success=False)

    def test_mode_and_format_errors_publish_neither_output(self):
        _,paths=self.counts()
        for extra in (['--format','unknown'],['--rows','bad'],['--reference-id','18s'],['--format','legacy','--rows','all']):
            self.execute('mrna',paths['mrna'],extra,success=False)
        self.execute('wrong',paths['mrna'],success=False)

    def test_duplicate_sample_names_rejected_in_unified_output(self):
        _,paths=self.counts()
        self.execute('mrna',[paths['mrna'][0],paths['mrna'][1],paths['mrna'][1]],success=False)

    def test_joint_threshold_still_allows_one_strong_tr(self):
        _,paths=self.counts(peaks=(100,5,5))
        for mode in ('mrna','rrna'):
            out,_,_=self.execute(mode,paths[mode],['--rows','candidates'])
            row=next(r for r in rows(out)[1] if r['Position']=='60')
            self.assertEqual(row['TotalScore'],'107');self.assertEqual(row['Candidate'],'true')


if __name__=='__main__':unittest.main()
