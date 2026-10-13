#!/usr/bin/env python3
"""Checkpoint boundary/failure tests; real alignment parity is tested separately."""
import json
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from unittest.mock import patch

import align_chunked
import chunk_checkpoint
from provenance import fingerprint


class Checkpoints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='nju-checkpoint-')
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.root = Path(cls.tmp.name)
        print('CASE_ROOT', cls.root, flush=True)

    def setUp(self):
        self.folder = self.root/self._testMethodName
        self.folder.mkdir()
        self.bowtie = self.folder/'tools';self.bowtie.mkdir()
        for name in ('bowtie2','bowtie2-align-s','bowtie2-align-l'):
            p=self.bowtie/name;p.write_text('#!/bin/sh\nexit 0\n');p.chmod(0o755)
        env=patch.dict(os.environ,PATH=str(self.bowtie)+os.pathsep+os.environ['PATH'])
        env.start();self.addCleanup(env.stop)
        self.native = self.folder/'fake_native.py'
        self.native.write_text('''#!/usr/bin/env python3
import json,sys
from pathlib import Path
a=dict(zip(sys.argv[2::2],sys.argv[3::2]))
lines=Path(a['--r1']).read_text().splitlines()
n=len(lines)//4
with Path(a['--output']).open('w') as f:
 for i in range(n):
  f.write(lines[4*i][1:]+'\\tENST1.1\\t1\\t10=\\t'+lines[4*i+1]+'\\n')
Path(a['--log']).write_text('fake batch complete\\n')
s={k:n for k in ['sam_records','sam_selected','quality_passed','reference_unique','coordinate_unique','read_groups','selected_reads']}
s.update(sam_materialized=False,screen_kind='none',screen_multimap=None)
Path(a['--stats']).write_text(json.dumps(s))
''')
        self.native.chmod(0o755)
        self.index = self.folder/'index'
        for part in ('1','2','3','4','rev.1','rev.2'):
            (self.folder/f'index.{part}.bt2').write_text('fixture index')
        self.exons = self.folder/'exons'
        self.exons.write_text('fixture annotation')
        self.r1,self.r2 = [self.folder/f'R{m}.fq' for m in (1,2)]
        for f in (self.r1,self.r2):
            f.write_text(''.join(f'@r{i}\nACGTACGTAC\n+\nJJJJJJJJJJ\n' for i in range(3)))
        self.original = self.options('original')
        align_chunked.run(self.original)
        self.chunks = Path(str(self.original.output)+'.chunks')
        self.checkpoint = self.folder/'checkpoint'
        self.export_options = SimpleNamespace(chunks=self.chunks,accepted=self.original.output,
            log=self.original.log,native=self.native,index=self.index,exons=self.exons,
            r1=self.r1,r2=self.r2,output=self.checkpoint,batches=1)

    def options(self, name, checkpoint=None):
        p=self.folder/name;p.mkdir()
        return SimpleNamespace(native=self.native,index=self.index,exons=self.exons,
             r1=self.r1,r2=self.r2,output=p/'accepted.tsv',log=p/'bowtie2.log',
             stats=p/'selection.json',threads=8,chunk_pairs=2,resume_checkpoint=checkpoint)

    def load(self):
        return chunk_checkpoint.load(self.checkpoint,self.native,self.index,self.exons,self.r1,self.r2)

    def test_resume_preserves_prefix_exactly_once(self):
        chunk_checkpoint.export(self.export_options)
        resumed=self.options('resumed',self.checkpoint)
        align_chunked.run(resumed)
        self.assertEqual(resumed.output.read_bytes(),self.original.output.read_bytes())
        stats=json.loads(resumed.stats.read_text())
        self.assertEqual((stats['input_pairs'],stats['selected_reads'],stats['reused_pairs']),(3,3,2))
        rows=[json.loads(l) for l in Path(str(resumed.output)+'.chunks/batches.jsonl').read_text().splitlines()]
        self.assertEqual([x['first_pair'] for x in rows],[1,3])

    def test_mixed_index_formats_rejected_before_alignment(self):
        (self.folder / 'index.1.bt2').rename(self.folder / 'index.1.bt2l')
        options = self.options('mixed')
        with self.assertRaisesRegex(ValueError, 'one format'):
            align_chunked.run(options)
        self.assertFalse(options.output.exists())

    def test_native_path_without_directory_uses_current_directory(self):
        options = self.options('relative')
        options.native = Path(self.native.name)
        previous = Path.cwd()
        try:
            os.chdir(self.folder)
            align_chunked.run(options)
        finally:
            os.chdir(previous)
        self.assertEqual(options.output.read_bytes(), self.original.output.read_bytes())

    def test_all_input_already_committed(self):
        self.export_options.batches=None
        chunk_checkpoint.export(self.export_options)
        resumed=self.options('resumed',self.checkpoint)
        align_chunked.run(resumed)
        self.assertEqual(resumed.output.read_bytes(),self.original.output.read_bytes())
        self.assertEqual(json.loads(resumed.stats.read_text())['reused_pairs'],3)

    def test_uncommitted_partial_tail_not_reused(self):
        progress=json.loads((self.chunks/'progress.json').read_text())
        progress.update(completed_batches=1,completed_pairs=2,accepted_reads=2)
        (self.chunks/'progress.json').write_text(json.dumps(progress))
        first=(self.chunks/'batches.jsonl').read_text().splitlines()[0]
        (self.chunks/'batches.jsonl').write_text(first+'\n{"batch":')
        with self.original.output.open('a') as f:f.write('truncated uncommitted row')
        self.original.log.write_text('\n=== batch 0, pairs 1-2 ===\nfake batch complete\n\n=== batch broken')
        self.export_options.batches=None
        chunk_checkpoint.export(self.export_options)
        self.assertEqual(len((self.checkpoint/'accepted.tsv').read_text().splitlines()),2)
        self.assertNotIn('broken',(self.checkpoint/'bowtie2.log').read_text())
        self.load()

    def test_corrupt_checkpoint_rejected_without_publishing(self):
        chunk_checkpoint.export(self.export_options)
        with (self.checkpoint/'accepted.tsv').open('a') as f:f.write('bad\n')
        resumed=self.options('resumed',self.checkpoint)
        with self.assertRaisesRegex(ValueError,'content changed'):
            align_chunked.run(resumed)
        self.assertFalse(resumed.output.exists())
        self.assertFalse(resumed.stats.exists())

    def test_changed_input_rejected(self):
        chunk_checkpoint.export(self.export_options)
        with self.r1.open('a') as f:f.write('bad\n')
        with self.assertRaisesRegex(ValueError,'prerequisite changed'):self.load()

    def test_wrong_metadata_counter_rejected(self):
        chunk_checkpoint.export(self.export_options)
        f=self.checkpoint/'checkpoint.json';d=json.loads(f.read_text());d['counts']['selected_reads']+=1
        f.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,'metadata/ledger mismatch'):self.load()

    def test_missing_committed_accepted_row_rejected(self):
        self.original.output.write_text(self.original.output.read_text().splitlines()[0]+'\n')
        with self.assertRaisesRegex(ValueError,'accepted prefix'):
            chunk_checkpoint.export(self.export_options)
        self.assertFalse((self.checkpoint/'checkpoint.json').exists())

    def test_duplicate_id_at_resume_boundary_rejected(self):
        # The bad boundary is present before the run; only the valid prefix commits.
        for f in (self.r1,self.r2):f.write_text(f.read_text().replace('@r2\n','@r1\n'))
        interrupted=self.options('interrupted')
        with self.assertRaisesRegex(ValueError,'Adjacent duplicate'):
            align_chunked.run(interrupted)
        work=Path(str(interrupted.output)+'.chunks')
        self.export_options.chunks=work
        self.export_options.accepted=work/'accepted.partial'
        self.export_options.log=work/'bowtie2.partial.log'
        chunk_checkpoint.export(self.export_options)
        resumed=self.options('resumed',self.checkpoint)
        with self.assertRaisesRegex(ValueError,'Adjacent duplicate'):
            align_chunked.run(resumed)
        self.assertFalse(resumed.output.exists())

    def test_mate_mismatch_in_skipped_prefix_rejected(self):
        self.r2.write_text(self.r2.read_text().replace('@r0\n','@different\n'))
        with self.assertRaisesRegex(ValueError,'prerequisite changed'):
            chunk_checkpoint.export(self.export_options)

    def test_sequence_edit_before_export_rejected_even_with_restored_mtime(self):
        for p in (self.r1,self.r2):
            before=p.stat();p.write_text(p.read_text().replace('ACGTACGTAC','CCGTACGTAC',1))
            os.utime(p,ns=(before.st_atime_ns,before.st_mtime_ns))
        with self.assertRaisesRegex(ValueError,'prerequisite changed'):
            chunk_checkpoint.export(self.export_options)
        self.assertFalse(self.checkpoint.exists())

    def test_reference_edit_before_export_rejected(self):
        self.exons.write_text('different annotation')
        with self.assertRaisesRegex(ValueError,'prerequisite changed'):
            chunk_checkpoint.export(self.export_options)

    def test_active_bowtie_changed_on_path_rejected(self):
        chunk_checkpoint.export(self.export_options)
        other=self.folder/'other_tools';shutil.copytree(self.bowtie,other)
        with patch.dict(os.environ,PATH=str(other)+os.pathsep+os.environ['PATH']):
            with self.assertRaisesRegex(ValueError,'active Bowtie2 differs'):self.load()
            self.export_options.output=self.folder/'other_checkpoint'
            with self.assertRaisesRegex(ValueError,'active Bowtie2 differs'):
                chunk_checkpoint.export(self.export_options)

    def test_alias_to_same_bowtie_allowed(self):
        chunk_checkpoint.export(self.export_options)
        alias=self.folder/'alias';alias.mkdir();(alias/'bowtie2').symlink_to(self.bowtie/'bowtie2')
        with patch.dict(os.environ,PATH=str(alias)+os.pathsep+os.environ['PATH']):self.load()

    def test_large_index_backend_edit_rejected(self):
        chunk_checkpoint.export(self.export_options)
        (self.bowtie/'bowtie2-align-l').write_text('#!/bin/sh\nexit 1\n')
        with self.assertRaisesRegex(ValueError,'tool prerequisite changed'):self.load()

    def test_added_index_family_rejected(self):
        chunk_checkpoint.export(self.export_options)
        (self.folder/'index.1.bt2l').write_text('different index family')
        with self.assertRaisesRegex(ValueError,'index file set changed'):self.load()

    def test_identical_frozen_native_copy_allowed(self):
        chunk_checkpoint.export(self.export_options)
        copied=self.folder/'frozen-native';shutil.copy2(self.native,copied)
        resumed=self.options('resumed',self.checkpoint);resumed.native=copied
        align_chunked.run(resumed)
        self.assertEqual(resumed.output.read_bytes(),self.original.output.read_bytes())

    def test_input_changed_during_batch_never_commits_or_publishes(self):
        execute=align_chunked.subprocess.run
        def changed(*args,**kwargs):
            result=execute(*args,**kwargs)
            self.r1.write_text(self.r1.read_text().replace('ACGTACGTAC','CCGTACGTAC',1))
            return result
        new=self.options('changed')
        with patch.object(align_chunked.subprocess,'run',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'prerequisite changed'):align_chunked.run(new)
        self.assertFalse(new.output.exists());self.assertFalse(new.stats.exists())
        self.assertEqual(Path(str(new.output)+'.chunks/batches.jsonl').read_text(),'')

    def test_ledger_binding_required(self):
        p=self.chunks/'batches.jsonl';rows=p.read_text().splitlines();row=json.loads(rows[0])
        row.pop('prerequisites_sha256');rows[0]=json.dumps(row);p.write_text('\n'.join(rows)+'\n')
        with self.assertRaisesRegex(ValueError,'not bound to run-start'):
            chunk_checkpoint.export(self.export_options)

    def test_legacy_run_and_checkpoint_rejected(self):
        chunk_checkpoint.export(self.export_options)
        p=self.checkpoint/'checkpoint.json';meta=json.loads(p.read_text());meta['format_version']=1
        p.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError,'v2 run-start'):self.load()
        (self.chunks/'prerequisites.json').unlink()
        with self.assertRaisesRegex(ValueError,'Missing run-start'):
            chunk_checkpoint.export(self.export_options)


if __name__ == '__main__':
    unittest.main()
