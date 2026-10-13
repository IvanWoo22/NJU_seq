"""Reference construction: public-source formats, sequence checks and cache integrity."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from build_references import (build_arabidopsis_rrna, build_bacterial, download,
                              extract_sources, main as build_main)
from prepare_human_reference import fasta, references


class ReferenceBuilders(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='nju reference ')
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def gz(self, name, text):
        path = self.root / name
        with gzip.open(path, 'wt') as dst:
            dst.write(text)
        return path

    def inputs(self, species='human'):
        if species == 'human':
            tx, gene = 'ENST00001.2', 'ENSG00001.3'
            header = f'{tx}|{gene}|a|b|c|d|4|protein_coding|'
            parent, identity, feature = tx, tx, 'transcript'
        else:
            tx, gene = 'AT1G01010.1', 'AT1G01010'
            header = f'{tx} cdna gene:{gene} gene_biotype:protein_coding transcript_biotype:protein_coding'
            parent, identity, feature = f'transcript:{tx}', f'transcript:{tx}', 'mRNA'
        transcripts = self.gz('transcripts.fa.gz', f'>{header}\nACGT\n>noncoding\nTTTT\n')
        # Minus-strand exons at 1:1-2 and 1:7-8 reconstruct AC + GT -> ACGT.
        annotation = self.gz('annotation.gff3.gz',
            f'##gff-version 3\n1\t.\t{feature}\t1\t8\t.\t-\t.\tID={identity};Parent=gene:{gene}\n'
            f'1\t.\texon\t1\t2\t.\t-\t.\tParent={parent},other;gene_id={gene}\n'
            f'1\t.\texon\t7\t8\t.\t-\t.\tParent={parent};gene_id={gene}\n'
            '1\t.\texon\t3\t4\t.\t+\t.\tParent=other\n')
        genome = self.gz('genome.fa.gz', '>1\nACAAAAGT\n')
        return transcripts, annotation, genome

    def test_coding_selection_multi_parent_exons_and_no_overwrite(self):
        for species in ('human', 'arabidopsis'):
            with self.subTest(species=species):
                paths = extract_sources(species, *self.inputs(species), self.root / species)
                self.assertEqual(len(list(fasta(paths[0]))), 1)
                self.assertEqual(len(paths[1].read_text().splitlines()), 2)
                self.assertEqual(paths[2].read_text(), '>1\nACAAAAGT\n')
                self.assertIn('Parent=transcript:AT' if species == 'arabidopsis' else 'Parent=ENST',
                              paths[1].read_text())
                with self.assertRaisesRegex(ValueError, 'Refusing existing'):
                    extract_sources(species, *self.inputs(species), self.root / species)

    def test_missing_exons_and_mismatched_genes_fail(self):
        tx, ann, genome = self.inputs()
        with gzip.open(ann, 'rt') as src:
            original = src.read()
        for number, bad in enumerate([original.replace('ENSG00001.3', 'ENSG00009.3'),
                                     '\n'.join(x for x in original.splitlines() if '\texon\t' not in x)]):
            self.gz(ann.name, bad)
            out = self.root / f'bad{number}'
            with self.assertRaises(ValueError):
                extract_sources('human', tx, ann, genome, out)
            self.assertFalse(out.exists())

    def test_download_cache_checksums_and_offline_mode(self):
        upstream = self.root / 'upstream.gz'
        upstream.write_bytes(b'fixture source')
        source = dict(file='cached.gz', url=upstream.as_uri(),
                      sha256=hashlib.sha256(upstream.read_bytes()).hexdigest())
        cache = self.root / 'cache'
        with self.assertRaisesRegex(ValueError, 'Missing cached'):
            download(source, cache, offline=True)
        result = download(source, cache)
        upstream.unlink()
        self.assertEqual(download(source, cache, offline=True), result)
        result.write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'Checksum mismatch'):
            download(source, cache, offline=True)

    @unittest.skipUnless(shutil.which('bowtie2-build') and shutil.which('samtools'), 'reference tools required')
    def test_bacterial_taxonomy_and_arabidopsis_genome_validation(self):
        ssu = self.gz('ssu.fa.gz', '>s Bacteria;A;Lactobacillaceae;Genus\nACGUACGU\n'
                     '>skip Bacteria;A;NotLactobacillaceae;Genus\nTTTT\n')
        lsu = self.gz('lsu.fa.gz', '>l Bacteria;Enterobacteriaceae;Genus\nGGCC\n')
        bac = self.root / 'bacterial'
        build_bacterial([ssu, lsu], bac, 1)
        self.assertEqual(list(fasta(bac / 'bac_rna.fa')),
                         [('s Bacteria;A;Lactobacillaceae;Genus', 'ACGTACGT'),
                          ('l Bacteria;Enterobacteriaceae;Genus', 'GGCC')])
        mrna, exons, genome = extract_sources('arabidopsis', *self.inputs('arabidopsis'),
                                             self.root / 'plant source')
        rrna = self.root / 'rrna.fa'
        rrna.write_text('>18s\nACGTACGT\n')
        args = argparse.Namespace(out=self.root / 'plant ref', mrna=mrna, exons=exons,
                                  genome=genome, bac_index=bac / 'bac_rna', rrna=[rrna],
                                  copy_inputs=True, threads=1, label='fixture')
        references(args, species='arabidopsis',
                   prefixes={'mrna': 'ath_protein_coding', 'rrna': 'ath_rrna_total', 'bac': 'bac'})
        audit = json.loads((args.out / 'audit.json').read_text())
        self.assertEqual(audit['negative_strand_transcripts'], 1)
        self.assertEqual(audit['genome_sequence_mismatches'], [])
        self.assertTrue((args.out / 'ath_protein_coding.1.bt2').is_file())
        args.out = self.root / 'plant without bacteria'
        args.bac_index = None
        references(args, species='arabidopsis',
                   prefixes={'mrna': 'ath_protein_coding', 'rrna': 'ath_rrna_total', 'bac': 'bac'})
        self.assertEqual(json.loads((args.out / 'audit.json').read_text())['bacterial_index'], [])
        self.assertFalse(list(args.out.glob('bac*')))
        args.out = self.root / 'bad sequence'
        genome.write_text('>1\nTTTTTTTT\n')
        with self.assertRaisesRegex(ValueError, 'sequence mismatches'):
            references(args, species='arabidopsis')
        self.assertFalse(args.out.exists())

    @unittest.skipUnless(shutil.which('bowtie2-build') and shutil.which('samtools'), 'reference tools required')
    def test_human_build_without_bacterial_index(self):
        paths = self.inputs()
        manifest = self.root / 'sources.json'
        manifest.write_text(json.dumps({'human': dict(release='fixture', sources={
            key: dict(file=p.name, url=p.as_uri(), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for key, p in zip(('transcripts', 'annotation', 'genome'), paths)})}))
        out, source = self.root / 'public ref', self.root / 'source'
        args = ['build_references.py', 'human', '--without-bacterial', '--download-dir', str(self.root),
                '--source-out', str(source), '--out', str(out), '--threads', '1', '--offline', '--copy-inputs']
        with patch('build_references.MANIFEST', manifest), patch.object(sys, 'argv', args):
            build_main()
        self.assertEqual(json.loads((out / 'audit.json').read_text())['bacterial_index'], [])
        self.assertFalse(list(out.glob('bac*')))
        self.assertTrue((out / 'mrna.1.bt2').is_file())
        self.assertTrue((out / 'rrna.1.bt2').is_file())
        # The lower-level CLI also permits omitting --bac-index.
        custom = self.root / 'custom ref'
        result = subprocess.run([sys.executable, str(ROOT / 'scripts/prepare_human_reference.py'),
            '--mrna', str(source / 'hsa_protein_coding.fa'), '--exons', str(source / 'hsa_exon.info'),
            '--genome', str(source / 'GRCh38.primary_assembly.genome.fa'), '--out', str(custom),
            '--threads', '1'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads((custom / 'audit.json').read_text())['bacterial_index'], [])
        self.assertFalse(list(custom.glob('bac*')))

    def test_rrna_edits_and_public_genome_coordinates(self):
        canonical = self.root / 'canonical'
        canonical.mkdir()
        (canonical / '18s.fa').write_text('>18s\nACGT\n')
        bed = self.root / 'loci.bed'
        bed.write_text('Chr4\t1\t4\tstable_locus\n')
        digest = lambda text: hashlib.sha256(text.encode()).hexdigest()
        expected = '>18s\nATGT\n>stable_locus\nTGC\n'
        profile = self.root / 'profile.json'
        profile.write_text(json.dumps(dict(canonical={'18s': dict(
            input_sequence_sha256=digest('ACGT'), output_sequence_sha256=digest('ATGT'),
            edits=[dict(start=1, end=2, before='C', after='T')])},
            intervals=dict(file=bed.name, sha256=digest(bed.read_text())),
            sequence_count=2, fasta_sha256=digest(expected))))
        genome = self.gz('col.fa.gz', '>accession OriSeqID=Chr4\nATGC\n')
        out = self.root / 'total.fa'
        build_arabidopsis_rrna(genome, out, profile, canonical)
        self.assertEqual(out.read_text(), expected)
        bed.write_text('Chr4\t0\t4\twrong\n')
        with self.assertRaisesRegex(ValueError, 'interval checksum'):
            build_arabidopsis_rrna(genome, self.root / 'bad.fa', profile, canonical)


if __name__ == '__main__':
    unittest.main()
