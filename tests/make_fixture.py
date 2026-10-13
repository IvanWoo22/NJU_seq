#!/usr/bin/env python3
"""Rebuild the tiny synthetic fixture. ENST/ENSG identifiers here are invented."""
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parent / 'data/synthetic'


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(6096)
    sequences = {name: ''.join(rng.choices('ACGT', k=length))
                 for name, length in [('mrna', 500), ('rrna', 240), ('bac', 240)]}
    tx, gene = 'ENST90000000001.1', 'ENSG90000000001.1'
    header = f'{tx}|{gene}|synthetic|synthetic|TEST|TEST|500|protein_coding|'
    for name, sequence in sequences.items():
        title = header if name == 'mrna' else name
        (ROOT / f'{name}.fa').write_text(f'>{title}\n{sequence}\n')
    (ROOT / 'genome.fa').write_text('>chr1\n' + 'N' * 100 + sequences['mrna'] + 'N' * 100 + '\n')
    (ROOT / 'exon.info').write_text(f'chr1\t101\t600\t+\tParent={tx};gene_id={gene}\n')
    manifest = ['sample\trole\tr1\tr2\talignment\tumi']
    expected = {}
    for sample in ('NC', 'TR1', 'TR2', 'TR3'):
        n = [1, 10, 40] if sample == 'NC' else [60, 1, 10]
        # Keys are the 1-based transcript candidate positions (read end minus one).
        counts = dict(zip((119, 118, 249), n))
        inserts = []
        for candidate, number in counts.items():
            fragment = sequences['mrna'][candidate - 35:candidate + 1]
            inserts.extend([fragment] * number)
        for name, number in [('bac', 2), ('rrna', 3)]:
            inserts.extend([sequences[name][50:86]] * number)
        reads = [[], []]
        for i, insert in enumerate(inserts):
            reverse = insert.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            for mate, sequence in enumerate([insert + 'AGATCGGAAGAGCACA', reverse + 'GATCGTCGGACTGTAG']):
                reads[mate].append(f'@{sample}_{i:04d}/{mate + 1}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
        for mate in (1, 2):
            (ROOT / f'{sample}.R{mate}.fq').write_text(''.join(reads[mate - 1]))
        manifest.append(f'{sample}\t' + ('NC' if sample == 'NC' else 'TR') +
                        f'\t{sample}.R1.fq\t{sample}.R2.fq\t-\tnone')
        expected[sample] = dict(raw_pairs=len(inserts), bac_removed=2, host_removed=3,
                                accepted=sum(n), genome_end_counts={str(k + 100): v for k, v in counts.items()})
    (ROOT / 'samples.tsv').write_text('\n'.join(manifest) + '\n')
    (ROOT / 'truth.json').write_text(json.dumps(dict(kind='synthetic; not biological human data',
                                                   samples=expected, candidate_positions=[219]), indent=2) + '\n')


if __name__ == '__main__':
    main()
