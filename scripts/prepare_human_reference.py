#!/usr/bin/env python3
"""Build human indexes from matching GENCODE-style transcripts, exons and genome."""
import argparse
import collections
import gzip
import hashlib
import json
import mmap
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def fingerprint(path, digest=True):
    st = path.stat()
    return dict(path=str(path.resolve()), bytes=st.st_size, mtime_ns=st.st_mtime_ns,
                **({'sha256': sha256(path)} if digest else {}))


def fasta(path, rna=False):
    """Read DNA reference FASTA; reject malformed headers and duplicate IDs."""
    name, chunks, seen = None, [], set()
    alphabet = frozenset('ACGTRYSWKMBDHVN' + ('U' if rna else ''))
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        for number, line in enumerate(f, 1):
            if line.startswith('>'):
                if name is not None:
                    if not chunks:
                        raise ValueError(f'{path}:{number}: empty FASTA record {name}')
                    yield name, ''.join(chunks).upper().replace('U', 'T')
                name, chunks = line[1:].strip(), []
                if not name:
                    raise ValueError(f'{path}:{number}: empty FASTA header')
                identity = name.split()[0]
                if identity in seen:
                    raise ValueError(f'{path}:{number}: duplicate FASTA ID {identity}')
                seen.add(identity)
            else:
                sequence = line.strip()
                if not sequence:
                    continue
                if '>' in sequence:
                    raise ValueError(f'{path}:{number}: FASTA header must start a new line')
                if name is None:
                    raise ValueError(f'{path}:{number}: sequence before FASTA header')
                if set(sequence.upper()) - alphabet:
                    raise ValueError(f'{path}:{number}: invalid DNA FASTA sequence character')
                chunks.append(sequence)
    if name is None:
        raise ValueError(f'{path}: empty FASTA file')
    if not chunks:
        raise ValueError(f'{path}: empty FASTA record {name}')
    yield name, ''.join(chunks).upper().replace('U', 'T')


def start_package(out):
    out = out.absolute()
    stage = out.with_name(out.name + '.partial')
    if out.exists() or out.is_symlink() or stage.exists() or stage.is_symlink():
        raise ValueError(f'Refusing existing output or staging directory: {out}')
    out.parent.mkdir(parents=True, exist_ok=True)
    stage.mkdir()
    return out, stage


def include_source(source, destination, copy_inputs):
    """Copy a source or link to it relative to the reference directory."""
    source = source.resolve()
    if copy_inputs:
        shutil.copyfile(source, destination)
    else:
        destination.symlink_to(os.path.relpath(source, destination.parent))


def references(args, species='human', prefixes=None):
    prefixes = prefixes or {'mrna': 'mrna', 'rrna': 'rrna', 'bac': 'bac'}
    out, stage = start_package(args.out)
    records = list(fasta(args.mrna))
    sequences = {}
    genes = {}
    for header, seq in records:
        tx, gene = transcript_identity(header, species)
        if tx in sequences:
            raise ValueError(f'Duplicate or unsupported transcript ID: {tx}')
        sequences[tx], genes[tx] = seq, gene
    del records
    exons = collections.defaultdict(set)
    coordinates = {}
    par_gene_aliases = set()
    with args.exons.open() as f:
        for line in f:
            chrom, a, b, strand, attributes = line.rstrip('\n').split('\t')
            attr = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
            tx = attr.get('Parent', '').removeprefix('transcript:')
            if tx not in sequences:
                continue
            if attr.get('gene_id') != genes[tx]:
                # Supplied GENCODE-style exon table uses the unsuffixed ENSG
                # attribute for chrY PAR transcripts; the FASTA adds _PAR_Y.
                if not (tx.endswith('_PAR_Y') and chrom == 'chrY'
                        and genes[tx].endswith('_PAR_Y')
                        and attr.get('gene_id') == genes[tx].removesuffix('_PAR_Y')):
                    raise ValueError(f'FASTA/annotation gene mismatch: {tx}')
                par_gene_aliases.add(tx)
            coord = (chrom, strand)
            if tx in coordinates and coordinates[tx] != coord:
                raise ValueError(f'Conflicting chromosome/strand: {tx}')
            coordinates[tx] = coord
            a, b = int(a), int(b)
            if a < 1 or b < a or strand not in ('+', '-'):
                raise ValueError(f'Invalid exon: {tx}')
            exons[tx].add((a, b))
    if sequences.keys() != exons.keys():
        raise ValueError('Some FASTA transcripts have no annotation')
    genome = args.genome.resolve()
    include_source(genome, stage / 'genome.fa', args.copy_inputs)
    with (stage / 'faidx.log').open('w') as log:
        subprocess.run(['samtools', 'faidx', str(stage/'genome.fa')],
                       stdout=log, stderr=log, check=True)
    fai = {}
    for line in (stage / 'genome.fa.fai').read_text().splitlines():
        name, length, offset, bases, width = line.split('\t')[:5]
        fai[name] = tuple(map(int, (length, offset, bases, width)))
    complement = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
    mismatches = []
    with genome.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as genome_map:
        def get_sequence(chrom, a, b):
            length, offset, bases, width = fai[chrom]
            if a < 1 or b > length:
                raise ValueError(f'Genome interval out of bounds: {chrom}:{a}-{b}')
            begin = offset + ((a-1)//bases)*width + (a-1)%bases
            end = offset + ((b-1)//bases)*width + (b-1)%bases + 1
            return genome_map[begin:end].replace(b'\n', b'').replace(b'\r', b'').decode().upper()
        with (stage/'exon.info').open('w') as dst:
            for tx in sorted(sequences):
                chrom, strand = coordinates[tx]
                blocks = sorted(exons[tx])
                if any(b1 >= a2 for (_, b1), (a2, _) in zip(blocks, blocks[1:])):
                    raise ValueError(f'Overlapping exons in transcript: {tx}')
                seq = ''.join(get_sequence(chrom, a, b) for a, b in blocks)
                if strand == '-':
                    seq = seq.translate(complement)[::-1]
                if seq != sequences[tx]:
                    mismatches.append([tx, len(sequences[tx]), len(seq),
                                       sum(a != b for a, b in zip(seq, sequences[tx]))])
                for a, b in blocks:
                    gene = genes[tx].removesuffix('_PAR_Y') if tx in par_gene_aliases else genes[tx]
                    parent = f'transcript:{tx}' if species == 'arabidopsis' else tx
                    dst.write(f'{chrom}\t{a}\t{b}\t{strand}\tParent={parent};gene_id={gene}\n')
    audit = dict(transcripts=len(sequences), fasta_gene_ids=len(set(genes.values())),
                 annotation_gene_ids=len({g.removesuffix('_PAR_Y') if tx in par_gene_aliases else g
                                          for tx, g in genes.items()}),
                 exon_rows=sum(map(len, exons.values())),
                 par_y_transcripts=sum(t.endswith('_PAR_Y') for t in sequences),
                 par_y_fasta_gene_aliases=len(par_gene_aliases),
                 negative_strand_transcripts=sum(s == '-' for c, s in coordinates.values()),
                 genome_sequence_mismatches=mismatches,
                 source=[fingerprint(p) for p in (args.mrna, args.exons, genome)],
                 annotation_scope='All supplied protein_coding transcripts; basic is not substituted',
                 rrna_source=[str(p.resolve()) for p in args.rrna],
                 reference_label=args.label,
                 input_storage='copy' if args.copy_inputs else 'relative_symlink',
                 species=species,
                 annotation_scope_note='Matching transcript FASTA, exon table and genome; source release is supplied by the caller')
    (stage/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    if mismatches:
        raise ValueError(f'{len(mismatches)} sequence mismatches; see {stage}/audit.json')
    # Normalize newlines and retain original record order, including pipe headers.
    with (stage/'mrna.fa').open('w') as dst:
        for header, sequence in fasta(args.mrna):
            dst.write(f'>{header}\n{sequence}\n')
    rrna_ids = set()
    with (stage/'rrna.fa').open('w') as dst:
        for source in args.rrna:
            for header, sequence in fasta(source):
                identity = header.split()[0]
                if identity in rrna_ids:
                    raise ValueError(f'Duplicate rRNA ID across files: {identity}')
                rrna_ids.add(identity)
                dst.write(f'>{header}\n{sequence}\n')
    parts = ('1', '2', '3', '4', 'rev.1', 'rev.2')
    extension = next((ext for ext in ('bt2', 'bt2l')
                      if all(Path(f'{args.bac_index}.{part}.{ext}').is_file() for part in parts)), None)
    if extension is None:
        raise ValueError('Bacterial index requires six matching .bt2 or .bt2l components')
    for part in parts:
        src = Path(f'{args.bac_index}.{part}.{extension}').resolve()
        include_source(src, stage/f'{prefixes["bac"]}.{part}.{extension}', args.copy_inputs)
    for name in ('rrna', 'mrna'):
        cmd = ['bowtie2-build', '--threads', str(args.threads),
               str(stage/f'{name}.fa'), str(stage/prefixes[name])]
        print('Building', name, flush=True)
        with (stage/f'{name}.build.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
    audit['rrna_references'] = [fingerprint(p) for p in args.rrna]
    audit['rrna_ids'] = sorted(rrna_ids)
    audit['bacterial_index'] = [fingerprint(p) for p in sorted(stage.glob(f'{prefixes["bac"]}.*.bt2*'))]
    audit['derived'] = [dict(fingerprint(p), name=p.name, path=str(out/p.name))
                        for p in sorted(stage.iterdir())
                        if p.is_file() and not p.is_symlink() and p.name != 'audit.json']
    (stage/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    stage.rename(out)
    print(json.dumps({k: audit[k] for k in ('transcripts', 'annotation_gene_ids', 'exon_rows',
                                          'par_y_transcripts', 'genome_sequence_mismatches')}))


def transcript_identity(header, species='human'):
    """Validate a selected protein-coding FASTA header and return transcript/gene IDs."""
    if species == 'human':
        fields = header.split('|')
        if len(fields) < 8 or fields[7] != 'protein_coding':
            raise ValueError(f'Unexpected protein coding header: {header}')
        tx, gene = fields[:2]
        tx_pattern, gene_pattern = r'ENST\d+\.\d+(?:_PAR_Y)?', r'ENSG\d+\.\d+(?:_PAR_Y)?'
    elif species == 'arabidopsis':
        fields = dict(item.split(':', 1) for item in header.split() if ':' in item)
        tx, gene = header.split()[0], fields.get('gene', '')
        if fields.get('transcript_biotype') != 'protein_coding':
            raise ValueError(f'Unexpected protein coding header: {header}')
        tx_pattern, gene_pattern = r'AT[1-5MC]G\d+\.\d+', r'AT[1-5MC]G\d+'
        if tx.rsplit('.', 1)[0] != gene:
            raise ValueError(f'Transcript/gene mismatch in FASTA: {header}')
    else:
        raise ValueError(f'Unsupported species: {species}')
    if not re.fullmatch(tx_pattern, tx) or not re.fullmatch(gene_pattern, gene):
        raise ValueError(f'Unsupported transcript/gene ID: {tx}, {gene}')
    return tx, gene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('mrna', 'exons', 'genome', 'bac-index', 'out'):
        parser.add_argument('--' + option, type=Path, required=True)
    parser.add_argument('--rrna', type=Path, action='append',
                        help='Repeat for each FASTA; default: repository data/hsa_rrna/{18s,28s,5-8s}.fa')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--copy-inputs', action='store_true',
                        help='Copy genome and bacterial index for a self-contained package; '
                             'default: relative links (move sources together with the package)')
    parser.add_argument('--label', default='unspecified release', help='Reference source/release; never inferred')
    args = parser.parse_args()
    args.rrna = args.rrna or [ROOT/'data/hsa_rrna'/f'{n}.fa' for n in ('18s', '28s', '5-8s')]
    if args.threads < 1:
        parser.error('threads must be positive')
    if args.genome.suffix == '.gz':
        parser.error('--genome requires an uncompressed, line-wrapped FASTA')
    for source in [args.mrna, args.exons, args.genome, *args.rrna]:
        if not source.is_file():
            parser.error(f'Missing reference input: {source}')
    references(args)


if __name__ == '__main__':
    main()
