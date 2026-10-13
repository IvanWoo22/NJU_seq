#!/usr/bin/env python3
"""Download pinned sources and build bacterial, human or Arabidopsis references."""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
from urllib.request import Request, urlopen

from prepare_human_reference import (ROOT, fasta, references, sha256, start_package,
                                     transcript_identity)

MANIFEST = ROOT / 'config/reference_sources.json'
FAMILIES = {'Enterobacteriaceae', 'Lactobacillaceae', 'Streptococcaceae'}


def download(source, directory, offline=False):
    """Reuse only verified downloads; publish new files after checking their digest."""
    path = directory / source['file']
    if path.exists():
        if sha256(path) != source['sha256']:
            raise ValueError(f'Checksum mismatch: {path}; use a clean download directory')
        print(f'Using {path}', flush=True)
        return path
    if offline:
        raise ValueError(f'Missing cached download: {path}')
    directory.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + '.partial')
    for attempt in range(3):
        print(f'Downloading {source["url"]} (attempt {attempt + 1}/3)', flush=True)
        try:
            request = Request(source['url'], headers={'User-Agent': 'NJU-seq reference builder'})
            with urlopen(request, timeout=60) as response, partial.open('wb') as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            if sha256(partial) != source['sha256']:
                raise ValueError(f'Checksum mismatch: {partial}')
            partial.rename(path)
            return path
        except (OSError, ValueError) as error:
            if attempt == 2:
                raise
            print(f'{error}; retrying', flush=True)
            time.sleep(2)


def source_record(path, definition):
    return dict(file=path.name, url=definition['url'], sha256=sha256(path),
                bytes=path.stat().st_size)


def write_sources(out, profile, paths, species):
    report = dict(release=profile['release'],
                  downloads=[source_record(paths[key], source)
                             for key, source in profile['sources'].items()],
                  recipe_sha256=sha256(MANIFEST))
    if species == 'arabidopsis':
        report['rrna_definitions'] = [dict(file=p.name, sha256=sha256(p)) for p in (
            ROOT / 'config/arabidopsis_rrna.json', ROOT / 'config/arabidopsis_rrna.intervals.bed')]
    (out / 'sources.json').write_text(json.dumps(report, indent=2) + '\n')


def attrs(text):
    return dict(item.split('=', 1) for item in text.split(';') if '=' in item)


def extract_sources(species, transcripts, annotation, genome, out):
    """Select coding transcripts and convert their GFF3 exons to five columns."""
    out, stage = start_package(out)
    prefix = 'hsa' if species == 'human' else 'ath'
    mrna_name = f'{prefix}_protein_coding.fa'
    exons_name = f'{prefix}_exon.info'
    genome_name = 'GRCh38.primary_assembly.genome.fa' if species == 'human' else 'ath.fa'
    genes = {}
    with (stage / mrna_name).open('w') as dst:
        for header, sequence in fasta(transcripts):
            coding = (len(header.split('|')) >= 8 and header.split('|')[7] == 'protein_coding'
                      if species == 'human' else
                      'transcript_biotype:protein_coding' in header.split())
            if not coding:
                continue
            tx, gene = transcript_identity(header, species)
            if tx in genes:
                raise ValueError(f'Duplicate transcript: {tx}')
            genes[tx] = gene
            dst.write(f'>{header}\n{sequence}\n')
    if not genes:
        raise ValueError('No protein_coding transcripts found')
    annotated, parents, exon_count = set(), {}, 0
    opener = gzip.open if annotation.suffix == '.gz' else open
    with opener(annotation, 'rt') as src, (stage / exons_name).open('w') as dst:
        for number, line in enumerate(src, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9:
                raise ValueError(f'{annotation}:{number}: expected nine GFF3 columns')
            chrom, _, feature, start, end, _, strand, _, attributes = row
            info = attrs(attributes)
            identity = info.get('ID', '').removeprefix('transcript:')
            if identity in genes and feature in ('transcript', 'mRNA'):
                parent = info.get('Parent', '').removeprefix('gene:')
                if parent != genes[identity]:
                    raise ValueError(f'GFF3 transcript/gene mismatch: {identity}')
                if identity in parents and parents[identity] != parent:
                    raise ValueError(f'Conflicting GFF3 parent: {identity}')
                parents[identity] = parent
            if feature != 'exon':
                continue
            for parent in info.get('Parent', '').split(','):
                tx = parent.removeprefix('transcript:')
                if tx not in genes:
                    continue
                if int(start) < 1 or int(end) < int(start) or strand not in ('+', '-'):
                    raise ValueError(f'{annotation}:{number}: invalid exon coordinates')
                gene = info.get('gene_id', genes[tx])
                par_alias = (tx.endswith('_PAR_Y') and chrom == 'chrY'
                             and gene == genes[tx].removesuffix('_PAR_Y'))
                if gene != genes[tx] and not par_alias:
                    raise ValueError(f'GFF3 exon/gene mismatch: {tx}')
                parent = f'transcript:{tx}' if species == 'arabidopsis' else tx
                dst.write(f'{chrom}\t{start}\t{end}\t{strand}\tParent={parent};gene_id={gene}\n')
                annotated.add(tx)
                exon_count += 1
    if annotated != genes.keys() or parents.keys() != genes.keys():
        raise ValueError('Some selected transcripts have no matching GFF3 transcript/exon records')
    with gzip.open(genome, 'rb') as src, (stage / genome_name).open('wb') as dst:
        shutil.copyfileobj(src, dst)
    report = dict(species=species, transcripts=len(genes), exon_rows=exon_count,
                  selection='transcript biotype = protein_coding',
                  inputs=[dict(file=p.name, sha256=sha256(p))
                          for p in (transcripts, annotation, genome)])
    (stage / 'sources.json').write_text(json.dumps(report, indent=2) + '\n')
    stage.rename(out)
    return out / mrna_name, out / exons_name, out / genome_name


def build_bacterial(sources, out, threads, expected=None):
    out, stage = start_package(out)
    counts, seen = [], set()
    with (stage / 'bac_rna.fa').open('w') as dst:
        for source in sources:
            count = 0
            for header, sequence in fasta(source, rna=True):
                taxonomy = header.split(maxsplit=1)
                if len(taxonomy) < 2 or not FAMILIES.intersection(taxonomy[1].split(';')):
                    continue
                identity = taxonomy[0]
                if identity in seen:
                    raise ValueError(f'Duplicate bacterial rRNA ID: {identity}')
                seen.add(identity)
                dst.write(f'>{header}\n{sequence}\n')
                count += 1
            counts.append(count)
    if not seen:
        raise ValueError('No bacterial sequences matched the selected families')
    report = dict(families=sorted(FAMILIES), source_counts=counts, sequences=len(seen),
                  fasta_sha256=sha256(stage / 'bac_rna.fa'),
                  inputs=[dict(file=p.name, sha256=sha256(p)) for p in sources])
    if expected and (counts != expected['source_counts']
                     or report['fasta_sha256'] != expected['fasta_sha256']):
        raise ValueError('Bacterial reference differs from the pinned sequence profile')
    with (stage / 'bowtie2-build.log').open('w') as log:
        subprocess.run(['bowtie2-build', '--threads', str(threads),
                        str(stage / 'bac_rna.fa'), str(stage / 'bac_rna')],
                       stdout=log, stderr=log, check=True)
    (stage / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    stage.rename(out)
    print(json.dumps(report | {'inputs': [p.name for p in sources]}), flush=True)


def build_arabidopsis_rrna(genome, out, profile_path=ROOT / 'config/arabidopsis_rrna.json',
                           canonical_dir=ROOT / 'data/ath_rrna'):
    """Recreate the Col_NC31/Col_TR31 rRNA profile from small edits and public loci."""
    profile = json.loads(profile_path.read_text())
    bed = profile_path.parent / profile['intervals']['file']
    if sha256(bed) != profile['intervals']['sha256']:
        raise ValueError('rRNA interval checksum mismatch')
    records = []
    for identity, item in profile['canonical'].items():
        raw = list(fasta(canonical_dir / f'{identity}.fa'))
        if len(raw) != 1 or raw[0][0].split()[0] != identity:
            raise ValueError(f'Expected one canonical rRNA record: {identity}')
        sequence = raw[0][1]
        if hashlib.sha256(sequence.encode()).hexdigest() != item['input_sequence_sha256']:
            raise ValueError(f'Canonical rRNA input checksum mismatch: {identity}')
        chunks, cursor = [], 0
        for edit in item['edits']:
            a, b = edit['start'], edit['end']
            if not (cursor <= a <= b <= len(sequence)) or sequence[a:b] != edit['before']:
                raise ValueError(f'Invalid rRNA sequence edit: {identity}')
            chunks.extend((sequence[cursor:a], edit['after']))
            cursor = b
        sequence = ''.join(chunks) + sequence[cursor:]
        if hashlib.sha256(sequence.encode()).hexdigest() != item['output_sequence_sha256']:
            raise ValueError(f'Canonical rRNA output checksum mismatch: {identity}')
        records.append((identity, sequence))
    intervals = collections.defaultdict(list)
    with bed.open() as src:
        for line in src:
            chrom, start, end, identity = line.rstrip('\n').split('\t')
            intervals[chrom].append((int(start), int(end), identity))
    extracted = {}
    for header, sequence in fasta(genome):
        fields = dict(x.split('=', 1) for x in header.split() if '=' in x)
        chrom = fields.get('OriSeqID', header.split()[0])
        if chrom in extracted:
            raise ValueError(f'Duplicate genome chromosome: {chrom}')
        if chrom not in intervals:
            continue
        pieces = []
        for a, b, identity in intervals[chrom]:
            if not 0 <= a < b <= len(sequence):
                raise ValueError(f'Invalid rRNA interval: {chrom}:{a}-{b}')
            pieces.append((identity, sequence[a:b]))
        extracted[chrom] = pieces
    if intervals.keys() != extracted.keys():
        raise ValueError('rRNA interval chromosome missing from genome')
    for chrom in intervals:
        records.extend(extracted[chrom])
    if len(records) != profile['sequence_count']:
        raise ValueError('Unexpected rRNA sequence count')
    with out.open('x') as dst:
        for header, sequence in records:
            dst.write(f'>{header}\n{sequence}\n')
    if sha256(out) != profile['fasta_sha256']:
        raise ValueError('rRNA sequence checksum mismatch; check the Col-PEK genome version')
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('species', choices=('bacterial', 'human', 'arabidopsis'))
    parser.add_argument('--download-dir', type=Path, default=ROOT / 'references/downloads')
    parser.add_argument('--out', type=Path, help='New reference directory')
    parser.add_argument('--source-out', type=Path, help='New host source directory')
    parser.add_argument('--bac-index', type=Path, default=ROOT / 'references/bacterial/bac_rna')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--copy-inputs', action='store_true',
                        help='Copy the host genome and bacterial index into the output')
    parser.add_argument('--offline', action='store_true', help='Require verified cached downloads')
    parser.add_argument('--download-only', action='store_true', help='Download sources without building')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    manifest = json.loads(MANIFEST.read_text())
    profile = manifest[args.species]
    args.out = args.out or ROOT / 'references' / ('human_rebuilt' if args.species == 'human' else args.species)
    args.source_out = args.source_out or ROOT / 'references' / f'{args.species}_source'
    if not args.download_only:
        outputs = [args.out] + ([] if args.species == 'bacterial' else [args.source_out])
        for path in outputs:
            if any(p.exists() or p.is_symlink() for p in (path, path.with_name(path.name + '.partial'))):
                parser.error(f'Refusing existing output or staging directory: {path}')
        for tool in ['bowtie2-build'] + ([] if args.species == 'bacterial' else ['samtools']):
            if not shutil.which(tool):
                parser.error(f'{tool} is not on PATH; activate the nju-seq environment')
        if args.species != 'bacterial':
            parts = ('1', '2', '3', '4', 'rev.1', 'rev.2')
            if not any(all(Path(f'{args.bac_index}.{p}.{ext}').is_file() for p in parts)
                       for ext in ('bt2', 'bt2l')):
                parser.error('Build the bacterial index first with: build_references.py bacterial')
    paths = {key: download(source, args.download_dir, args.offline)
             for key, source in profile['sources'].items()}
    if args.download_only:
        return
    if args.species == 'bacterial':
        build_bacterial([paths['ssu'], paths['lsu']], args.out, args.threads, profile['expected'])
        write_sources(args.out, profile, paths, args.species)
        return
    args.mrna, args.exons, args.genome = extract_sources(
        args.species, paths['transcripts'], paths['annotation'], paths['genome'], args.source_out)
    args.label = profile['release']
    if args.species == 'human':
        args.rrna = [ROOT / 'data/hsa_rrna' / f'{n}.fa' for n in ('18s', '28s', '5-8s')]
        prefixes = None
    else:
        args.rrna = [build_arabidopsis_rrna(paths['rrna_genome'], args.source_out / 'ath_rrna_total.fa')]
        prefixes = {'mrna': 'ath_protein_coding', 'rrna': 'ath_rrna_total', 'bac': 'bac'}
    references(args, species=args.species, prefixes=prefixes)
    write_sources(args.out, profile, paths, args.species)


if __name__ == '__main__':
    main()
