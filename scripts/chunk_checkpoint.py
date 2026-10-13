#!/usr/bin/env python3
"""Export only committed mRNA batches and verify immutable resume prerequisites."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

from provenance import fingerprint, save, sha, now

COUNTERS = ('sam_records', 'sam_selected', 'quality_passed', 'reference_unique',
            'coordinate_unique', 'read_groups', 'selected_reads')
VERSION = 2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def file_identity(path):
    # ctime also detects same-size edits followed by restoration of mtime.
    path = path.absolute()
    return dict(fingerprint(path), ctime_ns=path.stat().st_ctime_ns)


def bowtie_files():
    command = shutil.which('bowtie2')
    require(command is not None, 'bowtie2 is not on PATH')
    command = Path(command).resolve()
    backends = sorted({p.resolve() for pattern in ('bowtie2-align-s*', 'bowtie2-align-l*')
                       for p in command.parent.glob(pattern) if p.is_file()})
    require(backends, 'Bowtie2 alignment backends not found')
    return command, [command, *backends]


def prerequisites(native, index, exons, r1, r2):
    indexes = sorted(p for p in index.parent.glob(index.name+'.*.bt2*') if p.is_file())
    parts = ('1', '2', '3', '4', 'rev.1', 'rev.2')
    require(any({p.name for p in indexes} == {f'{index.name}.{part}.{ext}' for part in parts}
                for ext in ('bt2', 'bt2l')),
            'Expected six mRNA index files of one format (.bt2 or .bt2l)')
    command, files = bowtie_files()
    return dict(format_version=VERSION, index=str(index.resolve()),
                inputs=[file_identity(p) for p in (r1, r2)],
                references=[file_identity(p) for p in (exons, *indexes)],
                native=dict(identity=file_identity(native), sha256=sha(native)),
                bowtie2=str(command),
                tools=[dict(identity=file_identity(p), sha256=sha(p)) for p in files])


def prerequisite_digest(snapshot):
    return hashlib.sha256(json.dumps(snapshot, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def verify_prerequisites(snapshot, native, index, exons, r1, r2, *, hashes=True):
    require(snapshot.get('format_version') == VERSION, 'Unsupported prerequisite manifest')
    require(str(index.resolve()) == snapshot['index'], 'Checkpoint index mismatch')
    require([str(p.resolve()) for p in (r1, r2)] == [x['target'] for x in snapshot['inputs']],
            'Checkpoint FASTQ paths mismatch')
    require(str(exons.resolve()) == snapshot['references'][0]['target'], 'Checkpoint annotation mismatch')
    indexes=sorted(str(p.resolve()) for p in index.parent.glob(index.name+'.*.bt2*') if p.is_file())
    require(indexes==sorted(e['target'] for e in snapshot['references'][1:]),
            'Checkpoint index file set changed')
    for entry in snapshot['inputs'] + snapshot['references']:
        require(file_identity(Path(entry['path'])) == entry, 'Checkpoint prerequisite changed')
    # Resolve the executable that native align-select will actually invoke now.
    command, files = bowtie_files()
    require(str(command) == snapshot['bowtie2'] and
            [str(p) for p in files] == [t['identity']['target'] for t in snapshot['tools']],
            'Checkpoint active Bowtie2 differs from run start')
    # A supervisor may freeze an identical native executable in a new run package.
    if file_identity(native) != snapshot['native']['identity'] or hashes:
        require(sha(native) == snapshot['native']['sha256'], 'Checkpoint native binary changed')
    for item in [snapshot['native'], *snapshot['tools']]:
        path = Path(item['identity']['path'])
        require(file_identity(path) == item['identity'], 'Checkpoint tool prerequisite changed')
        if hashes:
            require(sha(path) == item['sha256'], 'Checkpoint tool content changed')


def ledger(rows, native, index, exons, digest):
    counts = dict.fromkeys(COUNTERS, 0)
    pairs = 0
    binaries = {}
    expected = sha(native)
    for number, row in enumerate(rows):
        require(row.get('prerequisites_sha256') == digest,
                'Checkpoint ledger is not bound to run-start prerequisites')
        require(row['batch'] == number and row['first_pair'] == pairs + 1 and
                isinstance(row['pairs'], int) and row['pairs'] > 0, 'Noncontiguous checkpoint ledger')
        a = row['argv']
        require(len(a) == 20 and a[1] == 'align-select', 'Unexpected checkpoint command')
        options = dict(zip(a[2::2], a[3::2]))
        require(set(options) == {'--mode','--index','--r1','--r2','--threads','--exons',
                                 '--output','--log','--stats'} and options['--mode'] == 'mrna',
                'Unexpected checkpoint alignment policy')
        require(Path(options['--index']).resolve() == index.resolve() and
                Path(options['--exons']).resolve() == exons.resolve(), 'Checkpoint reference differs')
        if a[0] not in binaries:
            binaries[a[0]] = sha(Path(a[0]))
        require(binaries[a[0]] == expected, 'Checkpoint native binary differs')
        stats = row['stats']
        require(not stats['sam_materialized'] and stats['screen_kind'] == 'none' and
                stats['screen_multimap'] is None, 'Checkpoint filtering policy differs')
        for key in COUNTERS:
            require(isinstance(stats[key], int) and stats[key] >= 0, 'Invalid checkpoint counter')
            counts[key] += stats[key]
        require(stats['selected_reads'] <= row['pairs'], 'Invalid accepted count')
        pairs += row['pairs']
    return pairs, counts


def export(args):
    manifest_path = args.chunks/'prerequisites.json'
    require(manifest_path.is_file(), 'Missing run-start prerequisites; legacy batches cannot be certified')
    snapshot = json.loads(manifest_path.read_text())
    verify_prerequisites(snapshot, args.native, args.index, args.exons, args.r1, args.r2)
    digest = prerequisite_digest(snapshot)
    progress_path = args.chunks / 'progress.json'
    progress = json.loads(progress_path.read_text())
    require(progress.get('prerequisites_sha256') == digest, 'Progress/prerequisites mismatch')
    n = progress['completed_batches'] if args.batches is None else args.batches
    require(0 < n <= progress['completed_batches'], 'Checkpoint must use committed batches')
    rows = []
    with (args.chunks/'batches.jsonl').open() as f:
        for _ in range(n):
            rows.append(json.loads(next(f)))
    pairs, counts = ledger(rows, args.native, args.index, args.exons, digest)
    if n == progress['completed_batches']:
        require(pairs == progress['completed_pairs'] and
                counts['selected_reads'] == progress['accepted_reads'], 'Progress/ledger mismatch')
    source_paths = [manifest_path, progress_path, args.chunks/'batches.jsonl', args.accepted, args.log]
    before = [fingerprint(p) for p in source_paths]
    args.output.mkdir(parents=True, exist_ok=False)
    with args.accepted.open('rb') as src, (args.output/'accepted.tsv').open('xb') as dst:
        for _ in range(counts['selected_reads']):
            line = src.readline()
            fields = line.rstrip(b'\r\n').split(b'\t')
            require(line.endswith(b'\n') and len(fields) == 5 and len(fields[4]) >= 10,
                    'Truncated or invalid checkpoint accepted prefix')
            dst.write(line)
    headers = 0
    with args.log.open('rb') as src, (args.output/'bowtie2.log').open('xb') as dst:
        for line in src:
            if line.startswith(b'=== batch '):
                if headers == n:
                    break
                m = re.fullmatch(rb'=== batch (\d+), pairs (\d+)-(\d+) ===\n', line)
                require(m is not None, 'Invalid checkpoint log boundary')
                row = rows[headers]
                require(tuple(map(int, m.groups())) ==
                        (headers, row['first_pair'], row['first_pair']+row['pairs']-1),
                        'Checkpoint log/ledger mismatch')
                headers += 1
            dst.write(line)
    require(headers == n, 'Incomplete checkpoint log')
    with (args.output/'batches.jsonl').open('x') as dst:
        for row in rows:
            dst.write(json.dumps(row)+'\n')
    require(before == [fingerprint(p) for p in source_paths],
            'Checkpoint source changed during export; stop producer first')
    verify_prerequisites(snapshot, args.native, args.index, args.exons, args.r1, args.r2)
    meta = dict(format_version=VERSION, created_utc=now(), completed_batches=n, completed_pairs=pairs,
                counts=counts, source_files=before, source_progress=progress,
                prerequisites=snapshot, prerequisites_sha256=digest,
                files_sha256={name:sha(args.output/name) for name in
                              ('accepted.tsv','bowtie2.log','batches.jsonl')})
    save(args.output/'checkpoint.json', meta)
    print(json.dumps(dict(output=str(args.output), batches=n, pairs=pairs,
                          accepted=counts['selected_reads'])), flush=True)


def load(checkpoint, native, index, exons, r1, r2):
    meta = json.loads((checkpoint/'checkpoint.json').read_text())
    require(meta.get('format_version') == VERSION,
            'Unsupported checkpoint: v2 run-start prerequisites are required; restart legacy runs')
    snapshot = meta['prerequisites']
    digest = prerequisite_digest(snapshot)
    require(digest == meta['prerequisites_sha256'], 'Checkpoint prerequisite digest mismatch')
    verify_prerequisites(snapshot, native, index, exons, r1, r2)
    require(set(meta['files_sha256']) == {'accepted.tsv','bowtie2.log','batches.jsonl'},
            'Invalid checkpoint file manifest')
    for name, digest in meta['files_sha256'].items():
        require(sha(checkpoint/name) == digest, f'Checkpoint content changed: {name}')
    rows = [json.loads(line) for line in (checkpoint/'batches.jsonl').read_text().splitlines()]
    pairs, counts = ledger(rows, native, index, exons, meta['prerequisites_sha256'])
    require(len(rows) == meta['completed_batches'] and pairs == meta['completed_pairs'] and
            counts == meta['counts'], 'Checkpoint metadata/ledger mismatch')
    with (checkpoint/'accepted.tsv').open('rb') as f:
        require(sum(1 for _ in f) == counts['selected_reads'], 'Checkpoint accepted count mismatch')
    return meta


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('chunks','accepted','log','native','index','exons','r1','r2','output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--batches', type=int)
    export(p.parse_args())
