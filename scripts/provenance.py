#!/usr/bin/env python3
"""Shared provenance helpers for portable NJU-seq tools."""
from datetime import datetime, timezone
import hashlib
import json

def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, data):
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def fingerprint(path):
    stat = path.stat()
    return dict(path=str(path), target=str(path.resolve()), bytes=stat.st_size,
                mtime_ns=stat.st_mtime_ns, inode=stat.st_ino, device=stat.st_dev)
