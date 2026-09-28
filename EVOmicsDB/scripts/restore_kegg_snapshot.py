#!/usr/bin/env python3
"""Validate and restore missing snapshot tables from a user-authorized local copy.

No downloads are performed. Historical KEGG SHA256 values remain authoritative.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'resources/kegg_joint_snapshot'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def restore(destination, source=None, check_only=False):
    manifest = json.loads((destination / 'manifest.json').read_text())
    pending = []
    for name, expected in manifest['files_sha256'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Invalid manifest filename: ' + name)
        target = destination / relative
        if target.exists():
            if not target.is_file() or digest(target) != expected:
                raise ValueError('Existing snapshot SHA256 mismatch: ' + name)
        elif check_only or source is None:
            raise FileNotFoundError('External snapshot table required: ' + name)
        else:
            candidate = source / relative
            if not candidate.is_file() or digest(candidate) != expected:
                raise ValueError('Source table missing or SHA256 mismatch: ' + name)
            pending.append((candidate, target, expected))
    # Validate the entire set before copying anything. Never replace existing files.
    for candidate, target, expected in pending:
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as out, candidate.open('rb') as incoming:
            shutil.copyfileobj(incoming, out)
        if digest(target) != expected:
            target.unlink()
            raise ValueError('Copied table SHA256 mismatch: ' + target.name)
    return len(pending)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, help='User-authorized version-matched local tables')
    parser.add_argument('--destination', type=Path, default=DEFAULT)
    parser.add_argument('--check', action='store_true', help='Check only; never copy files')
    args = parser.parse_args()
    if not args.check and args.source is None:
        parser.error('--source is required unless --check is used')
    copied = restore(args.destination.resolve(), args.source.resolve() if args.source else None, args.check)
    print(f'PASS: all required snapshot hashes match; restored {copied} files.')


if __name__ == '__main__':
    main()
