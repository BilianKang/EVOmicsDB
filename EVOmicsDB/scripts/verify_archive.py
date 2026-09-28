#!/usr/bin/env python3
"""Verify listed release files against SHA256SUMS.txt, without dependencies.

Run from any directory. Local outputs, restored inputs and Git metadata are not
part of the checksum list and do not affect verification of the release files.
"""
from pathlib import Path, PurePosixPath
import hashlib
import re
import sys

ROOT = Path(__file__).resolve().parents[1]


def verify(root: Path) -> tuple[int, list[str]]:
    root = root.resolve()
    errors = []
    names = set()
    checked = 0
    try:
        lines = (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        return 0, [f"Cannot read SHA256SUMS.txt: {exc}"]
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{64})  (.+)", line)
        if match is None:
            errors.append(f"Malformed checksum line {number}; expected SHA256, two spaces, relative path")
            continue
        expected, name = match.groups()
        relative = PurePosixPath(name)
        if (not relative.parts or relative.is_absolute() or ".." in relative.parts
                or str(relative) != name or "\\" in name):
            errors.append("Unsafe or noncanonical path: " + name)
            continue
        if name in names:
            errors.append("Duplicate checksum path: " + name)
            continue
        names.add(name)
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            errors.append("Path leaves package: " + name)
            continue
        if not path.is_file():
            errors.append("Missing: " + name)
            continue
        try:
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
        except OSError as exc:
            errors.append(f"Cannot read {name}: {exc}")
            continue
        if digest.hexdigest() != expected.lower():
            errors.append("Mismatch: " + name)
        checked += 1
    if not names:
        errors.append("SHA256SUMS.txt contains no valid file records")
    return checked, errors


def main() -> int:
    checked, errors = verify(ROOT)
    print("\n".join(errors) if errors else f"PASS: {checked} files match SHA-256.")
    return int(bool(errors))


if __name__ == "__main__":
    sys.exit(main())
