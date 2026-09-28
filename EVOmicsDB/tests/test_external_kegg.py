"""Test snapshot restoration using artificial tables; no KEGG requests or data."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('snapshot_restore', ROOT / 'scripts/restore_kegg_snapshot.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    source = root / 'authorized'
    target = root / 'local'
    source.mkdir()
    target.mkdir()
    original = {'one.tsv': b'ID\tvalue\nTEST_A\t1\n', 'two.tsv': b'ID\tvalue\nTEST_B\t2\n'}
    manifest = {'files_sha256': {name: hashlib.sha256(data).hexdigest() for name, data in original.items()}}
    (target / 'manifest.json').write_text(json.dumps(manifest))
    for name, data in original.items():
        (source / name).write_bytes(data)
    try:
        module.restore(target, check_only=True)
    except FileNotFoundError:
        pass
    else:
        raise AssertionError('Missing external input was accepted')
    (source / 'two.tsv').write_bytes(b'wrong version')
    try:
        module.restore(target, source)
    except ValueError:
        pass
    else:
        raise AssertionError('Incorrect snapshot was accepted')
    assert not (target / 'one.tsv').exists(), 'Partial copying happened before complete validation'
    (source / 'two.tsv').write_bytes(original['two.tsv'])
    assert module.restore(target, source) == 2
    assert module.restore(target, check_only=True) == 0
    (target / 'one.tsv').write_bytes(b'local modification')
    try:
        module.restore(target, source)
    except ValueError:
        pass
    else:
        raise AssertionError('Unexpected local content was silently overwritten')

print('PASS: absent external inputs, wrong SHA256, pre-copy validation, successful restore and no silent overwrite.')
