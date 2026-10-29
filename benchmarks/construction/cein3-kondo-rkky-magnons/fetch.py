"""Refetch the pinned Zenodo record 10146787 (DOI 10.5281/zenodo.10146787, CC BY 4.0): 48 small .dat files, 1.0 MB.

    python3 -B fetch.py [DEST]     (default ~/.cache/ai-theorist/cein3-kondo-rkky-magnons/raw)
Each file is downloaded through the record's per-file API link and checked against members.json (md5).
"""
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/cein3-kondo-rkky-magnons/raw'
DEST.mkdir(parents=True, exist_ok=True)
for f in json.loads((HERE/'members.json').read_text())['files']:
    path = DEST/f['name']
    if not path.exists() or hashlib.md5(path.read_bytes()).hexdigest() != f['md5']:
        with urllib.request.urlopen(f['url'], timeout=120) as r:
            path.write_bytes(r.read())
    assert hashlib.md5(path.read_bytes()).hexdigest() == f['md5'], f['name']
print('verified', len(json.loads((HERE/'members.json').read_text())['files']), 'files in', DEST)
