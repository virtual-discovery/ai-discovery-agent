"""Download and verify the seven zips of Zenodo record 17951035 (CC BY 4.0).

    python3 -B benchmarks/construction/kagome135-correlation-origin/fetch.py [DEST]

Source: S. Enzner, G. Sangiovanni, '135 Kagome Wannier Hamiltonians', Zenodo (2025), DOI 10.5281/zenodo.17951035,
data of M. Crispino et al., 'Tunable electronic correlations in 135-kagome metals', arXiv:2512.22576.
Each zip is fetched through the record's per-file link, its md5 checked against the record, and its members
extracted into DEST/x and checked against members.md5. DEST defaults to
~/.cache/ai-theorist/kagome135-correlation-origin.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://zenodo.org/api/records/17951035/files/%s/content'
ZIPS = {'CsTiSb.zip': '3022b11da82d2a3dd91b57d2fe10d75b', 'CsCrSb.zip': 'ba28567cb06c56faaf9188c06ce5a1ca',
        'CsCrBi.zip': '95d1e17aaeb1c4b9392c096e27afd76f', 'DFT.zip': 'f8b8de95fa310d206ef01bdd27577501',
        'CsVSb.zip': '35a8b537aaf029def98e427f23785bed', 'CsCrBi_9.2.zip': '9c4a37b89801300f30eb2bd5d86167d3',
        'CsTiBi.zip': '1d92f3a6309e431f6debbbf1ed61df8f'}
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/kagome135-correlation-origin')
    dest.mkdir(parents=True, exist_ok=True)
    members = {}
    if (HERE/'members.md5').exists():
        members = dict(reversed(line.split(None, 1)) for line in (HERE/'members.md5').read_text().splitlines() if line.strip())
    for name, checksum in ZIPS.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(RECORD % name, target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
        with zipfile.ZipFile(target) as archive:
            for member in archive.namelist():
                if member.endswith('/'):
                    continue
                archive.extract(member, dest/'x')
                if members and md5(dest/'x'/member) != members[member]:
                    raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified {len(ZIPS)} zips and {len(members)} members in {dest}/x')


if __name__ == '__main__':
    main()
