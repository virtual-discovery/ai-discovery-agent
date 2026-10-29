"""Download and verify CWS_PRM.zip of Zenodo record 20131063 (CC BY 4.0).

    python3 -B benchmarks/construction/cu2wse4-orbital-edelstein/fetch.py [DEST]

Source: K. Nakazawa, T. Yamaguchi, A. Yamakage, "Orbital-related gyrotropic responses in Cu2WSe4 and chirality
indicator", Phys. Rev. Materials 10, 064603 (2026), arXiv:2512.16387; data DOI 10.5281/zenodo.20131063.

Downloads CWS_PRM.zip (24,610,495 bytes) through the record's per-file link, checks its md5, extracts the 34
members into DEST/x and checks each against members.md5. DEST defaults to
~/.cache/ai-theorist/cu2wse4-orbital-edelstein.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

URL = 'https://zenodo.org/api/records/20131063/files/CWS_PRM.zip/content'
MD5 = '9a50ca688ae85f597776e7d6e5931765'
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/cu2wse4-orbital-edelstein')
    dest.mkdir(parents=True, exist_ok=True)
    target = dest/'CWS_PRM.zip'
    if not target.exists() or md5(target) != MD5:
        urllib.request.urlretrieve(URL, target)
    if md5(target) != MD5:
        raise SystemExit('CWS_PRM.zip: md5 mismatch')
    members = dict(reversed(line.split(None, 1)) for line in (HERE/'members.md5').read_text().splitlines() if line.strip())
    with zipfile.ZipFile(target) as archive:
        for member in archive.namelist():
            archive.extract(member, dest/'x')
            if md5(dest/'x'/member) != members[member]:
                raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified CWS_PRM.zip and {len(members)} members in {dest}/x')


if __name__ == '__main__':
    main()
