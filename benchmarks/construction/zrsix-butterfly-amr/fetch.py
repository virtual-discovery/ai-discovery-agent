"""Download and verify the ZrSiX archive of Materials Cloud record tyan3-q3j67 (CC BY 4.0).

    python3 -B benchmarks/construction/zrsix-butterfly-amr/fetch.py [DEST]

Source: S. Zhang, O. V. Yazyev, "Magnetoresistance in ZrSiX (X = S, Se, Te) nodal-line semimetals",
Phys. Rev. B 112, 165147 (2025), arXiv:2504.12547; data DOI 10.24435/materialscloud:3q-r9.

Downloads 'ZrSiX data.zip' (8,378,201 bytes) and README.txt through the record's per-file links, checks their md5,
extracts every member except the VASP POTCARs (licensed, never staged) and checks each against members.md5.
DEST defaults to ~/.cache/ai-theorist/zrsix-butterfly-amr; members land in DEST/ZrSiX data/...
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/tyan3-q3j67/files'
FILES = {'data.zip': ('ZrSiX%20data.zip', 'cbc8f3896f69ebf7c1cf68ad9f1a1dc8'),
         'README.txt': ('README.txt', '199687afb593272f9f73fc873188a519')}
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/zrsix-butterfly-amr')
    dest.mkdir(parents=True, exist_ok=True)
    for name, (remote, checksum) in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{remote}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = dict(reversed(line.split(None, 1)) for line in (HERE/'members.md5').read_text().splitlines() if line.strip())
    n = 0
    with zipfile.ZipFile(dest/'data.zip') as archive:
        for member in archive.namelist():
            if member.endswith('/') or member.endswith('POTCAR'):
                continue
            archive.extract(member, dest)
            if md5(dest/member) != members[member]:
                raise SystemExit(f'{member}: md5 mismatch')
            n += 1
    print(f'verified 2 files and {n} members (POTCARs skipped) in {dest}')


if __name__ == '__main__':
    main()
