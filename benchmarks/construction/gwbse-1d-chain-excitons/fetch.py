"""Download and verify the Materials Cloud archive 42t0v-jam38 (Grillo et al., ACS Nano 20, 2664 (2026)).

    python3 benchmarks/construction/gwbse-1d-chain-excitons/fetch.py [DEST]

Downloads 1D_ACS_Nano_Exciton.zip (3,319,010 bytes) and README.txt through the archive's per-file
links, checks their md5, extracts the members listed in members.md5 and checks each member's md5.
DEST defaults to ~/.cache/ai-theorist/gwbse-1d-chain-excitons.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/42t0v-jam38/files'
FILES = {'1D_ACS_Nano_Exciton.zip': 'a3cdbff7bbca1465fa0d2a8ec75ccec4', 'README.txt': '781fc14a3db2c37c787b60dac65d45cf'}
HERE = Path(__file__).resolve().parent


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/gwbse-1d-chain-excitons')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = [line.split(None, 1) for line in (HERE/'members.md5').read_text().splitlines() if line.strip()]
    with zipfile.ZipFile(dest/'1D_ACS_Nano_Exciton.zip') as archive:
        for checksum, member in members:
            archive.extract(member, dest)
            if md5(dest/member) != checksum:
                raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified 2 files and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
