"""Download and verify the BaB2 analysis files of Materials Cloud Archive record 0kqen-w8h90.

    python3 -B benchmarks/construction/bab2-hole-doping-two-gap/fetch.py [DEST]

Source: Bercx, Ponce, Zhang, Trezza, Ghorbani Ghezeljehmeidan, Bastonero, Qiao, von Rohr, Pizzi, Chiavazzo,
Marzari, "Charting the landscape of Bardeen-Cooper-Schrieffer superconductors in experimentally known
compounds", PRX Energy 4, 033012 (2025); archive DOI 10.24435/materialscloud:x2-ap (CC BY 4.0).

Downloads bab2-analysis.zip (28,310,115 bytes), README.md and files_description.md through the archive's
per-file links, checks their md5, extracts the 381 members listed in members.md5 (the large Wannier,
ph.out and supercell force-constant files are not extracted) and checks each member's md5.
DEST defaults to ~/.cache/ai-theorist/bab2-hole-doping-two-gap; members land in DEST/x/.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/0kqen-w8h90/files'
FILES = {'bab2-analysis.zip': 'a63a59c6c913b4cab7d5050eb0f1710c',
         'README.md': '57c5c2ef290197ab4035051e52fcc5d1',
         'files_description.md': '405192de6ce3f5f9b25d2f376ba92694'}
HERE = Path(__file__).resolve().parent


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/bab2-hole-doping-two-gap')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = [line.split(None, 1) for line in (HERE/'members.md5').read_text().splitlines() if line.strip()]
    with zipfile.ZipFile(dest/'bab2-analysis.zip') as archive:
        for checksum, member in members:
            archive.extract(member, dest/'x')
            if md5(dest/'x'/member) != checksum:
                raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified {len(FILES)} files and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
