"""Download and verify the files of Materials Cloud Archive record rfcdy-wcy32 used by 2d-qsh-shc-plateau.

    python3 -B benchmarks/construction/2d-qsh-shc-plateau/fetch.py [DEST]

Source: J. Zhou, S. Ponce, J.-C. Charlier, "High-throughput calculations of spin Hall conductivity in non-magnetic 2D
materials", npj 2D Mater. Appl. 9, 39 (2025), arXiv:2501.02110; archive DOI 10.24435/materialscloud:yv-h7 (CC BY 4.0).

Downloads Data.7z (408,780,659 bytes), README.txt and files_description.md through the archive's per-file links and
checks their md5. Data.7z is one solid LZMA2 block, so no range read of single members is possible; the whole file is
downloaded (409 MB, within the 3 GB cap) and only the members listed in members.md5 are extracted (about 110 MB) with
p7zip (`7z`, host tool) and checked against their md5. DEST defaults to ~/.cache/ai-theorist/2d-qsh-shc-plateau;
members land in DEST/x/.
"""
import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/rfcdy-wcy32/files'
FILES = {'Data.7z': '1d261e780f092b8c7c1f2e051855c798',
         'README.txt': 'aedc173ccb0587b13fad1cbbcd5f58cc',
         'files_description.md': '10c58cb5ed5d38de19f445f2fca21c19'}
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/2d-qsh-shc-plateau')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = [line.split(None, 1) for line in (HERE/'members.md5').read_text().splitlines() if line.strip()]
    todo = [m for c, m in members if not (dest/'x'/m).exists() or md5(dest/'x'/m) != c]
    if todo:
        (dest/'members.lst').write_text('\n'.join(todo) + '\n')
        subprocess.run(['7z', 'x', '-y', f'-o{dest/"x"}', str(dest/'Data.7z'), f'@{dest/"members.lst"}'],
                       check=True, stdout=subprocess.DEVNULL)
    for checksum, member in members:
        if md5(dest/'x'/member) != checksum:
            raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified {len(FILES)} files and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
