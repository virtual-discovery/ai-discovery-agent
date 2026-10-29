"""Download and verify the Zenodo record 23319151 (Sugawara 2026, Mg2Si/Mg2Ge/Mg2Sn QSGW archive).

    python3 -B benchmarks/construction/mg2x-midir-forbidden-transition/fetch.py [DEST]

Downloads the archive zip (77,333,006 bytes) and its README/LICENSES/MD5SUMS through Zenodo's per-file
links, checks their md5 against the record, extracts the members listed in members.md5 into DEST/x and checks
each member against members.md5 (whose checksums are those of the archive's own MD5SUMS.txt).
DEST defaults to ~/.cache/ai-theorist/mg2x-midir-forbidden-transition (read by model.py).
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://zenodo.org/api/records/23319151/files'
ZIP = '20261012aa-paper1-zenodo-archive.zip'
FILES = {ZIP: '865ba37d739a8118a8a556cb95c643bb', 'README.md': '246ab9331d014b4665cb1629ca218586',
         'LICENSES.md': 'aed582f0574ab7043a3627ea9ac59b82', 'MD5SUMS.txt': '91946a35a975e8bfcda1d31e4a6db173'}
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/mg2x-midir-forbidden-transition')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = [line.split(None, 1) for line in (HERE/'members.md5').read_text().splitlines() if line.strip()]
    archive_sums = dict(reversed(line.split(None, 1)) for line in (dest/'MD5SUMS.txt').read_text().splitlines() if line.strip())
    with zipfile.ZipFile(dest/ZIP) as z:
        names = {n.lstrip('./'): n for n in z.namelist()}
        for checksum, member in members:
            if archive_sums.get(member) != checksum:
                raise SystemExit(f'{member}: not as listed in the archive MD5SUMS.txt')
            z.extract(names[member], dest/'x')
            if md5(dest/'x'/names[member]) != checksum:
                raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified {len(FILES)} files and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
