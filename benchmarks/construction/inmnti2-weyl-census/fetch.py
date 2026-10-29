"""Download and verify the InMnTi2 files of Materials Cloud Archive record hkprp-xw234.

    python3 -B benchmarks/construction/inmnti2-weyl-census/fetch.py [DEST]

Source: D. Grassano, L. Binci, N. Marzari, "Type-I antiferromagnetic Weyl semimetal InMnTi2",
Phys. Rev. Research 6, 013140 (2024), arXiv:2208.11412; archive DOI 10.24435/materialscloud:ph-3c (CC BY 4.0).

Downloads upload.zip (160,482,239 bytes) and files_description.md through the archive's per-file links and
checks their md5. Extracts the 96 members listed in members.md5 (QE inputs/outputs, Wannier90 .win and hr.dat,
the Kubo outputs, WannierTools node/chirality outputs, README) and checks each md5. The two Wannier90 .wout
files (579 MB and 49 MB) are not extracted; their header (first 400 lines) and final-state block are written
as <member>.excerpt by streaming. Surface spectral functions and Fermi-arc maps are not used.
DEST defaults to ~/.cache/ai-theorist/inmnti2-weyl-census; members land in DEST/x/.
"""
import collections
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/hkprp-xw234/files'
FILES = {'upload.zip': '27d1aa76a5cc384765108beb2ae05289',
         'files_description.md': '37d5be45010d2ae698b72c26c1dcd574'}
WOUT = ['hubbard/wannier/InMnTi2.wout', 'standard/wannier/InMnTi2.wout']
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def excerpt(archive, member, target):
    head, tail = [], collections.deque(maxlen=120)
    with archive.open(member) as f:
        for i, line in enumerate(f):
            line = line.decode(errors='replace')
            if i < 400:
                head.append(line)
            tail.append(line)
    tail = list(tail)
    start = max(i for i, line in enumerate(tail) if 'Final State' in line) if any('Final State' in l for l in tail) else 0
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(''.join(head) + '\n[... lines omitted ...]\n\n' + ''.join(tail[start:]))


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/inmnti2-weyl-census')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = [line.split(None, 1) for line in (HERE/'members.md5').read_text().splitlines() if line.strip()]
    with zipfile.ZipFile(dest/'upload.zip') as archive:
        for checksum, member in members:
            archive.extract(member, dest/'x')
            if md5(dest/'x'/member) != checksum:
                raise SystemExit(f'{member}: md5 mismatch')
        for member in WOUT:
            excerpt(archive, member, dest/'x'/(member + '.excerpt'))
    print(f'verified {len(FILES)} files and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
