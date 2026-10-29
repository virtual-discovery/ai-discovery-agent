"""Download and verify the files used from Materials Cloud Archive record zxf0b-rs406.

    python3 -B benchmarks/construction/metal-resistivity-te-vs-anharmonicity/fetch.py [DEST]

Source: A. Wang, J. Yin, F. A. Goudreault, M. Cote, O. Hellman, S. Ponce, "Opposite impact of thermal expansion and
phonon anharmonicity on the phonon-limited resistivity of elemental metals from first principles" (2026),
DOI 10.1103/yk57-bt6t; archive DOI 10.24435/materialscloud:4h-r1 (CC BY 4.0).

Downloads README.md, Acknowledgements.md and data_for_paper.tar.gz (237,561,636 bytes) through the archive's
per-file links (a tar.gz has no index, so the whole file is streamed once), checks their md5, then extracts only
the 150 members listed in members.md5 (EPW ibte inputs/outputs and scf inputs, q2r and TDEP force constants,
TDEP unit cells, READMEs) and checks each member's md5. The three 0.2-0.3 GB spectral-decomposition files are
never written to disk. DEST defaults to ~/.cache/ai-theorist/metal-resistivity-te-vs-anharmonicity; members
land in DEST/x/.
"""
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/zxf0b-rs406/files'
FILES = {'README.md': '089442c4cc8bfd798f062c178048897c',
         'Acknowledgements.md': '968797ce577c1c08e4ab2b2195796e48',
         'data_for_paper.tar.gz': '83486659df1249b85d0025c2231d06b6'}
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else
                Path.home()/'.cache/ai-theorist/metal-resistivity-te-vs-anharmonicity')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)   # follows the 302 to the object store
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    wanted = dict(reversed(line.split(None, 1)) for line in (HERE/'members.md5').read_text().splitlines() if line.strip())
    with tarfile.open(dest/'data_for_paper.tar.gz', 'r:gz') as archive:
        for member in archive:
            if member.name in wanted:
                archive.extract(member, dest/'x')
    for member, checksum in wanted.items():
        if md5(dest/'x'/member) != checksum:
            raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified {len(FILES)} files and {len(wanted)} members in {dest}')


if __name__ == '__main__':
    main()
