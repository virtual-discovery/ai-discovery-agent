"""Refetch the pinned source files for hydrogen-plasma-facing-materials.

    python3 benchmarks/construction/hydrogen-plasma-facing-materials/fetch.py [--dest DIR]

Source: A. Fedrigucci, N. Marzari, P. Ricci, "DFT calculations of surface binding and interstitial hydrogen
formation energies for plasma-facing materials", Materials Cloud Archive record gpxwh-dx645,
DOI 10.24435/materialscloud:9k-09, CC BY 4.0 (paper: PRX Energy 3, 043002 (2024), arXiv:2407.00858).

The tarball is 46.9 MB (765 MB unpacked, 490 Quantum ESPRESSO runs), well under the 3 GB cap, so it is fetched
whole and unpacked under DEST/raw. Checksums are the archive's own md5 values (record API, 2026-10-19).
"""
import argparse
import hashlib
import tarfile
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/gpxwh-dx645/files/%s/content'
FILES = {
    'PFMs_screening_QE_data.tar.xz': (46921676, '56aadb3e22b5a0d52ff75d5757798aaa'),
    'README.md': (2514, '61db996d231cb32b2fbf19328de17307'),
    'files_description.md': (276, '21acb466968792b84ff2903d878da8c5'),
}
DEFAULT = Path.home()/'.cache/ai-theorist/hydrogen-plasma-facing-materials'


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dest', type=Path, default=DEFAULT)
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    for name, (size, digest) in FILES.items():
        path = a.dest/name
        if not path.exists() or md5(path) != digest:
            urllib.request.urlretrieve(RECORD % name, path)
        assert path.stat().st_size == size and md5(path) == digest, f'checksum mismatch: {name}'
    raw = a.dest/'raw'
    if not (raw/'H-IFE').exists():
        with tarfile.open(a.dest/'PFMs_screening_QE_data.tar.xz') as t:
            t.extractall(raw)
    n = len(list(raw.rglob('pwscf.out')))
    assert n == 490, n
    print('fetched and verified;', n, 'QE runs under', raw)


if __name__ == '__main__':
    main()
