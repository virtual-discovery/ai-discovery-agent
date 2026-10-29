"""Refetch the pinned source files for gan-irradiation-defect-levels and unpack them.

    python3 benchmarks/construction/gan-irradiation-defect-levels/fetch.py [--dest DIR] [--pseudo]

Source: A. S. Hauck, M. Jin, B. R. Tuttle, "Structures and energies of GaN native point defects, divacancy, and
carbon substitutional from hybrid functional calculations", Zenodo record 21893631 (2026),
DOI 10.5281/zenodo.21893631, CC BY 4.0. gan_defects.zip is 25.1 MB, so it is fetched whole through Zenodo's
per-file link (well under the 3 GB cap) and unpacked into DEST/raw (115 MB, 67 calculation folders).

--pseudo also fetches the two pslibrary 1.0.0 pseudopotentials used for the construction-time DFPT dielectric
tensors (dfpt/run_dfpt.sh) into DEST/dfpt.
"""
import argparse
import hashlib
import urllib.request
import zipfile
from pathlib import Path

ZENODO = 'https://zenodo.org/api/records/21893631/files/%s/content'
FILES = {'gan_defects.zip': (25125074, '85d4a96b1bf0bcf032c89bf92dbb021f',
                             '5112c84e0d6cd8f5041c67da3fb1ece63a087dcd1d191641e5de1664e929fe25')}
PSLIB = 'https://pseudopotentials.quantum-espresso.org/upf_files/%s'
PSEUDO = {'Ga.pbe-dn-rrkjus_psl.1.0.0.UPF': (1710832, '4451a9dabd3471075670e2a6b0dab093c036eadc0d5ee43f05f4c688bf67ff44'),
          'N.pbe-n-rrkjus_psl.1.0.0.UPF': (773856, '311bcb45ca5ad404fbab64e88d1bbe2e762c4582d841747fb9fbf4eb7f533e4a')}
DEFAULT = Path.home()/'.cache/ai-theorist/gan-irradiation-defect-levels'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dest', type=Path, default=DEFAULT)
    p.add_argument('--pseudo', action='store_true')
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    for name, (size, md5, sha) in FILES.items():
        path = a.dest/name
        if not path.exists() or digest(path, 'md5') != md5:
            urllib.request.urlretrieve(ZENODO % name, path)
        assert path.stat().st_size == size, name
        assert digest(path, 'md5') == md5 and digest(path, 'sha256') == sha, f'checksum mismatch: {name}'
    with zipfile.ZipFile(a.dest/'gan_defects.zip') as z:
        z.extractall(a.dest/'raw')
    print('fetched and verified gan_defects.zip; unpacked under', a.dest/'raw')
    if a.pseudo:
        (a.dest/'dfpt').mkdir(exist_ok=True)
        for name, (size, sha) in PSEUDO.items():
            path = a.dest/'dfpt'/name
            if not path.exists():
                urllib.request.urlretrieve(PSLIB % name, path)
            assert path.stat().st_size == size and digest(path, 'sha256') == sha, name
        print('fetched pslibrary pseudopotentials into', a.dest/'dfpt')


if __name__ == '__main__':
    main()
