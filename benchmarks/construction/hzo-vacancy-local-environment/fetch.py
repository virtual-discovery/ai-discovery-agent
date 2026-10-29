"""Refetch the pinned source files for hzo-vacancy-local-environment and extract the members used.

    python3 benchmarks/construction/hzo-vacancy-local-environment/fetch.py [--dest DIR]

Source: S. Taper, B. Cai, N. Strkalj, B. Monserrat, G. Di Martino, Research Data supporting "Role of Local
Cation Environment on the Formation Energy of Oxygen Vacancies in Ferroelectric Hf1-xZrxO2", Apollo -
University of Cambridge Repository (2026), DOI 10.17863/CAM.121316, handle 1810/389433, CC BY 4.0.
Data.zip is 12.8 MB, so it is fetched whole (well under the 3 GB cap); only the members listed below are
extracted. Checksums: MD5 from the repository's bitstream records, SHA256 computed on download (2026-10-19).
"""
import argparse
import hashlib
import urllib.request
import zipfile
from pathlib import Path

API = 'https://api.repository.cam.ac.uk/server/api/core/bitstreams/%s/content'
FILES = {
    'Data.zip': ('62e6e6c4-bd26-4c0a-897e-b703211c9e48', 12791098, '70721b2c03a34c84e9a2deab61160afb',
                 'e6fd7f7aa4ed5ba2154fcf06dd9808d8c643637576b8b489887780942a0c4ea9'),
    'Readme.md': ('562efebd-5797-43d4-85eb-f176d5f76cfd', 15116, '77a9d035a0ac4ff76779ff70356b7934',
                  '89705f7d31c61d01e46aae37d378e4bf644e61ee74513a3cd95bde7b7c5304fe'),
}
PREFIX = 'Apollo deposit/'
MEMBERS = [
    'Main text/Structures/Fig2b_e_0.25_HZO.cif',
    'Main text/Structures/Fig2c_f_0.50_HZO.cif',
    'Main text/Structures/FIG2d_g_0.75_HZO.cif',
    'Main text/Fig_2/Data_Fig2b_g.csv',
    'SM material/FIGS3.csv',
    'SM material/FIGS4.txt',
    'SM material/FIGS4 3x3x3_supercell.cif',
    'SM material/FIGS5.txt',
    'SM material/FIGS12.csv',
]
DEFAULT = Path.home()/'.cache/ai-theorist/hzo-vacancy-local-environment'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dest', type=Path, default=DEFAULT)
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    for name, (uuid, size, md5, sha) in FILES.items():
        path = a.dest/name
        if not path.exists() or digest(path, 'md5') != md5:
            urllib.request.urlretrieve(API % uuid, path)
        assert path.stat().st_size == size, name
        assert digest(path, 'md5') == md5 and digest(path, 'sha256') == sha, f'checksum mismatch: {name}'
    with zipfile.ZipFile(a.dest/'Data.zip') as z:
        for member in MEMBERS:
            z.extract(PREFIX + member, a.dest/'raw')
    print('fetched and verified; extracted', len(MEMBERS), 'members under', a.dest/'raw'/PREFIX)


if __name__ == '__main__':
    main()
