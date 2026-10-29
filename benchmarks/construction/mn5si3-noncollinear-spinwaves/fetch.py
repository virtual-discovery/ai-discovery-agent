"""Refetch the pinned Materials Cloud files for mn5si3-noncollinear-spinwaves and unpack the members used.

    python3 -B fetch.py CACHE        # CACHE = ~/.cache/ai-theorist/mn5si3-noncollinear-spinwaves

Both records are CC BY 4.0. Downloads 126 MB; unpacks only the spin-wave/neutron folders (no DFT potentials).
"""
import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path

RECORDS = {
    'paxx3-a1k09': {'Mn5Si3.zip': 'bd4cfd1a5a489c22560ff0bf50f5432b', 'README.txt': '079bc7495b5a6fae9b89b2598aeede70'},
    '4tgt6-1p378': {'Mn5Si3.zip': '1c5a20826c98f71ab10322b50b9f2315', 'README.txt': '9ca1cb59f38d539421a7a8c8097e33e5'},
}
UNPACK = {'paxx3-a1k09': ('x22', ['Spin-wave calculations/*', 'DFT/inputcard']),
          '4tgt6-1p378': ('x24', ['Fig 1c/*', 'Fig 2/inputfiles/*', 'Fig 2/INS_data_Mn5Si3/*', 'Fig 2/spirit_output/spinconfig_MnSinoncol*final*.ovf',
                                  'Fig 2/outputfiles/dispersion_MnSiDisper*.dat', 'Fig 2/outputfiles/kpath_MnSiDisper*.dat',
                                  'Fig 2/inputcard_*.inp', 'Fig 2/input_*.cfg', 'Fig 3a/*', 'Fig 3c/*', 'Fig 4a/*', 'Fig 4b/*',
                                  'Fig 4c/*', 'Fig 5b/*', 'Fig 5c/*'])}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main(cache):
    cache = Path(cache)
    for rec, files in RECORDS.items():
        (cache/rec).mkdir(parents=True, exist_ok=True)
        for name, want in files.items():
            dest = cache/rec/name
            if not dest.exists() or md5(dest) != want:
                url = f'https://archive.materialscloud.org/api/records/{rec}/files/{name}/content'
                urllib.request.urlretrieve(url, dest)
            assert md5(dest) == want, (rec, name)
        out, members = UNPACK[rec]
        subprocess.run(['unzip', '-q', '-o', str(cache/rec/'Mn5Si3.zip'), *members, '-x', '__MACOSX/*', '-d', str(cache/out)],
                       check=False)
    print('ok')


if __name__ == '__main__':
    main(sys.argv[1])
