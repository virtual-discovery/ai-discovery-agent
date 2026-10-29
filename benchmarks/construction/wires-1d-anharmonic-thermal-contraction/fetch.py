"""Refetch the pinned Materials Cloud archive k45dq-ebf35 (DOI 10.24435/materialscloud:fj-20, CC BY 4.0).

    python3 -B fetch.py [DEST]   # default DEST: ~/.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction

Downloads README.txt and the 145.7 MB zip, checks the md5 published by the archive API, and unpacks only the
SSCHA/ and QHA/ folders (the ML/ potentials are not used) into DEST/raw/.
"""
import hashlib
import sys
import zipfile
from pathlib import Path
import requests

BASE = 'https://archive.materialscloud.org/api/records/k45dq-ebf35/files/'
FILES = {'README.txt': 'c866d88afdf77c93202d26d1c30e0b38',
         'Extreme-anharmonicity-thermal-contration-1D-wires.zip': '1de3fb44caf32b1ed403cbe0e7b687ce'}
DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction'


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for name, digest in FILES.items():
        path = DEST/name
        if not path.exists() or md5(path) != digest:
            with requests.get(BASE + name + '/content', stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(path, 'wb') as f:
                    for block in r.iter_content(1 << 20):
                        f.write(block)
        assert md5(path) == digest, f'md5 mismatch for {name}'
        print(name, digest, 'ok')
    with zipfile.ZipFile(DEST/'Extreme-anharmonicity-thermal-contration-1D-wires.zip') as z:
        members = [m for m in z.namelist() if m.split('/')[1:2] in (['SSCHA'], ['QHA']) and '.DS_Store' not in m]
        z.extractall(DEST/'raw', members)
    print(len(members), 'members unpacked into', DEST/'raw')


if __name__ == '__main__':
    main()
