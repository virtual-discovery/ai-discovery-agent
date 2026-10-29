"""Refetch the pinned archive files (Materials Cloud record bhv3r-a8q53, CC BY 4.0) and verify checksums.

    python3 fetch.py     # downloads into ~/.cache/ai-theorist/mpea-vacancy-diffusion-kmc/
"""
import hashlib
import urllib.request
from pathlib import Path

DEST = Path.home()/'.cache/ai-theorist/mpea-vacancy-diffusion-kmc'
BASE = 'https://archive.materialscloud.org/api/records/bhv3r-a8q53/files/'
FILES = {'all_data.json': '15df9136ed668320057131c8c59cad2a',      # 143,013,225 bytes
         'files_description.md': '8849580fa4e1bb9f72bca8d1e74e571f'}


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
            urllib.request.urlretrieve(BASE + name + '/content', path)
        assert md5(path) == digest, f'checksum mismatch for {name}'
        print(name, 'ok', path.stat().st_size)


if __name__ == '__main__':
    main()
