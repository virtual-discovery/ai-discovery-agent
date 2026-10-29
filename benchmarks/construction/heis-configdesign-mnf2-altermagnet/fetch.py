"""Refetch the pinned MnF2 archive from Materials Cloud record 3p6fv-8hd59 and unpack it without POTCARs.

    python3 benchmarks/construction/heis-configdesign-mnf2-altermagnet/fetch.py [CACHE_DIR]

CACHE_DIR defaults to ~/.cache/ai-theorist/heis-configdesign-mnf2-altermagnet. Only MnF2.tar.gz
(51,865,693 bytes, md5 4311cefc120956738f6d26440354f3fd, CC BY 4.0) is needed; the Fe archive is not.
VASP POTCAR files (licensed) are skipped during extraction and are never staged.
"""
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

URL = 'https://archive.materialscloud.org/api/records/3p6fv-8hd59/files/MnF2.tar.gz/content'
MD5 = '4311cefc120956738f6d26440354f3fd'
SIZE = 51865693
CACHE = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/heis-configdesign-mnf2-altermagnet')


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    archive = CACHE/'MnF2.tar.gz'
    if not archive.exists() or archive.stat().st_size != SIZE:
        urllib.request.urlretrieve(URL, archive)
    digest = hashlib.md5(archive.read_bytes()).hexdigest()
    assert digest == MD5, f'md5 mismatch: {digest}'
    with tarfile.open(archive) as tar:
        members = [m for m in tar.getmembers() if not m.name.endswith('/POTCAR') and '__pycache__' not in m.name]
        tar.extractall(CACHE, members=members)
    print('unpacked', len(members), 'members into', CACHE/'MnF2')


if __name__ == '__main__':
    main()
