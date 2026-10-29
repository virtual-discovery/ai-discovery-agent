"""Fetch the pinned Zenodo archive of Kundu et al., PRL 132, 116301 (2024) and verify checksums.

    python3 -B fetch.py [DEST]     (default ~/.cache/ai-theorist/tan-pressure-kappa-channels)

Writes DEST/IFCs-data-pressure-TaN.tar.gz (25.5 MB) and extracts it to DEST/raw/.
"""
import hashlib
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

URL = 'https://zenodo.org/api/records/10564835/files/IFCs-data-pressure-TaN.tar.gz/content'
MD5 = '0e91973a88688ca03ac776d379486f9c'
HERE = Path(__file__).parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/tan-pressure-kappa-channels'
    dest.mkdir(parents=True, exist_ok=True)
    tgz = dest/'IFCs-data-pressure-TaN.tar.gz'
    if not tgz.exists() or md5(tgz) != MD5:
        urllib.request.urlretrieve(URL, tgz)
    assert md5(tgz) == MD5, 'archive checksum mismatch'
    with tarfile.open(tgz) as t:
        t.extractall(dest/'raw')
    members = json.loads((HERE/'members.json').read_text())
    for name, digest in members.items():
        assert md5(dest/'raw'/name) == digest, name
    print('ok', len(members), 'members verified in', dest/'raw')


if __name__ == '__main__':
    main()
