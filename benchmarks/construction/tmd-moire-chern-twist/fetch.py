"""Download and verify the Zenodo record used by tmd-moire-chern-twist, and extract the members the build reads.

    python3 -B benchmarks/construction/tmd-moire-chern-twist/fetch.py [DEST]

Source: X.-W. Zhang, C. Wang, X. Liu, Y. Fan, T. Cao, D. Xiao, "Polarization-driven band topology evolution in twisted
MoTe2 and WSe2", Nat. Commun. 15, 4223 (2024), arXiv:2311.12776; data: Zenodo record 10888769
(DOI 10.5281/zenodo.10888769, CC BY 4.0). One file, source-data.tar.gz (20,803,943 bytes, md5 038b59d5...).
The ML/MD-VASP/POTCAR member (a licensed VASP potential) is never extracted. DEST defaults to
~/.cache/ai-theorist/tmd-moire-chern-twist; members land in DEST/x/source-data/.
"""
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

URL = 'https://zenodo.org/api/records/10888769/files/source-data.tar.gz/content'
MD5 = '038b59d589a8f254d76599e6b3e96f74'
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/tmd-moire-chern-twist')
    dest.mkdir(parents=True, exist_ok=True)
    tgz = dest/'source-data.tar.gz'
    if not tgz.exists() or md5(tgz) != MD5:
        urllib.request.urlretrieve(URL, tgz)
    if md5(tgz) != MD5:
        raise SystemExit('source-data.tar.gz: md5 mismatch')
    members = [line.split(None, 1) for line in (HERE/'members.md5').read_text().splitlines() if line.strip()]
    wanted = {m for _, m in members}
    with tarfile.open(tgz) as tar:
        for info in tar.getmembers():
            if info.name in wanted:
                tar.extract(info, dest/'x')
    for checksum, member in members:
        if md5(dest/'x'/member) != checksum:
            raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified source-data.tar.gz and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
