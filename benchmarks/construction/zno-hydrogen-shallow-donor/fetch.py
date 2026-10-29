"""Stream ZnO_ubiquitous.tar.gz once from the Materials Cloud Archive, verify its md5, and extract
only members below a size limit (vasprun.xml, CHGCAR, WAVECAR... are skipped).

    python3 -B fetch.py [DEST]   (default ~/.cache/ai-theorist/zno-hydrogen-shallow-donor/raw)
Writes DEST/listing.json with every member name and size.
"""
import hashlib, json, os, sys, tarfile, urllib.request
from pathlib import Path

URL = 'https://archive.materialscloud.org/api/records/c91vp-p9h03/files/ZnO_ubiquitous.tar.gz/content'
MD5 = 'd21e5cfb405c798dcd23ad41c9678614'
README_URL = 'https://archive.materialscloud.org/api/records/c91vp-p9h03/files/README.txt/content'
LIMIT = 12 * 2**20
SKIP = ('vasprun.xml', 'WAVECAR', 'CHGCAR', 'CHG', 'LOCPOT', 'PROCAR', 'ELFCAR', 'AECCAR')


class Hashing:
    def __init__(self, f):
        self.f, self.h, self.n = f, hashlib.md5(), 0
    def read(self, n=-1):
        b = self.f.read(n)
        self.h.update(b); self.n += len(b)
        return b


def open_with_retry(url, tries=240, wait=60):
    import time, urllib.error
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=600)
        except urllib.error.HTTPError as exc:
            if exc.code not in (500, 502, 503, 504) or i == tries - 1:
                raise
            print('HTTP', exc.code, 'retrying', time.strftime('%H:%M:%S'), flush=True)
            time.sleep(wait)


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/.cache/ai-theorist/zno-hydrogen-shallow-donor/raw'))
    dest.mkdir(parents=True, exist_ok=True)
    stream = Hashing(open_with_retry(URL))
    listing = []
    with tarfile.open(fileobj=stream, mode='r|gz') as tar:
        for m in tar:
            listing.append({'name': m.name, 'size': m.size, 'type': 'dir' if m.isdir() else 'file'})
            base = os.path.basename(m.name)
            if m.isfile() and m.size <= LIMIT and not base.startswith(SKIP) and '..' not in m.name:
                tar.extract(m, dest)
    while stream.read(1 << 20):
        pass
    md5 = stream.h.hexdigest()
    (dest/'listing.json').write_text(json.dumps({'md5': md5, 'bytes': stream.n, 'members': listing}, indent=0))
    print('md5', md5, 'expected', MD5, 'bytes', stream.n)
    assert md5 == MD5, 'checksum mismatch'


if __name__ == '__main__':
    main()
