"""Download the Materials Cloud record 1qt4h-tr585 (data.tar.xz, README.txt), verify checksums, unpack.

    python3 -B fetch.py [DEST]   (default ~/.cache/ai-theorist/sngeo2-coherent-miscibility)
"""
import hashlib, os, sys, tarfile, urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/1qt4h-tr585/files/'
FILES = {'data.tar.xz': ('ba7dcda1497afecf1863f2975a73eca0',
                         '98aef23f3b7adb9d7c707bd5cefe727d431f813946a4706be562c4b1365eaa26'),
         'README.txt': ('c44362bc1f34d2d6b172f9f2a982a24e',
                        '4afdd69ec91b4ec76d134291e817f25cf4ad5b2238b19d20b90510d56ef0adbe')}


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/.cache/ai-theorist/sngeo2-coherent-miscibility'))
    dest.mkdir(parents=True, exist_ok=True)
    for name, (md5, sha) in FILES.items():
        path = dest/name
        if not path.exists():
            urllib.request.urlretrieve(RECORD + name + '/content', path)
        blob = path.read_bytes()
        assert hashlib.md5(blob).hexdigest() == md5, name
        assert hashlib.sha256(blob).hexdigest() == sha, name
        print(name, len(blob), 'bytes, checksums ok')
    out = dest/'x'
    if not out.exists():
        with tarfile.open(dest/'data.tar.xz') as tar:
            tar.extractall(out)
    print('unpacked into', out)


if __name__ == '__main__':
    main()
