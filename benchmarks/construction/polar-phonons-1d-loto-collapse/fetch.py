"""Refetch the pinned Materials Cloud archive 08gsg-gpd20 and verify every member used by the task.

    python3 -B fetch.py [DEST]     # default DEST: ~/.cache/ai-theorist/polar-phonons-1d-loto-collapse

Archive: N. Rivano, N. Marzari, T. Sohier, "Infrared-active phonons in one-dimensional materials and their
spectroscopic signatures", Materials Cloud Archive record 08gsg-gpd20 (DOI 10.24435/materialscloud:46-wj),
CC BY 4.0. The record holds two files, downloaded through the archive's per-file links (30.8 MB in total):

  materialscloud_data_Rivano.zip   30,776,903 bytes  md5 f3ee3048a131eb831eecee865d51e705
  files_description.md                    280 bytes  md5 0c399e8707d568e68028464c7cb683b0

The zip is unpacked into DEST/materialscloud_data_Rivano (134 MB; macOS metadata skipped). members.json records
the sha256 of every unpacked member; this script writes it on the first run and verifies it afterwards.
"""
import hashlib
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/polar-phonons-1d-loto-collapse'
API = 'https://archive.materialscloud.org/api/records/08gsg-gpd20/files/'
FILES = {'materialscloud_data_Rivano.zip': ('f3ee3048a131eb831eecee865d51e705', 30776903),
         'files_description.md': ('0c399e8707d568e68028464c7cb683b0', 280)}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(2 ** 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for name, (digest, size) in FILES.items():
        path = DEST/name
        if not (path.exists() and path.stat().st_size == size and md5(path) == digest):
            print('downloading', name)
            urllib.request.urlretrieve(API + name + '/content', path)
        if md5(path) != digest:
            raise SystemExit(f'md5 mismatch for {name}')
    record = HERE/'members.json'
    known = json.loads(record.read_text()) if record.exists() else {}
    out = {}
    with zipfile.ZipFile(DEST/'materialscloud_data_Rivano.zip') as z:
        for info in z.infolist():
            if info.is_dir() or info.filename.startswith('__MACOSX') or info.filename.endswith('.DS_Store'):
                continue
            data = z.read(info)
            path = DEST/info.filename
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists() or path.stat().st_size != len(data):
                path.write_bytes(data)
            digest = hashlib.sha256(data).hexdigest()
            if info.filename in known and known[info.filename]['sha256'] != digest:
                raise SystemExit(f'sha256 mismatch for {info.filename}')
            out[info.filename] = {'bytes': len(data), 'sha256': digest}
    if set(out) != set(known):
        record.write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print(f'{len(out)} members verified in {DEST}')


if __name__ == '__main__':
    main()
