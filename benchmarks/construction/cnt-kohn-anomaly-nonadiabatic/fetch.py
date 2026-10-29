"""Refetch the pinned Materials Cloud archive files for cnt-kohn-anomaly-nonadiabatic and verify them.

    python3 -B fetch.py [DEST]      (default ~/.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic)

Downloads the four record files (27 MB in total) through the archive's per-file API links, checks
their md5 against the record metadata, unpacks them into DEST/raw and checks every member listed in
members.md5. The prebuilt binary Fig1/coupling_calculator is not used (it is rebuilt from source).
"""
import hashlib
import sys
import zipfile
from pathlib import Path
import urllib.request

RECORD = 'https://archive.materialscloud.org/api/records/vs5fd-dbv98/files/{}/content'
FILES = {'README.md': '1c6dbfd0eb10564470d215e3eca226fe', 'dft_inputs.zip': 'c85fa7a8af0ffa84597e7fe76e1749b5',
         'Fig1.zip': '4c1760d1aaf80add6900cd65233ca571', 'Fig2.zip': 'b17a42c46d5a7856356d71eab691875c'}
HERE = Path(__file__).resolve().parent


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def main():
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic'
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        path = dest/name
        if not path.exists() or md5(path) != checksum:
            urllib.request.urlretrieve(RECORD.format(name), path)
        assert md5(path) == checksum, f'md5 mismatch for {name}'
        if name.endswith('.zip'):
            with zipfile.ZipFile(path) as z:
                z.extractall(dest/'raw')
    bad = []
    for line in (HERE/'members.md5').read_text().splitlines():
        checksum, member = line.split(maxsplit=1)
        if md5(dest/'raw'/member) != checksum:
            bad.append(member)
    assert not bad, f'member checksum mismatch: {bad}'
    print(f'verified {len(FILES)} files and {len((HERE/"members.md5").read_text().splitlines())} members in {dest}')


if __name__ == '__main__':
    main()
