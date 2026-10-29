"""Refetch the pinned Materials Cloud files for co3sn2s2-dmi-diffuse and unpack the scan arrays.

    python3 -B fetch.py [CACHE_DIR]

Default CACHE_DIR: ~/.cache/ai-theorist/co3sn2s2-dmi-diffuse. Downloads the six record files through the
archive API's per-file links into CACHE_DIR/raw, checks each md5, and extracts the .npy members (no __MACOSX
entries) of narrowscans.zip and widescans.zip into CACHE_DIR/x/{narrowscans,widescans}/.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/0cwzy-zz096/files/{}/content'
MD5 = {'narrowscans.zip': '1900407dd1f7b413ba3786980eb9103a',
       'widescans.zip': '5c96956ea4dace5ef9ca4ca488808de7',
       'README.md': 'd56203a3e02a037e8ee037b8fd4041d0',
       'readcombined.py': '237aa52e042fe890e98614bc636bc8a6',
       'SpinW_Co3Sn2S2_DMI.m': '83a27de47d4381f4933e04ab156bd6ea',
       'Co3Sn2S2-DMS-analysis.nb': 'b0e5102993736567538f35cd80a165e4'}
SIZES = {'narrowscans.zip': 1248776, 'widescans.zip': 10213141, 'README.md': 2064, 'readcombined.py': 812,
         'SpinW_Co3Sn2S2_DMI.m': 8578, 'Co3Sn2S2-DMS-analysis.nb': 7376932}


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


def main():
    cache = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/co3sn2s2-dmi-diffuse'
    raw = cache/'raw'
    raw.mkdir(parents=True, exist_ok=True)
    for name, digest in MD5.items():
        dest = raw/name
        if not dest.exists() or md5(dest) != digest:
            with urllib.request.urlopen(RECORD.format(name), timeout=600) as response:   # follows the 302 to storage
                dest.write_bytes(response.read())
        assert dest.stat().st_size == SIZES[name] and md5(dest) == digest, f'checksum mismatch for {name}'
        print(f'{name}: md5 ok')
    for archive in ('narrowscans.zip', 'widescans.zip'):
        with zipfile.ZipFile(raw/archive) as z:
            for member in z.namelist():
                if member.startswith('__MACOSX') or not member.endswith('.npy'):
                    continue
                target = cache/'x'/member
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(z.read(member))
    print('extracted to', cache/'x')


if __name__ == '__main__':
    main()
