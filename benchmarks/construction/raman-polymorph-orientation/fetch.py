"""Refetch the pinned Computational Raman Database release and extract the three SnF2 records.

    python3 -B fetch.py [DEST]

Downloads Computational_Raman_Database_v1.0.1.zip (152,753,693 bytes, md5 e2aaea03df70fda95dbd970c60456882)
from Materials Cloud Archive record 0m9r0-g3486 (DOI 10.24435/materialscloud:pg-h3, CC BY 4.0) into DEST
(default ~/.cache/ai-theorist/raman-polymorph-orientation), checks the md5, reads the member crd.json
(454,304,470 bytes, an ASE-database JSON document with 5099 records) and writes pool.json next to this script:
the archived records 3216 (mp-2843, C2/c), 2564 (mp-7456, P2_12_12_1) and 2390 (mp-7457, P4_12_12), unchanged.
The sha256 of each record's canonical JSON (sorted keys) is asserted.
"""
import hashlib
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

URL = ('https://archive.materialscloud.org/api/records/0m9r0-g3486/files/'
       'Computational_Raman_Database_v1.0.1.zip/content')
MD5 = 'e2aaea03df70fda95dbd970c60456882'
RECORDS = {'3216': 'mp-2843', '2564': 'mp-7456', '2390': 'mp-7457'}
SHA256 = {'3216': 'da176f4d3a76e91f399bd8ca8a1404f0d92c238d31555169d5b961917c860d4c', '2564': '849d4bfea4de74b2e6022c51bfbd4bf999ee1e1bccea4ccd59f3465c2999c376', '2390': 'e5735d7dc160445293f7135b0afd19a52d10710ce85306e1a26d42a6b06f7193'}
HERE = Path(__file__).resolve().parent


def canonical_sha(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()


def main():
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/raman-polymorph-orientation'
    dest.mkdir(parents=True, exist_ok=True)
    zpath = dest/'Computational_Raman_Database_v1.0.1.zip'
    if not zpath.exists():
        urllib.request.urlretrieve(URL, zpath)
    md5 = hashlib.md5()
    with open(zpath, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            md5.update(chunk)
    assert md5.hexdigest() == MD5, 'md5 mismatch for the database zip'
    with zipfile.ZipFile(zpath) as z:
        db = json.loads(z.read('crd.json'))
    pool = {}
    for key, mpid in RECORDS.items():
        rec = db[key]
        assert rec['key_value_pairs']['mpid'] == mpid
        if not SHA256[key].startswith('{'):
            assert canonical_sha(rec) == SHA256[key], f'record {key} changed'
        pool[key] = rec
        print(key, mpid, rec['key_value_pairs']['spacegroup'], canonical_sha(rec))
    (HERE/'pool.json').write_text(json.dumps(pool, sort_keys=True) + '\n')
    print('wrote', HERE/'pool.json')


if __name__ == '__main__':
    main()
