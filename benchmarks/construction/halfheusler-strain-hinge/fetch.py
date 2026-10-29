"""Download and verify the archive of Zenodo record 16927281 (CC BY 4.0).

    python3 -B benchmarks/construction/halfheusler-strain-hinge/fetch.py [DEST]

Source: S. K. Das, I. C. Fulga, R. Dhawan, H. C. Kandpal, J. van den Brink, J. I. Facio, "Electronic structure and
hinge states of strained half-Heusler compounds LiSbZn and LiBiZn", Phys. Rev. B 112, 125132 (2025), arXiv:2504.20190;
data DOI 10.5281/zenodo.16927281 ("Hinge electronic structure of strained half-Heuslers").

Downloads strained_hH_zenodo.zip (2,681,691 bytes) through the record's per-file link, checks its md5, extracts all
members into DEST/x/ and checks each member's md5 against MEMBERS below.  DEST defaults to
~/.cache/ai-theorist/halfheusler-strain-hinge.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

URL = 'https://zenodo.org/api/records/16927281/files/strained_hH_zenodo.zip/content'
ZIP_MD5 = 'f66b568ce17e4055122facfe02e01e24'
ZIP_SHA256 = '0d4ac9e4cb785c7f023050f4660655535f4b22e4963685475879b797f2a2c5bb'
MEMBERS = {
    'code/strained_hH.py': '1cf64e62695bc7f7b488fd75af947331',
    'wannier/+hamdata_libizn_eq_orig_setup.txt': 'c387330316fa084fa2383ddec25b87d8',
    'wannier/+hamdata_libizn_eq_strained_setup.txt': 'c341c104858f00961b8f13af8949fe95',
    'wannier/+hamdata_lisbzn_eq_orig_setup.txt': 'e46182ddfa7c6319d9e135a84f65a59c',
    'wannier/+hamdata_lisbzn_eq_strained_setup.txt': '1657fad07d1e43e9cc3d6952aa20dc08',
    'kwant_ham/LiBiZn_nostrain.ham': '0b5608916402527457ef4c0036854957',
    'kwant_ham/LiBiZn_strain.ham': '8aece76502dcd6b32f628b89dc04c630',
    'kwant_ham/LiSbZn_nostrain.ham': 'e50458c596cfef7205243c973802a4b2',
    'kwant_ham/LiSbZn_strain.ham': '65eac7f5885161ae0d1e58e4f106c273',
}


def digest(path, algo='md5'):
    h = hashlib.new(algo)
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/halfheusler-strain-hinge')
    dest.mkdir(parents=True, exist_ok=True)
    target = dest/'strained_hH_zenodo.zip'
    if not target.exists() or digest(target) != ZIP_MD5:
        urllib.request.urlretrieve(URL, target)
    if digest(target) != ZIP_MD5 or digest(target, 'sha256') != ZIP_SHA256:
        raise SystemExit('strained_hH_zenodo.zip: checksum mismatch')
    with zipfile.ZipFile(target) as archive:
        for member, checksum in MEMBERS.items():
            archive.extract(member, dest/'x')
            if digest(dest/'x'/member) != checksum:
                raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified strained_hH_zenodo.zip and {len(MEMBERS)} members in {dest}')


if __name__ == '__main__':
    main()
