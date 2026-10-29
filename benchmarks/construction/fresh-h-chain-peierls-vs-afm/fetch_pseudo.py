"""Refetch the pinned pseudopotential for fresh-h-chain-peierls-vs-afm and vendor it into the lab.

    python3 benchmarks/construction/fresh-h-chain-peierls-vs-afm/fetch_pseudo.py

Source: SSSP 1.3.0 PBE precision, Materials Cloud Archive record rcyfm-68h65
(DOI 10.24435/materialscloud:f3-ym; Prandini et al., npj Comput. Mater. 4, 72 (2018)).
Downloads the library tarball, its metadata JSON and LICENSE.txt into
~/.cache/ai-theorist/fresh-h-chain-peierls-vs-afm/, verifies the archive MD5s published by
the record API, extracts H_ONCV_PBE-1.0.oncvpsp.upf, checks it against the MD5 in the SSSP
metadata and the SHA256 pinned below, and copies it to environment/lab/ (lab image only; the
agent image never sees it). The SG15 family is licensed CC BY-SA 4.0 (see LICENSE.txt).
"""
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/rcyfm-68h65/files/{}/content'
FILES = ['SSSP_1.3.0_PBE_precision.tar.gz', 'SSSP_1.3.0_PBE_precision.json', 'LICENSE.txt']   # MD5s from the record API
UPF = 'H_ONCV_PBE-1.0.oncvpsp.upf'
UPF_MD5 = '1790becc920ee074925cf490c71280fe'          # from SSSP_1.3.0_PBE_precision.json
UPF_SHA256 = '69d631780c0178c5239ed3dcb4a3e22c13bd13cde9beb97e123d036949da5c09'
CACHE = Path.home()/'.cache/ai-theorist/fresh-h-chain-peierls-vs-afm'
HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]/'tasks/fresh-h-chain-peierls-vs-afm/environment/lab'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def api_md5():
    with urllib.request.urlopen('https://archive.materialscloud.org/api/records/rcyfm-68h65', timeout=60) as r:
        entries = json.load(r)['files']['entries']
    return {name: entries[name]['checksum'].split(':', 1)[1] for name in FILES}


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    md5 = api_md5()
    for name in FILES:
        target = CACHE/name
        if not target.exists() or digest(target, 'md5') != md5[name]:
            print('downloading', name)
            urllib.request.urlretrieve(RECORD.format(name), target)
        assert digest(target, 'md5') == md5[name], f'MD5 mismatch for {name}'
        print(name, target.stat().st_size, 'md5', md5[name])
    meta = json.loads((CACHE/'SSSP_1.3.0_PBE_precision.json').read_text())
    assert meta['H']['filename'] == UPF and meta['H']['md5'] == UPF_MD5 and meta['H']['pseudopotential'] == 'SG15'
    with tarfile.open(CACHE/'SSSP_1.3.0_PBE_precision.tar.gz') as tar:
        member = next(m for m in tar.getmembers() if m.name.endswith('/' + UPF) or m.name == UPF)
        with tar.extractfile(member) as src, open(CACHE/UPF, 'wb') as dst:
            shutil.copyfileobj(src, dst)
    assert digest(CACHE/UPF, 'md5') == UPF_MD5, 'UPF MD5 differs from the SSSP metadata'
    sha = digest(CACHE/UPF, 'sha256')
    assert UPF_SHA256 is None or sha == UPF_SHA256, 'UPF SHA256 differs from the pinned value'
    shutil.copy2(CACHE/UPF, LAB/UPF)
    print(UPF, 'md5', UPF_MD5, 'sha256', sha, '->', LAB/UPF)
    print('SSSP recommended cutoffs for H (Ry):', meta['H'].get('cutoff_wfc'), meta['H'].get('cutoff_rho'))


if __name__ == '__main__':
    main()
