"""Refetch the two pinned Ge pseudopotentials for fresh-ge-internal-strain and vendor them into the lab.

    python3 benchmarks/construction/fresh-ge-internal-strain/fetch_pseudos.py [--install]

1. ge_pbe_v1.4.uspp.F.UPF (GBRV 1.4, 3d in valence; the SSSP 1.3.0 PBE precision choice for Ge) from the SSSP
   library, Materials Cloud Archive record rcyfm-68h65 (DOI 10.24435/materialscloud:f3-ym; Prandini et al.,
   npj Comput. Mater. 4, 72 (2018)). The archive files are checked against the MD5s published by the record API,
   the UPF against the MD5 in the SSSP metadata and the SHA256 pinned below. License of the GBRV family: GPL v3
   (stated in the archive's LICENSE.txt); the rest of the record is CC BY 4.0.
2. Ge.pbe-n-rrkjus_psl.1.0.0.UPF (PSlibrary 1.0.0, A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014); 3d frozen in the
   core) from the Quantum ESPRESSO pseudopotential library. License: GNU GPL v2 or later. Checked against the
   SHA256 and MD5 pinned below.

Downloads go to ~/.cache/ai-theorist/fresh-ge-internal-strain/. --install copies the two UPF files and the
SSSP LICENSE.txt (as SSSP_LICENSE.txt) into environment/lab/ (lab image only; the agent image never sees them).
"""
import hashlib
import json
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/rcyfm-68h65'
ARCHIVE = ['SSSP_1.3.0_PBE_precision.tar.gz', 'SSSP_1.3.0_PBE_precision.json', 'LICENSE.txt']
PINNED_MD5 = {'SSSP_1.3.0_PBE_precision.tar.gz': 'fde94756886f32ada7bf597547557eb5',
              'SSSP_1.3.0_PBE_precision.json': '1692c5c9ce89e1c7c783f8f0eee0cbfa',
              'LICENSE.txt': 'c0f9df05fb4cd4921f4caa25606815ef'}
GBRV = ('ge_pbe_v1.4.uspp.F.UPF', '9c9eaa91e581c3f09632fb3098b2c6b2',
        '1b95046ab0f3228e0404719ca66d0309af5a62313be0b22171b963566c90430b')
PSL_URL = 'https://pseudopotentials.quantum-espresso.org/upf_files/'
PSL = ('Ge.pbe-n-rrkjus_psl.1.0.0.UPF', 'a23a7f3f7456476e534e54484b27d8e8',
       'b756d20fe98992cb3afef5fcde78515720a89336d412c13b0ceddf794d54cfcd')
CACHE = Path.home()/'.cache/ai-theorist/fresh-ge-internal-strain'
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-ge-internal-strain/environment/lab'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def download(url, target):
    print('downloading', url)
    with urllib.request.urlopen(url, timeout=600) as r, open(target, 'wb') as out:
        shutil.copyfileobj(r, out)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(RECORD, timeout=60) as r:
        entries = json.load(r)['files']['entries']
    for name in ARCHIVE:
        md5 = entries[name]['checksum'].split(':', 1)[1]
        assert md5 == PINNED_MD5[name], f'record MD5 of {name} changed: {md5}'
        target = CACHE/name
        if not target.exists() or digest(target, 'md5') != md5:
            download(f'{RECORD}/files/{name}/content', target)
        assert digest(target, 'md5') == md5, f'MD5 mismatch for {name}'
        print(name, target.stat().st_size, 'md5', md5)
    meta = json.loads((CACHE/'SSSP_1.3.0_PBE_precision.json').read_text())
    name, md5, sha = GBRV
    assert meta['Ge']['filename'] == name and meta['Ge']['md5'] == md5, meta['Ge']
    with tarfile.open(CACHE/'SSSP_1.3.0_PBE_precision.tar.gz') as tar:
        member = next(m for m in tar.getmembers() if Path(m.name).name == name)
        with tar.extractfile(member) as src, open(CACHE/name, 'wb') as dst:
            shutil.copyfileobj(src, dst)
    name, md5, sha = PSL
    if not (CACHE/name).exists() or digest(CACHE/name, 'sha256') != sha:
        download(PSL_URL + name, CACHE/name)
    for name, md5, sha in (GBRV, PSL):
        assert digest(CACHE/name, 'md5') == md5, f'MD5 mismatch for {name}'
        assert digest(CACHE/name, 'sha256') == sha, f'SHA256 mismatch for {name}'
        print(name, (CACHE/name).stat().st_size, 'md5', md5, 'sha256', sha)
    print('SSSP Ge entry:', meta['Ge'])
    if '--install' in sys.argv:
        for name, _, _ in (GBRV, PSL):
            shutil.copy2(CACHE/name, LAB/name)
        shutil.copy2(CACHE/'LICENSE.txt', LAB/'SSSP_LICENSE.txt')
        print('installed into', LAB)


if __name__ == '__main__':
    main()
