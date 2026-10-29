"""Refetch the two SSSP 1.3.0 PBE-precision pseudopotentials used by fresh-inp-polar-phonons.

    python3 benchmarks/construction/fresh-inp-polar-phonons/fetch_pseudos.py [--install]

Downloads SSSP_1.3.0_PBE_precision.tar.gz and its metadata JSON from Materials Cloud Archive record
rcyfm-68h65 (DOI 10.24435/materialscloud:f3-ym, CC BY 4.0) into ~/.cache/ai-theorist/fresh-inp-polar-phonons/,
checks the archive MD5s, extracts In and P, checks each UPF against the md5 in the SSSP metadata and against
the SHA256 pinned below. --install copies the two UPF files into the task's environment/lab/.
"""
import hashlib
import json
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/rcyfm-68h65/files/'
FILES = {'SSSP_1.3.0_PBE_precision.tar.gz': 'fde94756886f32ada7bf597547557eb5',
         'SSSP_1.3.0_PBE_precision.json': '1692c5c9ce89e1c7c783f8f0eee0cbfa',
         'LICENSE.txt': 'c0f9df05fb4cd4921f4caa25606815ef'}
UPF = {'In': ('In.pbe-dn-rrkjus_psl.0.2.2.UPF', None),
       'P': ('P.pbe-n-rrkjus_psl.1.0.0.UPF', None)}
SHA256 = json.loads((Path(__file__).with_name('pseudo_checksums.json')).read_text()) \
    if Path(__file__).with_name('pseudo_checksums.json').exists() else {}
CACHE = Path.home()/'.cache/ai-theorist/fresh-inp-polar-phonons'
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-inp-polar-phonons/environment/lab'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    for name, md5 in FILES.items():
        target = CACHE/name
        if not target.exists() or digest(target, 'md5') != md5:
            print('downloading', name)
            with urllib.request.urlopen(RECORD + name + '/content', timeout=600) as r, open(target, 'wb') as out:
                shutil.copyfileobj(r, out)
        assert digest(target, 'md5') == md5, name
    meta = json.loads((CACHE/'SSSP_1.3.0_PBE_precision.json').read_text())
    with tarfile.open(CACHE/'SSSP_1.3.0_PBE_precision.tar.gz') as tar:
        members = {Path(m.name).name: m for m in tar.getmembers() if m.isfile()}
        for element, (filename, _) in UPF.items():
            assert meta[element]['filename'] == filename, (element, meta[element])
            data = tar.extractfile(members[filename]).read()
            (CACHE/filename).write_bytes(data)
            assert hashlib.md5(data).hexdigest() == meta[element]['md5'], filename
            sha = hashlib.sha256(data).hexdigest()
            if SHA256:
                assert SHA256[filename] == sha, filename
            print(filename, 'md5', meta[element]['md5'], 'sha256', sha, 'cutoff_wfc', meta[element]['cutoff_wfc'],
                  'cutoff_rho', meta[element]['cutoff_rho'], 'pseudopotential', meta[element]['pseudopotential'])
    if '--install' in sys.argv:
        for filename, _ in UPF.values():
            shutil.copy2(CACHE/filename, LAB/filename)
        print('installed into', LAB)


if __name__ == '__main__':
    main()
