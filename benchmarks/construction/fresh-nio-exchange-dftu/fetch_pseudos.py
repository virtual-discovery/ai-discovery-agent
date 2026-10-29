"""Refetch the two SSSP 1.3.0 PBE-precision pseudopotentials used by fresh-nio-exchange-dftu.

    python3 -B benchmarks/construction/fresh-nio-exchange-dftu/fetch_pseudos.py [--install]

Downloads SSSP_1.3.0_PBE_precision.tar.gz, its metadata JSON and LICENSE.txt from Materials Cloud Archive
record rcyfm-68h65 (DOI 10.24435/materialscloud:f3-ym; Prandini et al., npj Comput. Mater. 4, 72 (2018))
into ~/.cache/ai-theorist/fresh-nio-exchange-dftu/, checks the archive MD5s (published by the record API),
extracts the Ni (GBRV-1.4, GPL v3) and O (Pslibrary 0.3.1, GPL v2+) files, checks each against the MD5 in
the SSSP metadata and the SHA256 pinned in pseudo_checksums.json. --install copies the two UPF files and the
archive LICENSE.txt (as SSSP_LICENSE.txt) into the task's environment/lab/ (lab image only).
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
UPF = {'Ni': 'ni_pbe_v1.4.uspp.F.UPF', 'O': 'O.pbe-n-kjpaw_psl.0.1.UPF'}
HERE = Path(__file__).resolve().parent
PINNED = json.loads((HERE/'pseudo_checksums.json').read_text())
CACHE = Path.home()/'.cache/ai-theorist/fresh-nio-exchange-dftu'
LAB = HERE.parents[1]/'tasks/fresh-nio-exchange-dftu/environment/lab'


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
        for element, filename in UPF.items():
            assert meta[element]['filename'] == filename, (element, meta[element])
            data = tar.extractfile(members[filename]).read()
            (CACHE/filename).write_bytes(data)
            assert hashlib.md5(data).hexdigest() == meta[element]['md5'] == PINNED[filename]['md5'], filename
            assert hashlib.sha256(data).hexdigest() == PINNED[filename]['sha256'], filename
            print(filename, 'md5', meta[element]['md5'], 'sha256', PINNED[filename]['sha256'],
                  'family', meta[element]['pseudopotential'], 'SSSP cutoffs (Ry)', meta[element]['cutoff_wfc'],
                  meta[element]['cutoff_rho'])
    if '--install' in sys.argv:
        for filename in UPF.values():
            shutil.copy2(CACHE/filename, LAB/filename)
        shutil.copy2(CACHE/'LICENSE.txt', LAB/'SSSP_LICENSE.txt')
        print('installed into', LAB)


if __name__ == '__main__':
    main()
