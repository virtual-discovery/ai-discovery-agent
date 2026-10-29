"""Refetch the Al pseudopotential used by fresh-al-surface-workfunction.

    python3 -B benchmarks/construction/fresh-al-surface-workfunction/fetch_pseudo.py [--install]

Downloads SSSP_1.3.0_PBE_precision.tar.gz, its metadata JSON and the record LICENSE from Materials Cloud
Archive record rcyfm-68h65 (DOI 10.24435/materialscloud:f3-ym, CC BY 4.0) into
~/.cache/ai-theorist/fresh-al-surface-workfunction/, checks the archive MD5s published by the record,
extracts Al.pbe-n-kjpaw_psl.1.0.0.UPF, and checks it against the md5 in the SSSP metadata and the SHA256
pinned below (identical to the file used by the v1 task t1-al-fcc-eos-fresh). --install copies the UPF
and the license into the task's environment/lab/.
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
UPF = 'Al.pbe-n-kjpaw_psl.1.0.0.UPF'
UPF_MD5 = 'cfc449ca30b5f3223ec38ddd88ac046d'
UPF_SHA256 = 'fd7b78921e6d0939095b681c328732668eaefd1dee6882fdbc03be463550cc97'
CACHE = Path.home()/'.cache/ai-theorist/fresh-al-surface-workfunction'
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-al-surface-workfunction/environment/lab'


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
    meta = json.loads((CACHE/'SSSP_1.3.0_PBE_precision.json').read_text())['Al']
    assert meta['filename'] == UPF and meta['md5'] == UPF_MD5, meta
    with tarfile.open(CACHE/'SSSP_1.3.0_PBE_precision.tar.gz') as tar:
        member = next(m for m in tar.getmembers() if m.isfile() and Path(m.name).name == UPF)
        data = tar.extractfile(member).read()
    assert hashlib.md5(data).hexdigest() == UPF_MD5 and hashlib.sha256(data).hexdigest() == UPF_SHA256
    (CACHE/UPF).write_bytes(data)
    print(UPF, 'md5', UPF_MD5, 'sha256', UPF_SHA256, 'bytes', len(data), 'SSSP metadata', json.dumps(meta))
    if '--install' in sys.argv:
        shutil.copy2(CACHE/UPF, LAB/UPF)
        shutil.copy2(CACHE/'LICENSE.txt', LAB/'SSSP_LICENSE.txt')
        print('installed into', LAB)


if __name__ == '__main__':
    main()
