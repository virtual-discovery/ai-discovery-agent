"""Refetch the SSSP 1.3.0 PBE-precision lithium pseudopotential used by fresh-li-bcc-soft-shear.

    python3 benchmarks/construction/fresh-li-bcc-soft-shear/fetch_pseudos.py [--install]

Downloads SSSP_1.3.0_PBE_precision.tar.gz, its metadata JSON and the archive LICENSE.txt from Materials Cloud
Archive record rcyfm-68h65 (DOI 10.24435/materialscloud:f3-ym, CC BY 4.0) into
~/.cache/ai-theorist/fresh-li-bcc-soft-shear/, checks the archive MD5s, extracts li_pbe_v1.4.uspp.F.UPF
(GBRV 1.4, GPL v3), checks it against the md5 in the SSSP metadata and against the SHA256 pinned in
pseudo_checksums.json. --install copies the UPF file and the archive LICENSE.txt (as SSSP_LICENSE.txt) into the
task's environment/lab/.
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
ELEMENT, FILENAME = 'Li', 'li_pbe_v1.4.uspp.F.UPF'
SHA256 = json.loads(Path(__file__).with_name('pseudo_checksums.json').read_text())
CACHE = Path.home()/'.cache/ai-theorist/fresh-li-bcc-soft-shear'
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-li-bcc-soft-shear/environment/lab'


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
    meta = json.loads((CACHE/'SSSP_1.3.0_PBE_precision.json').read_text())[ELEMENT]
    assert meta['filename'] == FILENAME, meta
    with tarfile.open(CACHE/'SSSP_1.3.0_PBE_precision.tar.gz') as tar:
        member = next(m for m in tar.getmembers() if m.isfile() and Path(m.name).name == FILENAME)
        data = tar.extractfile(member).read()
    (CACHE/FILENAME).write_bytes(data)
    assert hashlib.md5(data).hexdigest() == meta['md5'], FILENAME
    assert hashlib.sha256(data).hexdigest() == SHA256[FILENAME], FILENAME
    print(FILENAME, 'md5', meta['md5'], 'sha256', SHA256[FILENAME], 'family', meta['pseudopotential'],
          'cutoff_wfc', meta['cutoff_wfc'], 'cutoff_rho', meta['cutoff_rho'])
    if '--install' in sys.argv:
        shutil.copy2(CACHE/FILENAME, LAB/FILENAME)
        shutil.copy2(CACHE/'LICENSE.txt', LAB/'SSSP_LICENSE.txt')
        print('installed into', LAB)


if __name__ == '__main__':
    main()
