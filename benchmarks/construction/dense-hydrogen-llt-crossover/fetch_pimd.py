"""Fetch every PIMD thermodata.npy (nested tars PIMD/900K.tar and PIMD/1000K.tar inside PIMD.tar) by range reads.

    python3 -B fetch_pimd.py    # -> ~/.cache/ai-theorist/dense-hydrogen-llt-crossover/raw/PIMD/<T>/<rs>.<run>/thermodata.npy
SHA256 of each member is written to (or checked against) members_pimd.json.
"""
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import remotetar as RT

HERE = Path(__file__).resolve().parent
CACHE = Path(os.path.expanduser('~/.cache/ai-theorist/dense-hydrogen-llt-crossover'))


def main():
    outer = RT.listing('PIMD.tar')
    pinned_path = HERE/'members_pimd.json'
    pinned = json.loads(pinned_path.read_text()) if pinned_path.exists() else {}
    record = {}
    readme = [e for e in outer if e['name'] == 'PIMD/README'][0]
    data = RT.extract('PIMD.tar', readme)
    (CACHE/'raw/PIMD').mkdir(parents=True, exist_ok=True)
    (CACHE/'raw/PIMD/README').write_bytes(data)
    record['PIMD/README'] = hashlib.sha256(data).hexdigest()
    for e in outer:
        if not e['name'].endswith('K.tar'):
            continue
        inner = [x for x in RT.listing('PIMD.tar', e['offset'], e['offset'] + e['size']) if x['type'] == '0']

        def get(x):
            out = CACHE/'raw/PIMD'/x['name']
            if out.exists() and out.stat().st_size == x['size']:
                blob = out.read_bytes()
            else:
                blob = RT.extract('PIMD.tar', x)
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(blob)
            key = e['name'] + ':' + x['name']
            digest = hashlib.sha256(blob).hexdigest()
            if key in pinned and pinned[key] != digest:
                raise RuntimeError('SHA256 mismatch ' + key)
            return key, digest
        with ThreadPoolExecutor(3) as pool:
            for key, digest in pool.map(get, inner):
                record[key] = digest
        print(e['name'], len(inner), flush=True)
    pinned_path.write_text(json.dumps(dict(sorted(record.items())), indent=0) + '\n')


if __name__ == '__main__':
    main()
