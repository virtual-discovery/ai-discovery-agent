"""Fetch the small text members (input decks, LAMMPS logs, slurm outputs) of every classical-MD tar.

    python3 -B fetch_thermo.py [TARNAME ...]     # -> ~/.cache/ai-theorist/dense-hydrogen-llt-crossover/raw/<tar>/...

Source: Zenodo record 14790037 (Tenti et al., CC BY 4.0). Members are read by HTTP range requests
(remotetar.py); trajectories are never downloaded here. Listings are cached as list_<tar>.json; the
SHA256 of every fetched member is written to members_thermo.json (or checked against it when present).
"""
import hashlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import remotetar as RT

HERE = Path(__file__).resolve().parent
CACHE = Path(os.path.expanduser('~/.cache/ai-theorist/dense-hydrogen-llt-crossover'))
TARS = [f'N_{n}_T{t}K.tar' for n, ts in {128: [900, 950, 1000, 1100, 1200, 1300, 1400],
                                         256: [900, 950, 1000, 1100], 512: [900, 950, 1000, 1100],
                                         2048: [900, 950, 1000, 1100, 1200, 1400, 1600, 1800]}.items() for t in ts]


def listing(tar):
    path = CACHE/f'list_{tar}.json'
    if not path.exists() or path.stat().st_size == 0:
        path.write_text(json.dumps(RT.listing(tar)))
    return json.loads(path.read_text())


def wanted(e):
    n = e['name']
    return e['type'] == '0' and '/traj' not in n and (n.endswith(('input.lammps', 'log.lammps')) or '/slurm-' in n)


def main():
    tars = sys.argv[1:] or TARS
    pinned_path = HERE/'members_thermo.json'
    pinned = json.loads(pinned_path.read_text()) if pinned_path.exists() else {}
    record = dict(pinned)
    for tar in tars:
        entries = [e for e in listing(tar) if wanted(e)]

        def get(e):
            out = CACHE/'raw'/tar/e['name']
            key = f"{tar}:{e['name']}"
            if out.exists() and out.stat().st_size == e['size']:
                data = out.read_bytes()
            else:
                data = RT.extract(tar, e)
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(data)
            digest = hashlib.sha256(data).hexdigest()
            if key in pinned and pinned[key] != digest:
                raise RuntimeError('SHA256 mismatch ' + key)
            return key, digest
        with ThreadPoolExecutor(3) as pool:
            for key, digest in pool.map(get, entries):
                record[key] = digest
        print(tar, len(entries), 'members', flush=True)
    pinned_path.write_text(json.dumps(dict(sorted(record.items())), indent=0) + '\n')


if __name__ == '__main__':
    main()
