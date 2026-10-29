"""Fetch contiguous chunks of each run's nested traj/trajectory.tar (T <= 1100 K, all N) by range reads.

    python3 -B fetch_frames.py     # -> ~/.cache/ai-theorist/dense-hydrogen-llt-crossover/frames/<tar>/<rs>/chunk<k>.bin + meta.json

Six chunks per run, centred at relative positions (k + 0.5)/6 of the nested tar; each holds about ten
consecutive dumps (20 fs apart). The SHA256 of every chunk is written to (or checked against) members_frames.json.
Only these byte ranges are read; the trajectories (0.08-0.6 GB each) are never downloaded whole.
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
CHUNK = {128: 100000, 256: 190000, 512: 380000, 2048: 1200000}
NCHUNK = 6
TARS = [f'N_{n}_T{t}K.tar' for n in (128, 256, 512, 2048) for t in (900, 950, 1000, 1100)]


def main():
    pinned_path = HERE/'members_frames.json'
    pinned = json.loads(pinned_path.read_text()) if pinned_path.exists() else {}
    record = dict(pinned)
    tars = sys.argv[1:] or TARS
    for tar in tars:
        n = int(tar.split('_')[1])
        listing = json.loads((CACHE/f'list_{tar}.json').read_text())
        jobs = []
        for e in listing:
            if e['name'].endswith('/traj/trajectory.tar'):
                rs = e['name'].split('/')[1]
                for k in range(NCHUNK):
                    off = e['offset'] + int((k + 0.5) / NCHUNK * (e['size'] - CHUNK[n]))
                    off -= (off - e['offset']) % 512
                    end = min(off + CHUNK[n], e['offset'] + e['size'])
                    jobs.append((rs, k, e, off, end))

        def get(job):
            rs, k, e, off, end = job
            d = CACHE/'frames'/tar/rs
            out = d/f'chunk{k}.bin'
            key = f'{tar}:{e["name"]}:{off - e["offset"]}-{end - e["offset"]}'
            if out.exists() and out.stat().st_size == end - off:
                blob = out.read_bytes()
            else:
                blob = RT.fetch(tar, off, end - 1)
                d.mkdir(parents=True, exist_ok=True)
                out.write_bytes(blob)
            (d/f'chunk{k}.json').write_text(json.dumps({'base': e['offset'], 'offset': off, 'size': e['size']}))
            digest = hashlib.sha256(blob).hexdigest()
            if key in pinned and pinned[key] != digest:
                raise RuntimeError('SHA256 mismatch ' + key)
            return key, digest
        with ThreadPoolExecutor(3) as pool:
            for key, digest in pool.map(get, jobs):
                record[key] = digest
        print(tar, len(jobs), 'chunks', flush=True)
        pinned_path.write_text(json.dumps(dict(sorted(record.items())), indent=0) + '\n')


if __name__ == '__main__':
    main()
