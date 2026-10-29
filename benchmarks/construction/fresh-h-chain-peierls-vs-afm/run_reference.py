"""Run the reference calculations of reference_plan.PLAN through task_lab.run (in the runtime image).

    python3 benchmarks/construction/fresh-h-chain-peierls-vs-afm/run_reference.py [--seed EXPLORE.jsonl]

Writes reference_runs.jsonl in this folder (one line per tag; resumable). With --seed, rows of an
earlier exploration file whose request is identical to a planned request (same canonical request,
same task_lab code path) are copied instead of being recomputed; their tags (and source tags) are recorded in
reference_seeded.json. All other requests run via qe_batch.py (docker, cpuset 10,11,26,27, -np 2).
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from reference_plan import PLAN  # noqa: E402

OUT = HERE/'reference_runs.jsonl'


def canonical(req):
    return json.dumps(req, sort_keys=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--seed')
    p.add_argument('--cpuset', default='10,11,26,27')
    a = p.parse_args()
    done = {json.loads(l)['tag'] for l in OUT.read_text().splitlines() if l.strip()} if OUT.exists() else set()
    if a.seed:
        by_request = {}
        for line in Path(a.seed).read_text().splitlines():
            row = json.loads(line)
            if 'result' in row and 'band_gap_eV' in row['result']:   # rows made with the final parser only
                by_request.setdefault(canonical(row['request']), row)
        log = HERE/'reference_seeded.json'
        seeded = json.loads(log.read_text()) if log.exists() else {}
        with OUT.open('a') as f:
            for tag, req in PLAN.items():
                row = by_request.get(canonical(req))
                if tag not in done and row:
                    f.write(json.dumps({'tag': tag, 'request': req, 'charge': row['charge'], 'result': row['result'],
                                        'seeded_from': row['tag']}) + '\n')
                    seeded[tag] = row['tag']
                    done.add(tag)
        log.write_text(json.dumps(seeded, indent=1) + '\n')
        print(f'{len(seeded)} seeded rows in total')
    todo = [{'tag': t, 'request': r} for t, r in PLAN.items() if t not in done]
    print(f'{len(todo)} reference calculations to run')
    if todo:
        batch = HERE.parents[2]/'runs/build-fresh-h-chain-peierls-vs-afm/reference_batch.json'
        batch.parent.mkdir(parents=True, exist_ok=True)
        batch.write_text(json.dumps(todo, indent=1))
        subprocess.run([sys.executable, str(HERE/'qe_batch.py'), str(batch), str(OUT), '--cpuset', a.cpuset,
                        '--workdir', str(HERE.parents[2]/'runs/build-fresh-h-chain-peierls-vs-afm/reference_jobs')],
                       check=True)


if __name__ == '__main__':
    main()
