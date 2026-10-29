"""Run the reference calculations of reference_plan.PLAN through task_lab.run (in the runtime image).

    python3 -B benchmarks/construction/fresh-nio-exchange-dftu/run_reference.py [--cpuset 8,9,24,25]

Appends one line per tag to reference_runs.jsonl in this folder (resumable: finished tags are skipped).
All requests run via qe_batch.py (docker run --rm --cpuset-cpus CPUSET --user UID:GID, mpirun -np 2).
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
RUNS = HERE.parents[2]/'runs/build-fresh-nio-exchange-dftu'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cpuset', default='8,9,24,25')
    p.add_argument('--only', help='comma-separated tag prefixes to run')
    a = p.parse_args()
    done = {json.loads(l)['tag'] for l in OUT.read_text().splitlines() if l.strip()} if OUT.exists() else set()
    todo = [{'tag': t, 'request': r} for t, r in PLAN.items() if t not in done
            and (not a.only or any(t.startswith(x) for x in a.only.split(',')))]
    print(f'{len(todo)} reference calculations to run')
    if todo:
        batch = RUNS/'reference_batch.json'
        batch.parent.mkdir(parents=True, exist_ok=True)
        batch.write_text(json.dumps(todo, indent=1))
        subprocess.run([sys.executable, '-B', str(HERE/'qe_batch.py'), str(batch), str(OUT), '--cpuset', a.cpuset,
                        '--workdir', str(RUNS/'reference_jobs')], check=True)


if __name__ == '__main__':
    main()
