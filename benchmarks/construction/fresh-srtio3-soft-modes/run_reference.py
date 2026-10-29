"""Run the reference calculations of reference_plan through task_lab.run (runtime image, cpuset 10,11,26,27).

    python3 -B benchmarks/construction/fresh-srtio3-soft-modes/run_reference.py --stage A|B [--cpuset CPUS]

Appends one JSON line per tag to reference_runs.jsonl in this folder (resumable: tags already present are skipped).
Stage B needs the stage A results.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import reference_plan  # noqa: E402

OUT = HERE/'reference_runs.jsonl'
RUNS = HERE.parents[2]/'runs/build-fresh-srtio3-soft-modes'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['A', 'B'], required=True)
    p.add_argument('--cpuset', default='10,11,26,27')
    a = p.parse_args()
    done = {json.loads(l)['tag'] for l in OUT.read_text().splitlines() if l.strip()} if OUT.exists() else set()
    plan = reference_plan.plan(a.stage)
    todo = [{'tag': t, 'request': r} for t, r in plan.items() if t not in done]
    print(f'{len(todo)} reference calculations to run (stage {a.stage})', flush=True)
    if todo:
        RUNS.mkdir(parents=True, exist_ok=True)
        batch = RUNS/f'reference_batch_{a.stage}.json'
        batch.write_text(json.dumps(todo, indent=1))
        subprocess.run([sys.executable, '-B', str(HERE/'qe_batch.py'), str(batch), str(OUT), '--cpuset', a.cpuset,
                        '--workdir', str(RUNS/'reference_jobs'), '--name', f'sto-ref-{a.stage}'], check=True)


if __name__ == '__main__':
    main()
