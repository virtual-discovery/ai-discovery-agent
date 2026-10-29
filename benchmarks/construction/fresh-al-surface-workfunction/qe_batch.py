"""Run lab requests through the task's own QE code path inside the runtime image, on this builder's cpuset.

    python3 -B qe_batch.py BATCH.json [--via-run]

BATCH.json is a list of {"name": ..., "request": {...}}. Each request runs in a fresh container with
environment/lab mounted read-only; results go to ~/.cache/ai-theorist/fresh-al-surface-workfunction/runs/NAME/
(result.json plus the pw.x/pp.x inputs and outputs). Existing results are skipped. With --via-run the request is
passed through task_lab.validate() and task_lab.run() exactly as the lab service does; otherwise slab_qe.calculate()
is called directly (exploration requests outside the agent-visible ranges, e.g. reference-only thick slabs).
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

CPUSET = os.environ.get('CPUSET', '8,9,24,25')
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-al-surface-workfunction/environment/lab'
RUNS = Path.home()/'.cache/ai-theorist/fresh-al-surface-workfunction/runs'
IMAGE = 'ai-theorist-benchmark-runtime:20261007'
DIRECT = ("import json,sys; sys.path.insert(0,'/lab'); import slab_qe; r=slab_qe.calculate(json.loads(sys.argv[1]), sys.argv[2]); "
          "open(sys.argv[2]+'/result.json','w').write(json.dumps(r))")
VIA_RUN = ("import json,sys; sys.path.insert(0,'/lab'); import task_lab; q=json.loads(sys.argv[1]); c=task_lab.validate(q); "
           "r=task_lab.run(q, sys.argv[2]); r['_charge']=c; open(sys.argv[2]+'/result.json','w').write(json.dumps(r))")


def main():
    batch = json.loads(Path(sys.argv[1]).read_text())
    code = VIA_RUN if '--via-run' in sys.argv else DIRECT
    for item in batch:
        out = RUNS/item['name']
        if (out/'result.json').exists():
            continue
        out.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        proc = subprocess.run(['docker', 'run', '--rm', '--network', 'none', '--cpuset-cpus', CPUSET,
                               '--user', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp',
                               '-e', 'OMPI_MCA_btl_vader_single_copy_mechanism=none', '-e', 'PYTHONDONTWRITEBYTECODE=1',
                               '-e', 'LAB_MPI_RANKS=' + str(item.get('ranks', 2)),
                               *sum((['-e', f'{k}={v}'] for k, v in item.get('env', {}).items()), []),
                               '-v', f'{LAB}:/lab:ro', '-v', f'{RUNS}:/runs', IMAGE, 'python3', '-c', code,
                               json.dumps(item['request']), f'/runs/{item["name"]}'], capture_output=True, text=True)
        ok = (out/'result.json').exists()
        print(json.dumps({'name': item['name'], 'ok': ok, 'seconds': round(time.time() - t0, 1),
                          'stderr': None if ok else proc.stderr[-600:]}), flush=True)


if __name__ == '__main__':
    main()
