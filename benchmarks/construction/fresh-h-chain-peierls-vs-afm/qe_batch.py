"""Run a batch of lab SCF requests through task_lab.run, inside the runtime image.

    python3 qe_batch.py BATCH.json OUT.jsonl [--cpuset 10,11,26,27]

BATCH.json is a list of {"tag": str, "request": {...}}. Each request runs through the task's own
task_lab.normalize/run (the code path the agent uses), sequentially, with mpirun -np 2, in
docker run --rm --cpuset-cpus CPUSET --user UID:GID ai-theorist-benchmark-runtime:20261007.
Results are appended to OUT.jsonl (one line per request: tag, request, charge, result or error);
requests whose tag is already in OUT.jsonl are skipped, so the batch can be resumed.
"""
import argparse
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LAB = ROOT/'benchmarks/tasks/fresh-h-chain-peierls-vs-afm/environment/lab'
IMAGE = 'ai-theorist-benchmark-runtime:20261007'

INNER = r'''
import json, sys, time, traceback
sys.path.insert(0, '/lab')
import task_lab
batch = json.load(open(sys.argv[1]))
out = sys.argv[2]
done = set()
try:
    done = {json.loads(l)['tag'] for l in open(out) if l.strip()}
except FileNotFoundError:
    pass
for item in batch:
    if item['tag'] in done:
        continue
    row = {'tag': item['tag'], 'request': item['request']}
    try:
        row['charge'] = task_lab.validate(item['request'])
        row['result'] = task_lab.run(item['request'], sys.argv[3] + '/' + item['tag'])
    except Exception as exc:
        row['error'] = '%s: %s' % (type(exc).__name__, exc)
    with open(out, 'a') as f:
        f.write(json.dumps(row) + '\n')
    print(item['tag'], row.get('result', {}).get('wall_seconds'), row.get('error', ''), flush=True)
'''


def main():
    p = argparse.ArgumentParser()
    p.add_argument('batch')
    p.add_argument('out')
    p.add_argument('--cpuset', default='10,11,26,27')
    p.add_argument('--workdir', default=str(ROOT/'runs/build-fresh-h-chain-peierls-vs-afm/jobs'))
    a = p.parse_args()
    batch, out, work = Path(a.batch).resolve(), Path(a.out).resolve(), Path(a.workdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work/'inner.py').write_text(INNER)
    mounts = []
    for path in {batch.parent, out.parent, work}:
        mounts += ['-v', f'{path}:{path}']
    cmd = ['docker', 'run', '--rm', '--network', 'none', '--cpuset-cpus', a.cpuset, '--user', f'{os.getuid()}:{os.getgid()}',
           '-e', 'HOME=/tmp', '-e', 'OMPI_MCA_btl_vader_single_copy_mechanism=none', '-e', 'PYTHONDONTWRITEBYTECODE=1',
           '-v', f'{LAB}:/lab:ro', *mounts, IMAGE, 'python3', str(work/'inner.py'), str(batch), str(out), str(work)]
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
