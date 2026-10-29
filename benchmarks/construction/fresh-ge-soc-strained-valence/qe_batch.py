"""Run a batch of lab requests through task_lab.run, inside the runtime image (resumable).

    python3 qe_batch.py BATCH.json OUT.jsonl [--cpuset 12,13,28,29] [--workdir DIR]

BATCH.json is a list of {"tag": str, "request": {...}}. Each request goes through the task's own
task_lab.validate/run (the code path the agent uses), sequentially, with mpirun -np 2, in
docker run --rm --network none --cpuset-cpus CPUSET --user UID:GID ai-theorist-benchmark-runtime:20261007.
One JSON line per request is appended to OUT.jsonl (tag, request, charge, result or error); tags already
present in OUT.jsonl are skipped, so an interrupted batch can be resumed.
"""
import argparse
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LAB = ROOT/'benchmarks/tasks/fresh-ge-soc-strained-valence/environment/lab'
IMAGE = 'ai-theorist-benchmark-runtime:20261007'

INNER = r'''
import json, sys
sys.path.insert(0, '/lab')
import task_lab
batch = json.load(open(sys.argv[1]))
out = sys.argv[2]
try:
    done = {json.loads(l)['tag'] for l in open(out) if l.strip()}
except FileNotFoundError:
    done = set()
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
    print(item['tag'], row.get('charge'), row.get('result', {}).get('wall_seconds'), row.get('error', ''), flush=True)
'''


def main():
    p = argparse.ArgumentParser()
    p.add_argument('batch')
    p.add_argument('out')
    p.add_argument('--cpuset', default='12,13,28,29')
    p.add_argument('--workdir', default=str(ROOT/'runs/build-fresh-ge-soc-strained-valence/jobs'))
    a = p.parse_args()
    batch, out, work = Path(a.batch).resolve(), Path(a.out).resolve(), Path(a.workdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work/'inner.py').write_text(INNER)
    mounts = []
    for path in {batch.parent, out.parent, work}:
        mounts += ['-v', f'{path}:{path}']
    cmd = ['docker', 'run', '--rm', '--network', 'none', '--cpuset-cpus', a.cpuset, '--user', f'{os.getuid()}:{os.getgid()}',
           '-e', 'HOME=/tmp', '-e', 'OMPI_MCA_btl_vader_single_copy_mechanism=none', '-e', 'PYTHONDONTWRITEBYTECODE=1',
           '-v', f'{LAB}:/lab:ro', *mounts, IMAGE, 'python3', '-B', str(work/'inner.py'), str(batch), str(out), str(work)]
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
