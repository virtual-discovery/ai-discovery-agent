"""Run the oracle's own request chain through task_lab.run (stage C of the reference calculations).

    python3 -B oracle_chain.py [--cpuset 12,13,28,29]

The oracle chooses each geometry from earlier results (a0, the relaxed c), so its requests cannot be listed in advance.
This loop replays solution/solve.py offline (replay.py's stand-in client) until it asks for a request that is not yet in
reference_runs.jsonl, runs that request with qe_batch.py (runtime image, same task_lab.run), and repeats until the
replay completes. Tags are C_oracle_NN.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[2]/'runs/build-fresh-ge-soc-strained-valence'


def main():
    cpuset = sys.argv[sys.argv.index('--cpuset') + 1] if '--cpuset' in sys.argv else '12,13,28,29'
    for n in range(40):
        proc = subprocess.run([sys.executable, '-B', str(HERE/'replay.py')], capture_output=True, text=True)
        m = re.search(r'REPLAY: request not in reference runs: (\{.*?\})(?:\\n|\n|"|$)', proc.stdout + proc.stderr)
        if not m:
            print(proc.stdout[-3000:])
            return
        req = json.loads(m.group(1).replace('\\"', '"'))
        done = [json.loads(l)['tag'] for l in (HERE/'reference_runs.jsonl').read_text().splitlines() if l.strip()]
        tag = 'C_oracle_%02d' % sum(t.startswith('C_oracle_') for t in done)
        print('running', tag, json.dumps(req), flush=True)
        batch = RUNS/'oracle_chain_batch.json'
        batch.write_text(json.dumps([{'tag': tag, 'request': req}]))
        subprocess.run([sys.executable, '-B', str(HERE/'qe_batch.py'), str(batch), str(HERE/'reference_runs.jsonl'),
                        '--cpuset', cpuset, '--workdir', str(RUNS/'reference_jobs')], check=True)


if __name__ == '__main__':
    main()
