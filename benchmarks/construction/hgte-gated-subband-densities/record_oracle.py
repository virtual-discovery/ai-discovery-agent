"""Record the gate's oracle run (runs/check/<id>) into oracle_run.json.

    python3 -B record_oracle.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[2]/'runs/check/hgte-gated-subband-densities'
state = json.loads((RUN/'lab-results/state.json').read_text())
check = json.loads((HERE/'check.json').read_text())
trace = state['trace']
out = {'evaluations': len(trace), 'charge_used': state['used'], 'budget': state['budget'],
       'requests': [r['request'] for r in trace],
       'lab_job_seconds': round(sum(r['ended_epoch'] - r['started_epoch'] for r in trace), 2),
       'container_seconds': check['summary']['oracle']['container_seconds'],
       'accepted': state['submission']['success'], 'errors': state['submission']['report']['errors']}
(HERE/'oracle_run.json').write_text(json.dumps(out, indent=1) + '\n')
print(out['evaluations'], out['charge_used'], out['accepted'])
