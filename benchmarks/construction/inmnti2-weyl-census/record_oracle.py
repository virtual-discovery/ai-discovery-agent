"""Record the oracle's gate run (runs/check/inmnti2-weyl-census) into oracle_run.json for write_design.py.

    python3 -B benchmarks/construction/inmnti2-weyl-census/record_oracle.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE = HERE.parents[2]/'runs/check/inmnti2-weyl-census/lab-results/state.json'
state = json.loads(STATE.read_text())
sub = state['submission']
out = {'evaluations': len(state['trace']), 'charge_used': state['used'], 'budget': state['budget'],
       'requests': [{k: v for k, v in r['request'].items() if not k.startswith(('k_', 'seeds_', 'spheres'))} for r in state['trace']],
       'seconds': {r['id']: round(r['ended_epoch'] - r['started_epoch'], 1) for r in state['trace']},
       'submission': {k: v for k, v in sub['payload'].items() if k != 'evidence_ids'},
       'accepted': sub['success'], 'errors': {k: round(v, 6) for k, v in sub['report']['errors'].items()}}
(HERE/'oracle_run.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: out[k] for k in ('evaluations', 'charge_used', 'budget', 'accepted', 'errors', 'seconds')}, indent=1))
