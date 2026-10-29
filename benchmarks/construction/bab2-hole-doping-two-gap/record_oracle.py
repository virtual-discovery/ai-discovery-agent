"""Record the oracle's gate run (runs/check/bab2-hole-doping-two-gap) into oracle_run.json for write_design.py.

    cd benchmarks/construction/bab2-hole-doping-two-gap && python3 -B record_oracle.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE = HERE.parents[2]/'runs/check/bab2-hole-doping-two-gap/lab-results/state.json'
state = json.loads(STATE.read_text())
sub = state['submission']
out = {'evaluations': len(state['trace']), 'charge_used': state['used'], 'budget': state['budget'],
       'requests': [r['request'] for r in state['trace']], 'submission': {k: v for k, v in sub['payload'].items() if k != 'evidence_ids'},
       'accepted': sub['success'], 'errors': {k: round(v, 5) for k, v in sub['report']['errors'].items()}}
(HERE/'oracle_run.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: out[k] for k in ('evaluations', 'charge_used', 'budget', 'accepted', 'errors')}, indent=1))
