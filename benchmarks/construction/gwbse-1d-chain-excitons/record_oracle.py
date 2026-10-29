"""Record the oracle's gate run (runs/check/gwbse-1d-chain-excitons) into oracle_run.json for write_design.py.

    cd benchmarks/construction/gwbse-1d-chain-excitons && python3 -B record_oracle.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE = HERE.parents[2]/'runs/check/gwbse-1d-chain-excitons/lab-results/state.json'
state = json.loads(STATE.read_text())
sub = state['submission']
out = {'evaluations': len(state['trace']), 'charge_used': state['used'], 'budget': state['budget'],
       'requests': [r['request'] for r in state['trace']], 'submission': sub['payload'],
       'accepted': sub['success'], 'errors': {k: round(v, 4) for k, v in sub['report']['errors_eV'].items()}}
(HERE/'oracle_run.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: out[k] for k in ('evaluations', 'charge_used', 'budget', 'accepted', 'errors')}, indent=1))
