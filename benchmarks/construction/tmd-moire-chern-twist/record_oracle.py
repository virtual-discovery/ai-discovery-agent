"""Record the oracle's gate run (runs/check/tmd-moire-chern-twist) into oracle_run.json for write_design.py.

    python3 -B benchmarks/construction/tmd-moire-chern-twist/record_oracle.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE = HERE.parents[2]/'runs/check/tmd-moire-chern-twist/lab-results/state.json'
state = json.loads(STATE.read_text())
sub = state['submission']
out = {'evaluations': len(state['trace']), 'charge_used': state['used'], 'budget': state['budget'],
       'requests': [dict(r['request'], charge=r['charge']) for r in state['trace']],
       'total_job_seconds': round(sum(r['ended_epoch'] - r['started_epoch'] for r in state['trace']), 1),
       'submission': {k: v for k, v in sub['payload'].items() if k != 'evidence_ids'},
       'accepted': sub['success'], 'errors_meV': sub['report'].get('errors_meV')}
(HERE/'oracle_run.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: out[k] for k in ('evaluations', 'charge_used', 'budget', 'accepted', 'errors_meV', 'total_job_seconds')}, indent=1))
