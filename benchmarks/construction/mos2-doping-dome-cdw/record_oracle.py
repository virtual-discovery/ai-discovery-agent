"""Copy the gate's oracle run (runs/check/mos2-doping-dome-cdw) into oracle_run.json for write_design.py.

    python3 -B record_oracle.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[2]/'runs/check/mos2-doping-dome-cdw'


def main():
    state = json.loads((RUN/'lab-results/state.json').read_text())
    check = json.loads((HERE/'check.json').read_text())['summary']['oracle']
    out = {'evaluations': len(state['trace']), 'charge_used': state['used'], 'budget': state['budget'],
           'requests': [r['request']['kind'] for r in state['trace']],
           'errors': state['submission']['report']['errors'], 'container_seconds': check.get('container_seconds')}
    (HERE/'oracle_run.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out))


if __name__ == '__main__':
    main()
