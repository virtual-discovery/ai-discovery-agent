"""Copy the gate's oracle run summary (runs/check/hydride-coulomb-mustar) into oracle_run.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
work = os.path.join(HERE, '../../../runs/check/hydride-coulomb-mustar')
state = json.load(open(os.path.join(work, 'lab-results/state.json')))
rc = open(os.path.join(work, 'oracle.rc')).read()
sub = state['submission']
out = {'evaluations': len(state['trace']), 'charge_used': state['used'], 'budget': state['budget'],
       'oracle_seconds': int(rc.split('seconds=')[1].split()[0]),
       'requests': [{'id': r['id'], 'kind': r['request']['kind'], 'charge': r['charge'],
                     'seconds': round(r['ended_epoch'] - r['started_epoch'], 1),
                     **{k: r['request'][k] for k in ('system', 'pressure_GPa', 'solver') if k in r['request']}}
                    for r in state['trace']],
       'payload': sub['payload'], 'report': sub['report']}
json.dump(out, open(os.path.join(HERE, 'oracle_run.json'), 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('evaluations', 'charge_used', 'budget', 'oracle_seconds')}))
