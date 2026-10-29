"""Offline replay of solve.py and the ablations against the reference calculations (seconds, no QE).

    .venv/bin/python -B benchmarks/construction/fresh-ge-internal-strain/replay.py

A stand-in lab_client.py serves every request the oracle makes from reference_runs.jsonl (matched on the
canonical request of task_lab.normalize, charged with task_lab.validate, budget enforced) and grades the
submission with task_lab.grade. Then, as check_task does, every ablation must fail, every robustness variation
and the 4-significant-digit rounding must pass. Requests missing from the reference runs abort the replay.
"""
import copy
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-ge-internal-strain'
LAB = TASK/'environment/lab'

CLIENT = r'''
import json, os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.environ['REPLAY_LAB'])
import task_lab
state_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'state.json')
state = json.load(open(state_path))
budget = json.load(open(os.path.join(os.environ['REPLAY_LAB'], 'public.json')))['solver_budget']
rows = {}
for line in open(os.environ['REPLAY_RUNS']):
    row = json.loads(line)
    if 'result' in row:
        try:
            rows.setdefault(json.dumps(task_lab.normalize(row['request']), sort_keys=True), row)
        except ValueError:
            pass
def view(rec):
    return dict(rec, cumulative_charge=state['used'], budget=budget)
endpoint = sys.argv[1]
payload = json.loads(sys.argv[2]) if len(sys.argv) > 2 else None
try:
    if endpoint == 'evaluate':
        charge = task_lab.validate(payload)
        if state['used'] + charge > budget:
            raise ValueError('Budget exhausted')
        key = json.dumps(task_lab.normalize(payload), sort_keys=True)
        if key not in rows:
            print(json.dumps({'error': 'REPLAY: request not in reference runs: ' + key})); sys.exit(1)
        state['used'] += charge
        rec = {'id': len(state['trace']) + 1, 'status': 'complete', 'charge': charge, 'request': payload,
               'result': rows[key]['result'], 'tag': rows[key]['tag']}
        state['trace'].append(rec)
        out = view(rec)
    elif endpoint.startswith('job/'):
        out = view(state['trace'][int(endpoint.split('/')[1]) - 1])
    elif endpoint == 'submit':
        ok, report = task_lab.grade(payload, state['trace'], [])
        state['submission'] = {'payload': payload, 'success': ok, 'report': report}
        out = {'episode_finished': True, 'status': 'complete', 'accepted': ok}
    elif endpoint == 'submission':
        out = {'episode_finished': True, 'status': 'complete', 'accepted': state['submission']['success']}
    else:
        out = {'used': state['used'], 'budget': budget}
except ValueError as exc:
    out = {'error': str(exc)}
json.dump(state, open(state_path, 'w'))
print(json.dumps(out))
'''


def main():
    app = Path(tempfile.mkdtemp(prefix='ge-replay-'))
    try:
        (app/'lab_client.py').write_text(CLIENT)
        (app/'state.json').write_text(json.dumps({'used': 0, 'trace': [], 'submission': None}))
        env = dict(os.environ, APP_DIR=str(app), REPLAY_LAB=str(LAB), REPLAY_RUNS=str(HERE/'reference_runs.jsonl'),
                   PYTHONDONTWRITEBYTECODE='1')
        proc = subprocess.run([sys.executable, '-B', str(TASK/'solution/solve.py')], env=env, capture_output=True, text=True)
        print(proc.stdout[-2000:], proc.stderr[-3000:])
        state = json.loads((app/'state.json').read_text())
        sub = state['submission']
        print('charge used', state['used'], 'evaluations', len(state['trace']))
        print('tags:', [r['tag'] for r in state['trace']])
        if not sub:
            raise SystemExit('no submission')
        print('oracle accepted:', sub['success'])
        print(json.dumps(sub['report'], indent=1))
        sys.path.insert(0, str(LAB))
        import task_lab
        spec = importlib.util.spec_from_file_location('ablations', TASK/'solution/ablations.py')
        ab = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ab)
        payload, trace = sub['payload'], state['trace']
        grade = lambda p, t: task_lab.grade(copy.deepcopy(p), copy.deepcopy(t), [])     # noqa: E731
        for dp, name, mutate in ab.ABLATIONS:
            p2, t2 = mutate(copy.deepcopy(payload), copy.deepcopy(trace))
            ok, rep = grade(p2, t2)
            failed = [k for k, v in rep.get('checks', {}).items() if not v]
            print(f'ablation {dp} {name}: passes={ok} failing={failed}')
        rounded = json.loads(json.dumps(payload), parse_float=lambda s: float('%.4g' % float(s)))
        print('rounding 4 sig passes:', grade(rounded, trace)[0])
        for name, mutate in ab.ROBUSTNESS:
            ok, rep = grade(*mutate(copy.deepcopy(payload), copy.deepcopy(trace)))
            print(f'robustness {name}: passes={ok}', '' if ok else [k for k, v in rep['checks'].items() if not v])
    finally:
        shutil.rmtree(app, ignore_errors=True)


if __name__ == '__main__':
    main()
