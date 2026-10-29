"""Development gate for a benchmark task. Writes benchmarks/construction/TASK/check.json.

    python -m benchmarks.tools.check_task TASK_ID [--cpuset 2-5] [--reuse-state]

Checks: package structure and shared files; design (>= 5 decision points, each with an
ablation); rubric (>= 4 numbered items); data size; hidden-value leaks into anything
the agent can read; then, inside the task's own lab image (the runtime's numpy/scipy):
the reference solution must be accepted, every ablation must fail grading, rounding the
submission to 4 significant digits and each listed robustness variation must still pass,
and grading must be deterministic.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import time
from benchmarks.verify import ROOT
from benchmarks.tools.release import TASKS, problems

TEXT = {'.json', '.csv', '.txt', '.md', '.py', '.dat', '.xyz', '.in', '.yaml', '.toml', '.tsv', '.out'}
MAX_FILE_MB, MAX_TASK_MB = 40, 60

RUN_SH = r'''set -u
export LAB_RESULTS=/w/lab-results LAB_URL=http://127.0.0.1:8000 APP_DIR=/w/app PYTHONDONTWRITEBYTECODE=1
if [ "${REUSE:-0}" != 1 ]; then
  python3 /lab/labkit.py > /w/labkit.log 2>&1 &
  for i in $(seq 1 150); do python3 -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/status')" 2>/dev/null && break; sleep 0.2; done
  start=$(date +%s)
  timeout "$ORACLE_TIMEOUT" python3 /w/solution/solve.py > /w/oracle.log 2>&1; echo "rc=$?" > /w/oracle.rc
  echo "seconds=$(( $(date +%s) - start ))" >> /w/oracle.rc
fi
python3 /w/inner.py > /w/inner.json 2> /w/inner.err
'''

INNER = r'''import copy, importlib.util, json, sys
sys.path.insert(0, '/lab')
import task_lab
state = json.load(open('/w/lab-results/state.json'))
sub = state.get('submission') or {}
out = {'oracle_submitted': bool(sub), 'oracle_status': sub.get('status'), 'oracle_success': bool(sub.get('success')),
       'oracle_report': sub.get('report'), 'charge_used': state['used'], 'budget': state['budget'],
       'evaluations': len(state['trace']), 'failed_evaluations': sum(r['status'] == 'failed' for r in state['trace'])}
if sub:
    trace = [r for r in state['trace'] if r['status'] != 'running']
    validation = [v['result'] for v in state['validation']]
    payload = sub['payload']
    spec = importlib.util.spec_from_file_location('ablations', '/w/solution/ablations.py')
    ab = importlib.util.module_from_spec(spec); spec.loader.exec_module(ab)
    def grade(p, t):
        try:
            ok, report = task_lab.grade(copy.deepcopy(p), copy.deepcopy(t), copy.deepcopy(validation))
            return bool(ok), None
        except Exception as exc:
            return False, '%s: %s' % (type(exc).__name__, str(exc)[:200])
    def rounded(x):
        if isinstance(x, float):
            return float('%.4g' % x)
        if isinstance(x, list):
            return [rounded(v) for v in x]
        if isinstance(x, dict):
            return {k: rounded(v) for k, v in x.items()}
        return x
    rows = []
    for dp, name, mutate in ab.ABLATIONS:
        p2, t2 = mutate(copy.deepcopy(payload), copy.deepcopy(trace))
        ok, err = grade(p2, t2)
        rows.append({'decision_point': dp, 'name': name, 'passes_grading': ok, 'changed': p2 != payload or t2 != trace,
                     'grade_exception': err})
    out['ablations'] = rows
    out['rounding_4sig_passes'] = grade(rounded(payload), trace)[0]
    out['robustness'] = [{'name': name, 'passes_grading': grade(*mutate(copy.deepcopy(payload), copy.deepcopy(trace)))[0]}
                         for name, mutate in getattr(ab, 'ROBUSTNESS', [])]
    first, second = grade(payload, trace), grade(payload, trace)
    out['deterministic'] = first == second and first[0] == bool(sub.get('success'))
print(json.dumps(out))
'''


def leaks(folder):
    hidden = json.loads((folder/'environment/lab/private/hidden_values.json').read_text())
    needles = set()
    for value in hidden:
        if isinstance(value, str) and len(value) >= 4:
            needles.add(value)
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            for text in {f'{value:.4g}', f'{value:.5g}', f'{value:.6g}', repr(value)}:
                if len(text.replace('-', '').replace('.', '')) >= 4:
                    needles.add(text)
    visible = [folder/'instruction.md', folder/'environment/lab/public.json', folder/'tests/public.json',
               *[p for p in (folder/'environment/agent').rglob('*') if p.is_file() and p.suffix in TEXT]]
    found = []
    for path in visible:
        text = path.read_text(errors='ignore')
        found += [f'{path.relative_to(folder)}: {n}' for n in needles if n in text]
    return found, len(hidden)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('task_id')
    p.add_argument('--cpuset', default='2-5')
    p.add_argument('--oracle-timeout', type=int, default=14400)
    p.add_argument('--reuse-state', action='store_true', help='Re-grade the previous oracle state without rerunning it')
    a = p.parse_args()
    folder = TASKS/a.task_id
    issues = problems(folder)
    design = json.loads((folder/'design.json').read_text()) if (folder/'design.json').exists() else {}
    dps = design.get('decision_points', [])
    if len(dps) < 5:
        issues.append(f'only {len(dps)} decision points (need >= 5)')
    rubric_items = re.findall(r'^\s*\d+\.\s+\S', (folder/'rubric.md').read_text(), re.M) if (folder/'rubric.md').exists() else []
    if len(rubric_items) < 4:
        issues.append(f'rubric has {len(rubric_items)} numbered items (need >= 4)')
    sizes = {str(f.relative_to(folder)): f.stat().st_size for f in folder.rglob('*') if f.is_file()}
    big = [k for k, v in sizes.items() if v > MAX_FILE_MB * 2**20]
    total_mb = sum(sizes.values()) / 2**20
    if big and 'size' not in design.get('exceptions', {}):
        issues.append(f'files over {MAX_FILE_MB} MB: {big}')
    if total_mb > MAX_TASK_MB and 'size' not in design.get('exceptions', {}):
        issues.append(f'task is {total_mb:.1f} MB (limit {MAX_TASK_MB})')
    found, n_hidden = leaks(folder)
    if n_hidden == 0:
        issues.append('private/hidden_values.json lists no hidden values')
    if found:
        issues.append(f'hidden values visible to the agent: {found[:10]}')
    work = ROOT/'runs/check'/a.task_id
    image = f'check-{a.task_id}-lab'.lower()
    inner = {}
    if not issues or a.reuse_state:
        build = subprocess.run(['docker', 'build', '-q', '-t', image, '-f', str(folder/'environment/Dockerfile.lab'),
                                str(folder/'environment')], capture_output=True, text=True)
        if build.returncode:
            issues.append('lab image build failed: ' + build.stderr[-500:])
        else:
            if not a.reuse_state:
                if work.exists():   # earlier runs may have left root-owned files
                    subprocess.run(['docker', 'run', '--rm', '-v', f'{work}:/w', image, 'find', '/w', '-mindepth', '1',
                                    '-delete'], capture_output=True)
                    shutil.rmtree(work)
                work.mkdir(parents=True)
                shutil.copytree(folder/'environment/agent', work/'app')
                shutil.copytree(folder/'solution', work/'solution')
            (work/'run.sh').write_text(RUN_SH)
            (work/'inner.py').write_text(INNER)
            started = time.time()
            subprocess.run(['docker', 'run', '--rm', '--network', 'none', '--cpuset-cpus', a.cpuset,
                            '--user', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp',
                            '-e', f'ORACLE_TIMEOUT={a.oracle_timeout}', '-e', f'REUSE={int(a.reuse_state)}',
                            '-v', f'{work}:/w', image, 'bash', '/w/run.sh'], capture_output=True, text=True)
            wall = time.time() - started
            try:
                inner = json.loads((work/'inner.json').read_text())
            except (OSError, json.JSONDecodeError):
                issues.append('inner check produced no result: ' + (work/'inner.err').read_text()[-800:]
                              if (work/'inner.err').exists() else 'inner check did not run')
            inner['container_seconds'] = round(wall)
            if inner and not inner.get('oracle_success'):
                issues.append(f"reference solution not accepted: {inner.get('oracle_report')}; see {work/'oracle.log'}")
            names = {r['name'] for r in inner.get('ablations', [])}
            for dp in dps:
                if dp.get('ablation') not in names:
                    issues.append(f"decision point {dp.get('id')} has no ablation in ABLATIONS")
            for row in inner.get('ablations', []):
                if row['passes_grading'] or not row['changed']:
                    issues.append(f"ablation {row['name']} does not fail grading" if row['passes_grading']
                                  else f"ablation {row['name']} does not change the submission")
            if inner.get('oracle_submitted') and not inner.get('rounding_4sig_passes') and 'rounding' not in design.get('exceptions', {}):
                issues.append('rounding the submission to 4 significant digits fails grading (tolerance too tight or undocumented precision)')
            issues += [f"robustness variation {r['name']} fails" for r in inner.get('robustness', []) if not r['passes_grading']]
            if inner.get('oracle_submitted') and not inner.get('deterministic'):
                issues.append('grading is not deterministic')
    summary = {'passed': not issues, 'issues': issues, 'decision_points': len(dps), 'rubric_items': len(rubric_items),
               'task_mb': round(total_mb, 2), 'hidden_values': n_hidden,
               'oracle': {k: inner.get(k) for k in ('oracle_success', 'charge_used', 'budget', 'evaluations',
                                                    'failed_evaluations', 'container_seconds')},
               'ablations': inner.get('ablations'), 'rounding_4sig_passes': inner.get('rounding_4sig_passes'),
               'robustness': inner.get('robustness')}
    out = ROOT/'benchmarks/construction'/a.task_id/'check.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({'passed': summary['passed'], 'summary': summary,
                               'checked_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}, indent=1) + '\n')
    print(json.dumps(summary, indent=1))
    raise SystemExit(0 if summary['passed'] else 1)


if __name__ == '__main__':
    main()
