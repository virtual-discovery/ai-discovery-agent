"""Benchmark task registry: task rows for the runner, verification and freezing.

    python -m benchmarks.tools.release verify [--development]
    python -m benchmarks.tools.release freeze TASK_ID [TASK_ID ...]
"""
import argparse
import hashlib
import json
import tomllib
from benchmarks.verify import ROOT

TASKS = ROOT/'benchmarks/tasks'
MANIFEST = ROOT/'evidence/release.json'
VENDORED = {'environment/lab/labkit.py': 'labkit/labkit.py',
            'environment/agent/lab_client.py': 'labkit/lab_client.py',
            'tests/verify.py': 'labkit/verify.py'}
REQUIRED = ['task.toml', 'instruction.md', 'rubric.md', 'SOURCE.md', 'design.json', 'environment/docker-compose.yaml',
            'environment/Dockerfile.agent', 'environment/Dockerfile.lab', 'environment/lab/task_lab.py',
            'environment/lab/public.json', 'solution/solve.sh', 'solution/solve.py', 'solution/ablations.py',
            'tests/test.sh', 'tests/public.json', *VENDORED]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest():
    return json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {'benchmark': 'ai-theorist', 'tasks': {}}


def checksum(folder):
    from harbor.models.task.task import Task
    return Task(folder).checksum


def task_rows(require_frozen=True):
    frozen = manifest()['tasks']
    rows = []
    for folder in sorted(p for p in TASKS.iterdir() if (p/'task.toml').exists()):
        meta = tomllib.loads((folder/'task.toml').read_text())['metadata']
        row = {'task_id': folder.name, 'table': 1, 'start': meta['start'], 'study': meta['source_study'],
               'catalog_id': meta.get('catalog_id'), 'frozen': folder.name in frozen}
        if require_frozen and not row['frozen']:
            continue
        rows.append(row)
    return rows


def problems(folder):
    """Structural problems that make a task unusable; empty when sound."""
    issues = [f'missing {name}' for name in REQUIRED if not (folder/name).exists()]
    for dest, src in VENDORED.items():
        if (folder/dest).exists() and sha256(folder/dest) != sha256(ROOT/'benchmarks/tools'/src):
            issues.append(f'{dest} differs from the shared copy')
    if (folder/'environment/lab/public.json').exists() and (folder/'tests/public.json').exists():
        lab = json.loads((folder/'environment/lab/public.json').read_text())
        tests = json.loads((folder/'tests/public.json').read_text())
        if lab['task_id'] != folder.name or tests['task_id'] != folder.name or lab['solver_budget'] != tests['solver_budget']:
            issues.append('task_id/solver_budget inconsistent between lab and tests public.json')
    if any(folder.rglob('__pycache__')):
        issues.append('__pycache__ present (changes the checksum)')
    return issues


def verify_release(require_frozen=True, only=None):
    """Check every selected task; `only` limits a development run to the tasks it will use."""
    frozen = manifest()['tasks']
    bad = {}
    rows = [r for r in task_rows(require_frozen) if only is None or r['task_id'] in only]
    for row in rows:
        folder = TASKS/row['task_id']
        issues = problems(folder)
        if row['frozen'] and checksum(folder) != frozen[row['task_id']]['checksum']:
            issues.append('checksum differs from evidence/release.json')
        if issues:
            bad[row['task_id']] = issues
    if bad:
        raise SystemExit('Benchmark problems: ' + json.dumps(bad, indent=1))
    print(f'Benchmark: {len(rows)} tasks verified ({len(frozen)} frozen)')


def freeze(task_ids):
    data = manifest()
    for task_id in task_ids:
        folder = TASKS/task_id
        issues = problems(folder)
        check = ROOT/'benchmarks/construction'/task_id/'check.json'
        if not check.exists() or not json.loads(check.read_text()).get('passed'):
            issues.append('no passing check.json from check_task')
        if issues:
            raise SystemExit(f'{task_id}: {issues}')
        meta = tomllib.loads((folder/'task.toml').read_text())['metadata']
        data['tasks'][task_id] = {'checksum': checksum(folder), 'start': meta['start'], 'catalog_id': meta.get('catalog_id'),
                                  'source_study': meta['source_study'], 'check': json.loads(check.read_text())['summary']}
    MANIFEST.write_text(json.dumps(data, indent=1) + '\n')
    print(f'Frozen {len(task_ids)}; manifest has {len(data["tasks"])} tasks')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    v = sub.add_parser('verify')
    v.add_argument('--development', action='store_true')
    f = sub.add_parser('freeze')
    f.add_argument('task_ids', nargs='+')
    a = p.parse_args()
    if a.command == 'verify':
        verify_release(require_frozen=not a.development)
    else:
        freeze(a.task_ids)


if __name__ == '__main__':
    main()
