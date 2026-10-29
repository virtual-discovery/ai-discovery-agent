"""Create a benchmark task skeleton with the shared lab files vendored in.

    python -m benchmarks.tools.new_task TASK_ID --start released --catalog-id ID \
        --source-group GROUP --title "Short title"
"""
import argparse
import json
from pathlib import Path
import shutil
from benchmarks.verify import ROOT

TOOLS = ROOT/'benchmarks/tools'
TASKS = ROOT/'benchmarks/tasks'
VENDORED = {'environment/lab/labkit.py': 'labkit/labkit.py',
            'environment/agent/lab_client.py': 'labkit/lab_client.py',
            'tests/verify.py': 'labkit/verify.py'}

STUBS = {
    'instruction.md': '# {title}\n\n<Scientific question and graded fields only; do not prescribe the workflow.>\n',
    'rubric.md': '# Physical-reasoning rubric\n\nEach item is a checkable statement the report must make explicitly.\n\n1. ...\n',
    'SOURCE.md': '# Source and transformations\n\n<Citation, archive URL, license, files used with SHA256, transformations, what is hidden.>\n',
    'environment/lab/task_lab.py': '"""Lab logic for {task_id}: see benchmarks/tools/labkit/labkit.py for the contract."""\n'
                                   'import json\nfrom pathlib import Path\n\nHERE = Path(__file__).parent\n'
                                   "PRIVATE = json.loads((HERE/'private/hidden.json').read_text())\n\n\n"
                                   'def validate(request):\n    raise ValueError("not implemented")\n\n\n'
                                   'def run(request, workdir):\n    raise NotImplementedError\n\n\n'
                                   'def grade(payload, trace, validation):\n    return False, {{}}\n',
    'environment/lab/private/hidden.json': '{{}}\n',
    'environment/lab/private/hidden_values.json': '[]\n',
    'solution/solve.py': '"""Reference solution for {task_id}; runs as the Harbor oracle and in check_task."""\n'
                         'import json\nimport os\nimport subprocess\nfrom pathlib import Path\n\n'
                         "APP = Path(os.environ.get('APP_DIR', '/app'))\n\n\n"
                         'def lab(endpoint, payload=None):\n'
                         "    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])\n"
                         '    out = subprocess.run(cmd, capture_output=True, text=True)\n'
                         '    result = json.loads(out.stdout)\n'
                         "    if 'error' in result:\n        raise RuntimeError(result['error'])\n"
                         '    return result\n',
    'solution/ablations.py': '"""Plausible wrong paths, one or more per decision point; each must fail grading.\n\n'
                             'ABLATIONS lists (decision_point_id, name, mutate) where mutate(payload, trace)\n'
                             'returns a modified (payload, trace) as the wrong path would have produced it.\n"""\n\n'
                             'ABLATIONS = []\n\n# Defensible variations that must still pass (rounding to 4 significant\n'
                             '# digits is checked automatically).\nROBUSTNESS = []\n',
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('task_id')
    p.add_argument('--start', choices=['released', 'fresh'], required=True)
    p.add_argument('--catalog-id', required=True)
    p.add_argument('--source-group', required=True)
    p.add_argument('--title', required=True)
    p.add_argument('--agent-timeout', type=int, default=10800)
    p.add_argument('--cpus', type=int, default=2)
    p.add_argument('--memory-mb', type=int, default=4096)
    a = p.parse_args()
    folder = TASKS/a.task_id
    if folder.exists():
        p.error(f'{folder} exists')
    shutil.copytree(TOOLS/'template', folder, ignore=shutil.ignore_patterns('*.tmpl'))
    fields = dict(task_id=a.task_id, title=a.title, start=a.start, source_group=a.source_group,
                  catalog_id=a.catalog_id, agent_timeout=a.agent_timeout, cpus=a.cpus, memory_mb=a.memory_mb)
    toml = (TOOLS/'template/task.toml.tmpl').read_text()
    for key, value in fields.items():
        toml = toml.replace('{%s}' % key, str(value))
    (folder/'task.toml').write_text(toml)
    for dest, src in VENDORED.items():
        (folder/dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(TOOLS/src, folder/dest)
    for dest, text in STUBS.items():
        (folder/dest).parent.mkdir(parents=True, exist_ok=True)
        (folder/dest).write_text(text.format(**fields))
    (folder/'environment/lab/public.json').write_text(json.dumps(
        {'task_id': a.task_id, 'solver_budget': 40, 'reserved_validation_charge': 0}, indent=1) + '\n')
    (folder/'tests/public.json').write_text(json.dumps(
        {'task_id': a.task_id, 'solver_budget': 40, 'min_report_chars': 400}, indent=1) + '\n')
    (folder/'design.json').write_text(json.dumps({
        'task_id': a.task_id, 'catalog_id': a.catalog_id, 'start': a.start, 'source_group': a.source_group,
        'title': a.title, 'decision_points': [], 'graded_fields': [], 'length_estimate': {},
        'runtime_dependencies': [], 'data_provenance': {}, 'exceptions': {}}, indent=1) + '\n')
    construction = ROOT/'benchmarks/construction'/a.task_id
    construction.mkdir(parents=True, exist_ok=True)
    if not (construction/'README.md').exists():
        (construction/'README.md').write_text(f'# Construction of {a.task_id}\n\nScripts that download the pinned archive and produce the staged and private data.\n')
    print(f'Created {folder.relative_to(ROOT)} and {construction.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
