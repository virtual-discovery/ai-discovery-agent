"""Run a list of lab requests through environment/lab/task_lab.run (the agent's code path) inside the runtime image.

    python3 run_batch.py JOBS.json        (called by run_batch.sh inside the container)

JOBS.json is a list of [name, request]. Each result goes to /runs/<name>/result.json; existing results are skipped.
Reference runs may use meshes beyond the lab's offered ones (k64, trigonal k56): run() does not validate, so the
same pw.x input writer and parser are used.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, '/lab')
import task_lab  # noqa: E402

for name, request in json.loads(Path(sys.argv[1]).read_text()):
    folder = Path('/runs')/name
    if (folder/'result.json').exists():
        continue
    try:
        result = task_lab.run(request, folder)
        (folder/'result.json').write_text(json.dumps(result, indent=1))
        print(name, json.dumps({k: result[k] for k in ('free_energy_Ry', 'pressure_GPa', 'irreducible_k_points',
                                                        'symmetry_operations', 'wall_seconds')}), flush=True)
    except Exception as exc:
        print(name, 'FAILED', repr(exc), flush=True)
