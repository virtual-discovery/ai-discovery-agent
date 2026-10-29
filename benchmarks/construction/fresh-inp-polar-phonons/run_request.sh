#!/bin/bash
# usage: run_request.sh <name> '<request json>'
# Runs one lab request through environment/lab/qe_driver.py (the task_lab.run code path) inside the
# runtime image on this builder's cpuset; output in ~/.cache/ai-theorist/fresh-inp-polar-phonons/runs/<name>/.
set -eu
NAME=$1; REQ=$2
CPUSET=${CPUSET:-12,13,28,29}
LAB=$(cd "$(dirname "$0")/../../tasks/fresh-inp-polar-phonons/environment/lab" && pwd)
RUNS=$HOME/.cache/ai-theorist/fresh-inp-polar-phonons/runs
mkdir -p "$RUNS/$NAME"
docker run --rm --cpuset-cpus "$CPUSET" --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -e OMPI_MCA_btl_vader_single_copy_mechanism=none -e PYTHONDONTWRITEBYTECODE=1 \
  -v "$LAB":/lab:ro -v "$RUNS":/runs ai-theorist-benchmark-runtime:20261007 \
  python3 -c "import json,sys; sys.path.insert(0,'/lab'); import task_lab; r=task_lab.run(json.loads(sys.argv[1]), '/runs/$NAME'); open('/runs/$NAME/result.json','w').write(json.dumps(r,indent=1))" "$REQ"
python3 -c "import json;r=json.load(open('$RUNS/$NAME/result.json'));print(json.dumps({k:v for k,v in r.items() if k not in ('cell_vectors_A','force_constants_gamma_Ry_per_bohr2','stress_GPa')}))"
