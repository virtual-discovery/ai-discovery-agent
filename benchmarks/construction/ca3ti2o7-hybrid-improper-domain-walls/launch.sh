#!/bin/bash
# usage: launch.sh CPU TAG JSON  -- one relaxation in the runtime image with the lab's NEP code
C=$(cd "$(dirname "$0")" && pwd)
L=$C/../../tasks/ca3ti2o7-hybrid-improper-domain-walls/environment/lab
RAW=${RAW:-$HOME/.cache/ai-theorist/ca3ti2o7-hybrid-improper-domain-walls/raw}
mkdir -p $C/runs
docker run --rm --user $(id -u):$(id -g) --cpuset-cpus $1 -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
  -e PYTHONDONTWRITEBYTECODE=1 -e PYTHONPATH=/lab -v $C:/w -v $(cd $L && pwd):/lab:ro -v $RAW:/raw:ro -w /w \
  ai-theorist-benchmark-runtime:20261007 python3 run_wall.py "$3" > $C/runs/$2.log 2>&1
