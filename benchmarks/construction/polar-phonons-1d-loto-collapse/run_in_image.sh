#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository and the
# archive cache mounted (cache at /data, read-only).
#   bash benchmarks/construction/polar-phonons-1d-loto-collapse/run_in_image.sh make_task.py
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${PP1D_CACHE:-$HOME/.cache/ai-theorist/polar-phonons-1d-loto-collapse}"
docker run --rm --cpuset-cpus "${CPUSET:-2,3,18,19}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
  -v "$REPO:/repo" -v "$CACHE:/data:ro" -w /repo/benchmarks/construction/polar-phonons-1d-loto-collapse \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
