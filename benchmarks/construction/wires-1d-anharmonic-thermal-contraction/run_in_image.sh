#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository and cache mounted.
#   bash benchmarks/construction/wires-1d-anharmonic-thermal-contraction/run_in_image.sh make_task.py /cache/pool.json
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${WIRES_CACHE:-$HOME/.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction}"
docker run --rm --cpuset-cpus "${CPUSET:-2,3,18,19}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
  -v "$REPO:/repo" -v "$CACHE:/cache:ro" -w /repo/benchmarks/construction/wires-1d-anharmonic-thermal-contraction \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
