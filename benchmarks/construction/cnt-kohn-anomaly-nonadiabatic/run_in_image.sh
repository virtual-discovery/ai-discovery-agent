#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository and cache mounted.
#   bash benchmarks/construction/cnt-kohn-anomaly-nonadiabatic/run_in_image.sh make_task.py [args]
# The cache (fetch.py output) is mounted at /cache.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${CNT_CACHE:-$HOME/.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic}"
docker run --rm --cpuset-cpus "${CPUSET:-4,5,20,21}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS=4 -v "$REPO:/repo" -v "$CACHE:/cache" -w /repo/benchmarks/construction/cnt-kohn-anomaly-nonadiabatic \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
