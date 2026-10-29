#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/mnte-altermagnet-ahc-doping/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted at /data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${MNTE_CACHE:-$HOME/.cache/ai-theorist/mnte-altermagnet-ahc-doping}"
docker run --rm --cpuset-cpus "${CPUSET:-6,7,22,23}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
  -v "$REPO:/repo" -v "$CACHE:/data" -w /repo/benchmarks/construction/mnte-altermagnet-ahc-doping \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
