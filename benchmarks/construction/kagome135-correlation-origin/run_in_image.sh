#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/kagome135-correlation-origin/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted at /data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${KAG_CACHE:-$HOME/.cache/ai-theorist/kagome135-correlation-origin}"
docker run --rm --cpuset-cpus "${CPUSET:-2,3,18,19}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${THREADS:-1}" -e OPENBLAS_NUM_THREADS="${THREADS:-1}" -e WORKERS="${WORKERS:-4}" \
  -v "$REPO:/repo" -v "$CACHE:/data" -w /repo/benchmarks/construction/kagome135-correlation-origin \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
