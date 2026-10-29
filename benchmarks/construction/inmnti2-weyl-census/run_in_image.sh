#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/inmnti2-weyl-census/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted read-only at /data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${INMNTI2_CACHE:-$HOME/.cache/ai-theorist/inmnti2-weyl-census}"
docker run --rm --cpuset-cpus "${CPUSET:-6,7,22,23}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" -e OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS:-4}" \
  -v "$REPO:/repo" -v "$CACHE:/data" -w /repo/benchmarks/construction/inmnti2-weyl-census \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
