#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/cu2wse4-orbital-edelstein/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted at /data; the scratch folder (optional) at /scratch.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${CWS_CACHE:-$HOME/.cache/ai-theorist/cu2wse4-orbital-edelstein}"
docker run --rm --cpuset-cpus "${CPUSET:-6,7,22,23}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${THREADS:-1}" -e OPENBLAS_NUM_THREADS="${THREADS:-1}" -e CWS_DEGEN="${CWS_DEGEN:-5e-3}" \
  -e CWS_CACHE_DIR="${CWS_CACHE_DIR:-/data/cache}" -e WORKERS="${WORKERS:-4}" \
  -v "$REPO:/repo" -v "$CACHE:/data" -w /repo/benchmarks/construction/cu2wse4-orbital-edelstein \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
