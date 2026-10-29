#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/zrsix-butterfly-amr/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted at /data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${ZRSIX_CACHE:-$HOME/.cache/ai-theorist/zrsix-butterfly-amr}"
docker run --rm --cpuset-cpus "${CPUSET:-12,13,28,29}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${NTHREADS:-1}" -e OPENBLAS_NUM_THREADS="${NTHREADS:-1}" -e ZRSIX_NPROC="${NPROC:-4}" -e ZRSIX_TOP="${TOP:-160}" \
  -v "$REPO:/repo" -v "$CACHE:/data" -w /repo/benchmarks/construction/zrsix-butterfly-amr \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
