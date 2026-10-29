#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/2d-qsh-shc-plateau/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted at /data; the working directory is this folder.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${QSH_CACHE:-$HOME/.cache/ai-theorist/2d-qsh-shc-plateau}"
SCRATCH="${QSH_SCRATCH:-/tmp}"
docker run --rm --cpuset-cpus "${CPUSET:-2,3,18,19}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" -e OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS:-4}" \
  -v "$REPO:/repo" -v "$CACHE:/data" -v "$SCRATCH:/scratch" -w /repo/benchmarks/construction/2d-qsh-shc-plateau \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
