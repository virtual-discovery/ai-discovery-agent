#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/hgte-gated-subband-densities/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted at /data; the scratch folder at /scratch; cwd is this folder.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${HGTE_CACHE:-$HOME/.cache/ai-theorist/hgte-gated-subband-densities}"
SCRATCH="${HGTE_SCRATCH:-$CACHE/scratch}"
mkdir -p "$SCRATCH"
docker run --rm --cpuset-cpus "${CPUSET:-14,15,30,31}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" -e OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS:-4}" \
  -v "$REPO:/repo" -v "$CACHE:/data" -v "$SCRATCH:/scratch" -w /repo/benchmarks/construction/hgte-gated-subband-densities \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
