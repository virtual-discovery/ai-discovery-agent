#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/halfheusler-strain-hinge/run_in_image.sh script.py [args]
# The archive cache (fetch.py DEST) is mounted read-only at /data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${HH_CACHE:-$HOME/.cache/ai-theorist/halfheusler-strain-hinge}"
docker run --rm --cpuset-cpus "${CPUSET:-4,5,20,21}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" -e OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS:-4}" \
  -v "$REPO:/repo" -v "$CACHE:/data:ro" ${HH_SCRATCH:+-v "$HH_SCRATCH:/scratch"} -w /repo/benchmarks/construction/halfheusler-strain-hinge \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
