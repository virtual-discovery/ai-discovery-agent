#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset.
# usage: bash run_in_image.sh script.py [args...]   (data cache mounted at /data)
HERE=$(cd "$(dirname "$0")" && pwd)
CACHE=${CACHE:-$HOME/.cache/ai-theorist/hydride-coulomb-mustar}
REPO=$(cd "$HERE/../../.." && pwd)
exec docker run --rm --cpuset-cpus "${CPUSET:-4,5,20,21}" --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -e PYTHONDONTWRITEBYTECODE=1 -e OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} -e OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS:-4} \
  -v "$CACHE:/data" -v "$REPO:/repo" -w /repo/benchmarks/construction/hydride-coulomb-mustar \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
