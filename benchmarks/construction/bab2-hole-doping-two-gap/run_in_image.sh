#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/bab2-hole-doping-two-gap/run_in_image.sh make_task.py [args]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${BAB2_CACHE:-$HOME/.cache/ai-theorist/bab2-hole-doping-two-gap}"
docker run --rm --cpuset-cpus "${CPUSET:-4,5,20,21}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -v "$REPO:/repo" -v "$CACHE:/data:ro" -w /repo/benchmarks/construction/bab2-hole-doping-two-gap \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
