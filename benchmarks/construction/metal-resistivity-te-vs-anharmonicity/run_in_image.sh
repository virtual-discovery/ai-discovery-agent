#!/bin/bash
# Run a construction script inside the runtime image on the builder's cpuset, with the repository mounted.
#   bash benchmarks/construction/metal-resistivity-te-vs-anharmonicity/run_in_image.sh make_task.py [args]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${MRES_CACHE:-$HOME/.cache/ai-theorist/metal-resistivity-te-vs-anharmonicity}"
docker run --rm --cpuset-cpus "${CPUSET:-2,3,18,19}" --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
  -v "$REPO:/repo" -v "$CACHE:/data:ro" -w /repo/benchmarks/construction/metal-resistivity-te-vs-anharmonicity \
  ai-theorist-benchmark-runtime:20261007 python3 -B "$@"
