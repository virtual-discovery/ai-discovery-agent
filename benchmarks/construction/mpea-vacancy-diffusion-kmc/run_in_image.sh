#!/bin/bash
# Run construction steps inside the benchmark runtime image on the builder's cpuset (default 4,5,20,21).
#   ./run_in_image.sh python3 -B reference.py
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
docker run --rm --cpuset-cpus "${CPUSET:-4,5,20,21}" --user "$(id -u):$(id -g)" -e HOME=/tmp -e PYTHONDONTWRITEBYTECODE=1 \
  -v "$ROOT":"$ROOT" -v "$HOME/.cache/ai-theorist/mpea-vacancy-diffusion-kmc":/data -w "$HERE" \
  ai-theorist-benchmark-runtime:20261007 bash -c "gcc -O2 -o kmc ../../tasks/mpea-vacancy-diffusion-kmc/environment/lab/kmc.c -lm && $*"
