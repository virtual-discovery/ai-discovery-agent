#!/bin/bash
# usage: run_batch.sh JOBS.json   (JOBS.json inside this folder)
# Runs every request through environment/lab/task_lab.run in the runtime image on this builder's cpuset;
# results in ~/.cache/ai-theorist/fresh-li-bcc-soft-shear/runs/<name>/result.json.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
CPUSET=${CPUSET:-12,13,28,29}
LAB=$(cd "$HERE/../../tasks/fresh-li-bcc-soft-shear/environment/lab" && pwd)
RUNS=$HOME/.cache/ai-theorist/fresh-li-bcc-soft-shear/runs
mkdir -p "$RUNS"
docker run --rm --cpuset-cpus "$CPUSET" --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -e OMPI_MCA_btl_vader_single_copy_mechanism=none -e PYTHONDONTWRITEBYTECODE=1 \
  -v "$LAB":/lab:ro -v "$RUNS":/runs -v "$HERE":/c:ro ai-theorist-benchmark-runtime:20261007 \
  python3 -B /c/run_batch.py "/c/$1"
