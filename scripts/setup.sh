#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -c 'import sys; assert sys.version_info >= (3,12), "Python 3.12+ is required"'
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m benchmarks.verify
docker info >/dev/null
docker build -t ai-theorist-benchmark-runtime:20261007 -f benchmarks/runtime/Dockerfile benchmarks/runtime
echo 'Setup complete. Activate .venv, then run: python -m benchmarks.run --policy oracle --suite smoke --execute'
