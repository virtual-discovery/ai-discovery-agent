#!/bin/bash
# Rebuild the construction data of fresh-al-surface-workfunction from scratch (about 6 h on 2 cores).
#   CPUSET=8,9,24,25 bash run_reference.sh
set -eu
cd "$(dirname "$0")"
python3 -B fetch_pseudo.py --install
python3 -B reference_batch.py
python3 -B qe_batch.py reference_batch.json      # skips runs already present in the cache
python3 -B make_task.py
