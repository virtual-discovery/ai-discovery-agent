#!/bin/bash
cd "$(dirname "$0")"; export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 NSUPER={NSUPER}; PY=${PY:-python3}
for s in 01_marginals 02_make_pxc 03_numpy_conditional 06_tables; do echo "== $(date -u) $s"; $PY scripts/$s.py && echo "$s DONE" || { echo "$s FAILED"; exit 1; }; done
echo "REALSPACE FINISHED"
