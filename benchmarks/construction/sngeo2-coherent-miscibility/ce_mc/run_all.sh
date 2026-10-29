#!/bin/bash
# Cluster expansion + semi-grand-canonical Monte Carlo used to calibrate the incoherent-Tc window.
# Run inside the runtime image with the task cpuset, from this folder, after fetch.py:
#   docker run --rm --cpuset-cpus <cpus> --user $(id -u):$(id -g) -v $HOME/.cache/ai-theorist:$HOME/.cache/ai-theorist \
#       -v $PWD:/w -w /w -e HOME=$HOME -e PYTHONDONTWRITEBYTECODE=1 ai-theorist-benchmark-runtime:20261007 bash run_all.sh
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -B big.py            # 44 orbits (pairs <= 9 A, triplets <= 6.6 A, quadruplets <= 5.6 A) -> big.pkl
python3 -B lassoscan.py      # fit quality vs Lasso lambda
# Binder cumulant of x at the max-variance chemical potential (histogram reweighting), L x L x 1.5L cells
for lam in 0.1 0.3 1.0 3.0; do
  for L in 6 8; do
    python3 -B tc.py $lam $L $((L*3/2)) 2800,2700,2600,2500,2400,2300,2200,2000,1800,1600 1500 20 l${lam}_$L > logs/log_l${lam}_$L.txt &
  done
  wait
done
