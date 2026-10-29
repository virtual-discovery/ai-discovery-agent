#!/bin/bash
# Construction-time check of the isolated-atom cutoff dependence (W 30 vs 45 Ry, Ti 35 vs 60 Ry).
# Inputs are the archived single-atom inputs with only ecutwfc/ecutrho, outdir and pseudo_dir changed.
# Pseudopotentials (GBRV, as named in the archived inputs) are fetched from the GBRV site, md5-checked:
#   ti_pbesol_v1.4.uspp.F.UPF  193cf2c0c3c6a4da613d69453fa7e3e9
#   w_pbesol_v1.2.uspp.F.UPF   78d6efd6ad550bd5be97e9711072ceda
# Run from this folder inside ai-theorist-benchmark-runtime:20261007, e.g.
#   docker run --rm --cpuset-cpus <set> --user $(id -u):$(id -g) -e HOME=/tmp -v $PWD:/w -w /w \
#       ai-theorist-benchmark-runtime:20261007 bash run_cutoff_check.sh
# Wall time on 4 logical CPUs: W about 5-9 min each, Ti 35 Ry about 25 min, Ti 60 Ry about 40 min.
set -eu
export OMPI_MCA_btl_vader_single_copy_mechanism=none
mkdir -p pseudo
for f in ti_pbesol_v1.4.uspp.F.UPF w_pbesol_v1.2.uspp.F.UPF; do
  [ -s pseudo/$f ] || python3 -c "import urllib.request,sys; urllib.request.urlretrieve('https://www.physics.rutgers.edu/gbrv/$f', 'pseudo/$f')"
done
md5sum -c <<'EOF'
193cf2c0c3c6a4da613d69453fa7e3e9  pseudo/ti_pbesol_v1.4.uspp.F.UPF
78d6efd6ad550bd5be97e9711072ceda  pseudo/w_pbesol_v1.2.uspp.F.UPF
EOF
for f in W_30 W_45 Ti_35 Ti_60; do
  mpirun --oversubscribe -np 4 pw.x -in $f.in > $f.out 2>&1 || true
  grep -E "total energy" $f.out | tail -1
done
rm -rf out_* pseudo
