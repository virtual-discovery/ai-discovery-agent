#!/bin/bash
# Construction-time DFPT dielectric tensors of wurtzite GaN (PBE, QE 6.7) at the archive's HSE structure
# (a = 3.19483 A, c = 5.18861 A, u = 0.37747 from the relaxed 96-atom perfect supercell; ions not re-relaxed).
# Run inside ai-theorist-benchmark-runtime:20261007 with this folder as the working directory and /pseudo
# holding Ga.pbe-dn-rrkjus_psl.1.0.0.UPF and N.pbe-n-rrkjus_psl.1.0.0.UPF (pslibrary 1.0.0, sha256 in README):
#   docker run --rm --cpuset-cpus 2,3,18,19 --user $(id -u):$(id -g) -v $PWD:/w -v $PSEUDO:/pseudo -w /w \
#       ai-theorist-benchmark-runtime:20261007 bash run_dfpt.sh
set -euo pipefail
export OMPI_MCA_btl_vader_single_copy_mechanism=none OMP_NUM_THREADS=1
NP=${NP:-4}
MPI="mpirun --allow-run-as-root --use-hwthread-cpus -np $NP"
A_BOHR=$(python3 -c "print(3.1948296711/0.529177210903)")
COVERA=$(python3 -c "print(5.1886060722/3.1948296711)")
mkdir -p tmp out
cat > out/scf.in <<IN
&control
  calculation='scf', prefix='gan', outdir='./tmp', pseudo_dir='/pseudo', tprnfor=.true.
/
&system
  ibrav=4, celldm(1)=$A_BOHR, celldm(3)=$COVERA, nat=4, ntyp=2, ecutwfc=45, ecutrho=360
/
&electrons
  conv_thr=1.0d-11, mixing_beta=0.4
/
ATOMIC_SPECIES
Ga 69.723 Ga.pbe-dn-rrkjus_psl.1.0.0.UPF
N  14.007 N.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
Ga 0.3333333333 0.6666666667 0.0000000000
Ga 0.6666666667 0.3333333333 0.5000000000
N  0.3333333333 0.6666666667 0.3774692321
N  0.6666666667 0.3333333333 0.8774692321
K_POINTS automatic
8 8 5 0 0 0
IN
$MPI pw.x -nk $NP -in out/scf.in > out/scf.out
cat > out/ph.in <<IN
phonons at Gamma with dielectric tensor and Born charges
&inputph
  prefix='gan', outdir='./tmp', tr2_ph=1.0d-15, epsil=.true., fildyn='out/gan.dyn'
/
0.0 0.0 0.0
IN
$MPI ph.x -nk $NP -in out/ph.in > out/ph.out
cat > out/dynmat.in <<IN
&input fildyn='out/gan.dyn', asr='crystal', lperm=.true., filout='out/dynmat.modes' /
IN
dynmat.x < out/dynmat.in > out/dynmat.out
rm -rf tmp
echo done
