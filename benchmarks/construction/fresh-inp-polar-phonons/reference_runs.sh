#!/bin/bash
# Reference calculations for fresh-inp-polar-phonons, all through task_lab.run (run_request.sh).
# a0 = 5.9624 A is the Birch-Murnaghan minimum of the 6x6x6-shifted, 50 Ry SCF energies (eos_* runs).
cd "$(dirname "$0")"
R=./run_request.sh
run() { [ -f "$HOME/.cache/ai-theorist/fresh-inp-polar-phonons/runs/$1/result.json" ] || $R "$1" "$2" > /dev/null; }
A0=5.9624; AM=5.9028; AP=6.0220; AMH=5.9326; APH=5.9922; AEXP=5.8687
# EOS (6x6x6 shifted, 50 and 65 Ry)
for a in 5.88 5.92 5.96 6.00 6.04; do run eos_a${a}_k6s1 "{\"kind\":\"scf\",\"lattice_A\":$a,\"kmesh\":6,\"kshift\":1}"; done
for a in 5.88 5.92 5.96 6.00 6.04; do run eos_a${a}_k6s1_e65 "{\"kind\":\"scf\",\"lattice_A\":$a,\"kmesh\":6,\"kshift\":1,\"ecutwfc_Ry\":65}"; done
# equilibrium
run ref_a0_dfpt_k6s1 "{\"kind\":\"dfpt\",\"lattice_A\":$A0,\"kmesh\":6,\"kshift\":1}"
run ref_a0_diel_k8s1 "{\"kind\":\"dielectric\",\"lattice_A\":$A0,\"kmesh\":8,\"kshift\":1}"
run ref_a0_diel_k6s1 "{\"kind\":\"dielectric\",\"lattice_A\":$A0,\"kmesh\":6,\"kshift\":1}"
run ref_a0_diel_k4s1 "{\"kind\":\"dielectric\",\"lattice_A\":$A0,\"kmesh\":4,\"kshift\":1}"
run ref_a0_diel_k10s1 "{\"kind\":\"dielectric\",\"lattice_A\":$A0,\"kmesh\":10,\"kshift\":1}"
run ref_a0_diel_k8s0 "{\"kind\":\"dielectric\",\"lattice_A\":$A0,\"kmesh\":8,\"kshift\":0}"
run ref_a0_dfpt_k4s1 "{\"kind\":\"dfpt\",\"lattice_A\":$A0,\"kmesh\":4,\"kshift\":1}"
run ref_a0_dfpt_k8s0 "{\"kind\":\"dfpt\",\"lattice_A\":$A0,\"kmesh\":8,\"kshift\":0}"
# hydrostatic +-1%
for n in am ap; do
  a=$([ $n = am ] && echo $AM || echo $AP)
  run ref_${n}_dfpt_k6s1 "{\"kind\":\"dfpt\",\"lattice_A\":$a,\"kmesh\":6,\"kshift\":1}"
  run ref_${n}_diel_k8s1 "{\"kind\":\"dielectric\",\"lattice_A\":$a,\"kmesh\":8,\"kshift\":1}"
  run ref_${n}_diel_k4s1 "{\"kind\":\"dielectric\",\"lattice_A\":$a,\"kmesh\":4,\"kshift\":1}"
  run ref_${n}_diel_k8s0 "{\"kind\":\"dielectric\",\"lattice_A\":$a,\"kmesh\":8,\"kshift\":0}"
done
# strained cell
run ref_strain_dfpt_k6s1 "{\"kind\":\"dfpt\",\"lattice_A\":$A0,\"strain_zz\":-0.02,\"kmesh\":6,\"kshift\":1}"
run ref_strain_dfpt_k4s1 "{\"kind\":\"dfpt\",\"lattice_A\":$A0,\"strain_zz\":-0.02,\"kmesh\":4,\"kshift\":1}"
# experimental lattice parameter (wrong reference volume)
run ref_aexp_dfpt_k6s1 "{\"kind\":\"dfpt\",\"lattice_A\":$AEXP,\"kmesh\":6,\"kshift\":1}"
# robustness: denser meshes and smaller steps
for n in am ap; do
  a=$([ $n = am ] && echo $AM || echo $AP)
  run ref_${n}_diel_k10s1 "{\"kind\":\"dielectric\",\"lattice_A\":$a,\"kmesh\":10,\"kshift\":1}"
done
for n in amh aph; do
  a=$([ $n = amh ] && echo $AMH || echo $APH)
  run ref_${n}_dfpt_k6s1 "{\"kind\":\"dfpt\",\"lattice_A\":$a,\"kmesh\":6,\"kshift\":1}"
  run ref_${n}_diel_k8s1 "{\"kind\":\"dielectric\",\"lattice_A\":$a,\"kmesh\":8,\"kshift\":1}"
done
echo REFERENCE_RUNS_DONE
