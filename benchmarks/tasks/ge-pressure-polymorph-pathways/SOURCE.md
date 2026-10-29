# Source and transformations

- **Source:** A. Fantasia, F. Rovaris, O. Abou El Kheir, A. Marzegalli, D. Lanzoni, L. Pessina, P. Xiao, C. Zhou et al., *A NN-Potential for phase transformations in Ge* / "Development of a machine learning interatomic potential for exploring pressure-dependent kinetics of phase transitions in germanium" (2024). Materials Cloud Archive, record 238ad-1bg44, DOI 10.24435/materialscloud:r2-qc.
- **License:** CC BY 4.0 (whole archive).
- **File used:** `allfiles.zip` (329,636,512 bytes, md5 d5faa3ad7e44d8b797ff59d1d7ea66c8). It was read in place; only `training_set/{CD,bSn,ST12,BC8,R8}/config_*/OUTCAR` and `allfiles/Fantasy-frz.pb` were used.

## Construction

The scripts are in `benchmarks/construction/ge-pressure-polymorph-pathways/`:

- `extract.py` parses every OUTCAR into settings and ionic steps:
  - energy(σ→0) and TOTEN;
  - cell, volume and stress (kB, as printed);
  - positions and forces.
- `build_records.py` writes `environment/lab/records.json.gz` (922 records) and `environment/agent/catalog.json`.
  - Values are as printed in the OUTCARs, without renormalization.
  - No POTCAR or other licensed input is staged.
- `dump_weights.py` reads the frozen TensorFlow graph `Fantasy-frz.pb` with a minimal protobuf reader (`pbparse.py`, no TensorFlow). It writes the network weights to `environment/lab/potential/dp_weights.npz` and the hyperparameters to `dp_meta.json`.
- `environment/lab/dp_numpy.py` is a numpy re-implementation of the model, written for this task:
  - the hybrid se_e3 + se_e2_a descriptor, fitting net and atom_ener offset;
  - analytic forces and virial.
- `validate_dp.py` checks the evaluator; results are in `dp_validation.json`.
  - Its forces and virial match finite differences to 1e-9 eV/Å and 1e-7 eV.
  - On 300 random archived configurations, energy RMSE is 6.0 meV/atom (2–4 meV/atom for the five phases) and force RMSE 0.03–0.1 eV/Å. An error in the descriptor would give errors orders of magnitude larger.
- `reference.py` runs the potential through the lab's own code path (`task_lab.run_structure`):
  - fixed-cell ionic relaxations of every near-equilibrium frame;
  - cell + ion relaxations at −2…12 GPa.
- `make_task.py` computes the reference, ablations, robustness variants and calibration.

## What is hidden

The graded values are never served by the lab:

- the CD/β-Sn, β-Sn/ST12 and β-Sn/BC8 equal-enthalpy pressures;
- E(ST12) − E(BC8) and E(R8) − E(BC8) at 0 GPa;
- the 0 GPa order.

They come from Birch–Murnaghan fits to all near-equilibrium archived DFT frames. The ST12, BC8 and R8 frames are corrected by the released potential's internal-relaxation energy, because the archived ISIF = 6 relaxations kept the internal coordinates fixed.

## Scope

This is a retrospective analysis of released data, not new DFT. The internal-relaxation correction uses the authors' potential, not PBE itself (see the construction README for the uncertainty this adds).
