# Source and transformations

- **Source:** G. Tenti, B. Jäckl, K. Nakano, M. Rupp, M. Casula, *Hydrogen liquid-liquid transition from first principles and machine learning*, Phys. Rev. B (2025), doi:10.1103/pbrk-3zgd; arXiv:2502.02447.
- **Data:** Zenodo record 14790037, *Additional data for "Hydrogen liquid-liquid transition from first principles and machine learning"*, doi:10.5281/zenodo.14790037, https://zenodo.org/records/14790037.
- **License:** CC BY 4.0 (all files).

## Archive files read

All files were read by HTTP range requests into the remote tar files; no archive was downloaded whole.

- **Classical MD:** 23 tars `N_{128,256,512,2048}_T{...}K.tar` (Zenodo md5 listed in `benchmarks/construction/dense-hydrogen-llt-crossover/README.md`). From every run folder `T<T>K/rs<rs>/`:
  - `input.lammps`, `log.lammps` and every `slurm-*.out` (1491 members; SHA256 in `members_thermo.json`);
  - six contiguous byte ranges of the nested `traj/trajectory.tar` for every run at 900–1100 K (SHA256 of each range in `members_frames.json`).
- **PIMD:** `PIMD.tar` (README and all 175 `thermodata.npy` files inside the nested `PIMD/900K.tar` and `PIMD/1000K.tar`; SHA256 in `members_pimd.json`). Despite the suffix, these are text tables.
- **Not used:** the AIMD tars (`aimd_*K.tar`, about 3 GB each of QE output at a single density) and `mace_mod_loss_models.tar`.

## Transformations

The scripts are in `benchmarks/construction/dense-hydrogen-llt-crossover/`.

- `extract.py` parses every thermo block of each run's outputs and stitches them by step. Long runs were done as chained restart jobs, and `log.lammps` holds only the last job. A later job overrides an earlier one at the same step. The served series keeps every 5th thermo row (every 500 steps = 0.1 ps) of Temp, PotEng, TotEng and Press, unchanged apart from float32 storage of energies relative to a per-run reference that is added back. The job history (start command, step range, errors) is served with the input deck.
- `structure.py` parses the frames found in the trajectory ranges (about ten consecutive 20-fs dumps per range, six ranges per run). From them it computes:
  - max_k S(k)/N (|k| ≤ 5.5 Å⁻¹);
  - the fraction of atoms with a neighbour within 1.05 Å;
  - the stable molecular fraction over 80-fs windows;
  - g(r) per burst;
  - the mean-square displacement between sampled frames.
- PIMD: every 20th row of the columns nmove, H, E_pot, Temp(K), quantum kinetic (virial) and pressure (GPa). All 55 files under the 900K label are byte-identical to files under the 1000K label (their temperature column reads 1000 K). They are served as archived, under their labels.

## What is hidden

Nothing is withheld from purchase. The graded values are derived from the served data (`analysis.py`) and are never returned by the lab:

- the rs of the c_v maximum of the N = 2048 isotherms at 1100, 1200 and 1400 K, from runs with the first 10 ps discarded;
- the classification at 900, 950, 1000 and 1100 K (loop in the N = 2048 isotherm);
- the 900 K low-density phase (crystal: max_k S(k)/N ≈ 0.13 at N = 2048);
- the 1000 K PIMD crossover rs and the NQE crossover-pressure shift relative to classical N = 2048;
- that the 900 K NQE shift cannot be determined from the archive.

## Scope

This is a retrospective analysis of released MACE-PBE simulations. The graded numbers are specific to this model and protocol, not to experiment.
