# Source and transformations

- **Source:** E. Bousquet, E. Lelièvre-Berna, N. Qureshi, J.-R. Soh, N. A. Spaldin, A. Urru, X. H. Verbeek, S. F. Weber, *On the sign of the linear magnetoelectric coefficient in Cr2O3*, J. Phys.: Condens. Matter (2024), doi 10.1088/1361-648X/ad1a59 (arXiv:2309.02095).
- **Archive:** Materials Cloud Archive record xd06z-0pp29, https://archive.materialscloud.org/records/xd06z-0pp29, DOI 10.24435/materialscloud:ek-fp.
- **License:** CC BY 4.0 (all files).

## Files used

| File | Size (bytes) | md5 | sha256 |
|---|---|---|---|
| `DFT_calculations.tar.gz` | 213,571,029 | dd3a0fa4aa1e7b990d9fad58f5f02ca8 | 31e7515eef7f3594aa46b23fd270ef66283ca72168322d6212fe4ed36f0513bb |
| `README.txt` | 10,304 | 9c6dfac9cdfd791755a6bd86dc4d28af | 811a2732a32a914d7455ade75e9534ad30fcfb812b55479831d5e9bbe7c80a94 |
| `files_description.md` | 369 | cde11635fe27af49edafa52a5bc677e9 | 9828717c1e663d7c83c2d134297fac0677ba6c7a9e05fbf0a11c9ed7a16b64b6 |

## Construction

The scripts are in `benchmarks/construction/cr2o3-magnetoelectric-sign/`:

- `fetch.sh` downloads the three files, checks their md5 sums and unpacks the tarball.
- `extract.py` parses the archived outputs into `pool.json`, without recomputing anything:
  - VASP LCALCPOL runs (POSCAR, OUTCAR ion positions, dipole lines, energy) of all 180 displaced Born-charge folders plus the reference;
  - Elk `elk.in` and the final Moments block, SCF loop count and messages of the INFO.OUT of all 112 induced-moment runs, plus the last lines of the two MPI error logs stored in two of those folders;
  - VASP induced-moment runs: total magnetization line, last per-ion spin (and, when printed, orbital) tables, run status;
  - the two phonopy hdf5 files and the Γ-point block of the three band.yaml files; the structure, INCAR, KPOINTS, FORCE_SETS and bands.conf files; the QE input files.
- It asserts two archive inconsistencies that the task keeps as found: the `Disp_-0_005` OUTCARs belong to +0.05 Å runs (Cr1–O4) and +0.005 Å runs (O5, O6), although the POSCARs in those folders are displaced by −0.005 Å; and the Elk `Eu_3_x/Amplitude_0_03` run stopped after 129 SCF loops with "Convergence targets achieved", away from the linear trend of the other 27 amplitudes.
- `make_task.py` writes the catalog, the lab data, the private targets and the ablation values. `model.py` holds the reference formulas.
- No POTCAR or other licensed input is staged. The OUTCAR pseudopotential title and ZVAL lines are served.
- The archive README is staged in the catalog verbatim except that its arXiv DOI link is shortened to "arXiv:2309.02095" (the DOI string happens to contain a graded value at 4 significant digits).

## What is hidden

Nothing archived is withheld from purchase. Graded values are derived from the archive and never served:

- α⊥ of both domains (reference: Elk induced moments at ±0.005, ±0.01 and ±0.02 Å, Born charges from ±0.02 Å, the archived force-constant eigenvectors and eigenvalues);
- the mode fractions and dominant mode;
- the moment derivatives of Eu_3_x and Eu_4_x in the out-pointing domain;
- the Cr1 Born charge components;
- the domain of each dataset;
- the annealed domain.

## Scope

This is a retrospective analysis of released calculations; no DFT is rerun. The archive holds no outputs for the Quantum ESPRESSO calculations, so the QE comparison of the paper cannot be reproduced. The VASP orbital moments are printed only to 3 decimals and are not used.
