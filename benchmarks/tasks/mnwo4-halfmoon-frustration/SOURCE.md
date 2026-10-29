# Source and transformations

- **Source:** N. Natsch, T. N. Tošić, J.-R. Soh, N. A. Spaldin, *Paramagnetic half-moon shaped diffuse scattering arising from 3D magnetic frustration*, Materials Cloud Archive (2026), record p8297-g9n39, DOI 10.24435/materialscloud:s5-a8; paper arXiv:2503.23704.
- **License:** CC BY 4.0 (all files of the record).
- **Files used** (record checksums, verified on download by `fetch.py`):
  - `Fig_2.zip` (15,325,551 bytes, md5 aa3841d1f4b7be0718f4cf502ff12a9f): paper Fig. 1b,c run (model inputs `jfile`, `kfile`, `posfile`, `momfile`, `inpsd.dat`; 1.5 K `sqw.loopTemp.out`, `swdos.loopTemp.out`, `totenergy.loopTemp.out`).
  - `Fig_4.zip` (27,103,293 bytes, md5 9520e07c20ddde03d70f42770ec7890c): paper Fig. 3 runs, static S(q) maps at 14, 20, 30 and 60 K (`qfile`, `sq.loopTemp.out`).
  - `Fig_7.zip` (51,334,042 bytes, md5 c2080da562d727f4651f863e0c7a040e): paper Fig. 5b–f runs, 20 K maps of the reduced Hamiltonians.
  - `Fig_5_and_Fig_10.zip` (19,989,965 bytes, md5 eb1da6949f380e604eb8610030a2d924): paper SI Figs. 7/9 runs, 11 K static S(q) on the full grids of 36×24×36 and 48×24×24 supercells.
  - Not used: `Fig_9.zip` (SI Fig. 8b, a 48×24×24 repeat of the 20 K map) and `files_description.md`.

## What is author output and what is forward modelling

- **Author output (served or withheld unchanged except for rounding to 4 significant digits):** every spectrum, map, reduced-model map and zone slice is an archived UppASD output of the authors' model. The experimental neutron data of the paper are not in the archive and are not used; the authors' model is the system under study ("teacher").
- **Forward modelling by the task builders (construction only):** the reference propagation vector (Luttinger–Tisza maximum of J(q)) and the reference spiral–AF1 energy difference are computed from the archived `jfile`/`kfile` with `model.py`. Calibration fits (`calibrate.py`) use linear spin-wave theory and the self-consistent Gaussian approximation; a Monte Carlo check (`mc.c`) confirmed that equilibrium Monte Carlo of the archived Hamiltonian reproduces the archived maps no better than SCGA (block Pearson 0.92/0.87 vs 0.93/0.89 at 14/20 K), i.e. the archived maps are equilibrium-like and SCGA is adequate.

## Transformations

- `extract.py` reads the archived files, asserts that every run uses the same `jfile`, `kfile` and cell (apart from the reduced models), that the archived q grids reproduce the Cartesian q of the outputs, and that the map grids lie on q = (−l/2, k, l). It converts nothing except to arrays.
- The archived UppASD energies confirm the conventions given to the agent: the archived 1.5 K exchange energy (−0.13974 mRy/atom) equals −(1/N) Σ_i Σ_j J_ij ê_i·ê_j over all `jfile` rows (each bond twice) for the archived final configuration, and the anisotropy energy equals K1⟨(ê·n)²⟩ with K1 = −1.41×10⁻² mRy (easy axis). J and D in meV are the archived mRy values × 13.605693.
- `make_task.py` writes:
  - lab data: spectra truncated to E < 5 meV; the k ≥ 0 half of the 20/30/60 K maps (the full plane follows exactly from S(q) = S(−q)); reduced-model maps 5b and 5c; the 11 K zone data reshaped to (h, k, l) cubes;
  - the catalog (structure, bond lists of the archived `jfile` without values, easy axis, AF1 sign pattern of the archived `initialized_AF1`, protocols from the archived `inpsd.dat` files and the paper's Methods);
  - private targets (the J, D and energy values are stored for reference only and are not graded): block averages (7×11 points) of the archived 14 K |S| map; block averages (5×5 points) of the archived 20 K |S| maps of reduced models 5d and 5e; the archived J and D; the derived propagation vector and spiral–AF1 energy difference.

## What is hidden

- The exchange constants J1…J12 and D (archived `jfile`, `kfile`).
- The 14 K map (paper Fig. 3a) and the 20 K maps of reduced models 5d and 5e (paper Fig. 5d, e); the 5f map is withheld and not graded.

## Scope

Retrospective analysis of released simulation data. Graded: the propagation vector, the three withheld archived maps and the origin category. The exchange constants, D and the spiral–AF1 energy difference are requested but not graded, because without literature priors the purchasable data do not determine them uniquely (see the construction README).
