# Source and transformations

- **Paper:** M. Crispino, N. Witt, S. Enzner, T. Gorni, L. de' Medici, D. Di Sante, G. Sangiovanni, *Tunable electronic correlations in 135-kagome metals*, arXiv:2512.22576 (2025).
- **Data:** S. Enzner, G. Sangiovanni, *135 Kagome Wannier Hamiltonians*, Zenodo record 17951035, DOI 10.5281/zenodo.17951035, https://zenodo.org/records/17951035.
- **License:** CC BY 4.0 (all files).

## Files used (md5 from the record)

| File | Bytes | md5 |
|---|---|---|
| CsTiSb.zip | 4,977,439 | 3022b11da82d2a3dd91b57d2fe10d75b |
| CsVSb.zip | 3,869,479 | 35a8b537aaf029def98e427f23785bed |
| CsCrSb.zip | 4,255,787 | ba28567cb06c56faaf9188c06ce5a1ca |
| CsTiBi.zip | 3,919,667 | 1d92f3a6309e431f6debbbf1ed61df8f |
| CsCrBi.zip | 3,879,744 | 95d1e17aaeb1c4b9392c096e27afd76f |
| CsCrBi_9.2.zip | 3,895,438 | 9c4a37b89801300f30eb2bd5d86167d3 |
| DFT.zip | 819,205 | f8b8de95fa310d206ef01bdd27577501 |

Member checksums are in `benchmarks/construction/kagome135-correlation-origin/members.md5`.

## Transformations

The construction scripts are in `benchmarks/construction/kagome135-correlation-origin/`.

- **Staged unchanged** (`/app/archive/<folder>/`): `Fermi_level.dat`, `turn_*.win`, `turn_*.wout`, `turn_*_bands.dat.gnu` of each model folder, and `DFT/scf.in`.
- **Staged trimmed:** `DFT/scf.out` without its per-k eigenvalue listing (redundant with the models; its numbers collided by chance with graded values in the leak scan).
- **Not staged:** the pseudopotentials, `rel.in`/`rel.out` and the Wannier band files `turn_*_band.dat` (redundant with the models; chance numeric collisions).
- **Converted** (`/app/models/<folder>.npz`, and the same arrays in the lab): each `turn_*_hr.dat` as R, Wigner–Seitz degeneracies and H(R) (complex128, eV, archive energy scale). The conversion is lossless (hr values have 6 decimals); Hermiticity H(−R) = H(R)† holds exactly.
- **Lab evaluator** (`environment/lab/kag_ssmf.py`): our numpy implementation of Z2 slave-spin mean field with density-density Kanamori interaction and fully-localized-limit double counting, following the paper's supplement (Eqs. S1–S8). It is not the authors' code.

## What is hidden

Every graded number is derived by construction code from lab-path evaluations and from the models:

- Z̄ and Z_min of CsCr3Sb5 (pristine count 34, U = 6 eV, J = 0.18 U);
- the most renormalized orbital (dxz in the archive's site-local frame);
- the filling share of the Ti→Cr change in ln Z̄ (both counterfactual orders averaged);
- the Hund ratio Z̄(J = 0)/Z̄(J = 0.18 U) of CsCr3Sb5;
- Z̄ of both CsCr3Bi5 models at U = 4.5 eV;
- σ of dxz in CsCr3Bi5 (CsCr3Sb5 lattice) and CsTi3Sb5;
- the pristine electron counts (28, 31, 34, 28, 34, 34).

The service never returns any of these as such; it returns per-orbital Z for single requested (model, electrons, U, J/U).

## Scope

The paper's own Z values are not graded. Our evaluator reproduces the paper's Ti/V/Cr-Sb and CsCr3Bi5 ranges at the paper's nominal counts, but not its overdoped CsTi3Bi5 (see the construction README). The paper's electron counts (27/30/33) omit the Cs 6s electron that the archived DFT places in the 33-orbital manifold; the task defines the pristine count by the archive.
