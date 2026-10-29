# Source and transformations

- **Source:** D. W. Tam, R. Sibille, T. Mazet, B. Roessli, A. Piovano, B. Malaman, D. J. Gawryluk, L. Roduit, Y. Shi, X. Yuan, A. A. Turrini, L. Keller and Y. Soh, *Noncollinear kagome magnet driven by giant antisymmetric exchange*, Nat. Commun. (2026), DOI 10.1038/s41467-026-77000-6.
- **Dataset:** "Noncollinear kagome magnet driven by giant antisymmetric exchange", Materials Cloud Archive record 0cwzy-zz096 (2026.137), DOI 10.24435/materialscloud:xw-6f.
- **License:** CC BY 4.0 (all files of the record).

## Files used

| file | bytes | md5 |
|---|---|---|
| `narrowscans.zip` | 1,248,776 | 1900407dd1f7b413ba3786980eb9103a |
| `widescans.zip` | 10,213,141 | 5c96956ea4dace5ef9ca4ca488808de7 |
| `README.md` | 2,064 | d56203a3e02a037e8ee037b8fd4041d0 |
| `readcombined.py` | 812 | 237aa52e042fe890e98614bc636bc8a6 |
| `SpinW_Co3Sn2S2_DMI.m` | 8,578 | 83a27de47d4381f4933e04ab156bd6ea |

`Co3Sn2S2-DMS-analysis.nb` (7,376,932 bytes, md5 b0e5102993736567538f35cd80a165e4) was inspected but is not staged. It only imports and plots the arrays.

## Construction

The scripts are in `benchmarks/construction/co3sn2s2-dmi-diffuse/`. `fetch.py` refetches the files, and `make_task.py` writes every staged and private file.

- **Staged verbatim:** `README.md`, `readcombined.py` and `SpinW_Co3Sn2S2_DMI.m`, in `environment/agent/archive/`.
- **Narrow scans** (`narrowscans/<run>_Cut2D_Data.npy` and `_BinEdges.npy`): served whole, with their four channels unchanged. The integer channels are written as integers, the normalization is rounded to 8 significant digits (NaN written as null), and the separable bin edges are written as 1-D arrays.
- **Wide scans:** only the combined runs named in the README are used (3869, 3871, 3874 and 3875 for 1.5, 165, 180 and 300 K). They are served as square crops of ±0.48 Å⁻¹ around five reflections: (0 0 3), (1 0 1), (−1 0 2), (1 0 −2) and (1 0 4). The component runs 3870, 3872 and 3873 duplicate pixels of the combined files and are not served. The CSV exports are not served either: the archived `3870_Qys.csv` is a copy of `3870_Qxs.csv`, and `readcombined.py` writes lower bin edges + 0.015, which is the upper edge of a 0.015 Å⁻¹ bin rather than its centre.
- **Check:** `make_task.py` asserts that the annular profile computed from every served crop equals the profile from the full map.

## What is hidden

- **Narrow scans** at 170 K (run 3885), 175 K (3879) and 200 K (3884), and the (1 0 1) region of the 180 K wide scan (run 3874). Their annular profiles grade the submission. They are reduced exactly as the instruction defines them, with Poisson uncertainties from the archived counts.
- **Two properties of the archived Hamiltonian**, computed from the SpinW script: the classical ground-state canting angle and the zone-centre linear-spin-wave gap. The script's model is used as written: R-3m, a = 5.357 Å, c = 13.127 Å, Co S = 1 at (½, 0, ½), J1 = −56 meV and J1dmi on bond 1, and Jc1 = −3 meV on bond 3. The J1dmi matrix is attached to the bond Co(0, ½, ½) → Co(½, ½, ½), and images carry R J Rᵀ, transposed for reversed bonds.
  - The gap was checked against the linearised classical (Landau–Lifshitz) k = 0 frequencies.
  - The model also reproduces the IN8 magnon positions listed in the script to within about 0.02 r.l.u. at 32 meV.

## Scope

This is a retrospective analysis of released data and the released model. The SpinW internal choice of generator bond could not be verified without MATLAB, so the instruction fixes it. The opposite orientation (the bond Co(½, ½, ½) → Co(0, ½, ½)) gives 7.75° and 3.31 meV instead of 8.66° and 4.15 meV. The paper's measured canting (about 17°) and gap (2–3 meV, from earlier work) are not graded.
