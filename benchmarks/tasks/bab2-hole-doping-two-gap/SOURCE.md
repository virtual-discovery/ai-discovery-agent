# Source and transformations

- **Source:** M. Bercx, S. Poncé, Y. Zhang, G. Trezza, A. Ghorbani Ghezeljehmeidan, L. Bastonero, J. Qiao, F. O. von Rohr, G. Pizzi, E. Chiavazzo and N. Marzari, *Charting the landscape of Bardeen–Cooper–Schrieffer superconductors in experimentally known compounds*, PRX Energy 4, 033012 (2025), DOI 10.1103/sb28-fjc9.
- **Archive:** Materials Cloud Archive record 0kqen-w8h90 (v4, 2025-11-04), DOI 10.24435/materialscloud:x2-ap, https://archive.materialscloud.org/records/0kqen-w8h90.
- **License:** CC BY 4.0 (all files used).

## Files used

| File | Bytes | md5 |
|---|---|---|
| `bab2-analysis.zip` | 28,310,115 | a63a59c6c913b4cab7d5050eb0f1710c |
| `README.md` | 8,819 | 57c5c2ef290197ab4035051e52fcc5d1 |
| `files_description.md` | 860 | 405192de6ce3f5f9b25d2f376ba92694 |

Only 381 members of the zip are used; their md5 sums are in `benchmarks/construction/bab2-hole-doping-two-gap/members.md5`:

- for each `sensitivity/0`–`6` run: `epw/{scf.in, scf.out, iso.in, iso.out, epw_int.in, epw_int.out, bab2.a2f}`, all `bab2.imag_aniso_gap0_*` and `bab2.pade_aniso_gap0_*` files, `ph/{ph.in, scf.in, bab2.dyn0-21}` and `relax/{relax.in, relax.out}`;
- `331/{ph.in, scf_331.in, scf_331.out, bab2.dyn0-2}`;
- `331_cs/{ph.in, scf.in, scf.out, relax.in, relax.out, bab2.dyn0-2}`.

## Construction

The scripts are in `benchmarks/construction/bab2-hole-doping-two-gap/`:

- `fetch.py` downloads the three files through the archive's per-file links, checks them and extracts the 381 members.
- `extract.py` writes `pool.json`. Input files are copied verbatim, selected output lines verbatim, and numbers are parsed from the archived files without change: DFPT frequencies, the α²F table, the iso.out solver section, the anisotropic log lines, and the gap-distribution files. It asserts the folder dopings (from `tot_charge`) against the Fermi energies and electron counts.
- `make_task.py` writes the catalog, the lab data, the private targets and the ablation values. It reads the isotropic Eliashberg references from the lab's own solver (`environment/lab/eliashberg.py`). Converged tightly at the archived settings, that solver reproduces every converged archived EPW isotropic gap to within 0.021 meV (median 0.001 meV).

The lab serves the archived content under the archive's folder names, with the README text quoted verbatim. It adds one computation: new isotropic Migdal–Eliashberg solutions on an archived α²F.

## What is hidden

- **Never served:** the DFPT results of `sensitivity/3` (x = 0.08) and `sensitivity/5` (x = 0.065).
- **Graded values:** quantities derived from the archive, held in `private/hidden.json` and never served directly:
  - **x_c:** the zero of ω² at H, interpolated between the archived x = 0.06 and 0.065 runs;
  - **isotropic Tc:** at x = 0.065 and 0.10 (μ* = 0.10, 500 meV cutoff), from the leading eigenvalue of the linearized equations of the lab solver;
  - **anisotropic Tc:** at x = 0.065, from Δ_max² extrapolated from the archived 55 and 60 K;
  - **gap ratio:** at 25 K, from the cluster means of the archived imaginary-axis distribution;
  - **homogeneous 1/9 soft mode:** ω² at H extrapolated from the archived x = 0.08 and 0.10 runs;
  - **Cs soft mode:** the lowest archived frequency of `331_cs` at its folded H point.

## Scope

This is a retrospective analysis of released calculations; no DFT, DFPT or EPW calculation is redone. The tolerances let the defensible analyses in `solution/ablations.py` pass and make the plausible wrong paths fail.
