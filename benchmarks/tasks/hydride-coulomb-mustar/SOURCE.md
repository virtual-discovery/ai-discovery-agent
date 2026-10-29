# Source and transformations

- **Source:** K. Bozier and C. J. Pickard, *Enhanced Coulomb pseudopotentials in high-pressure hydrides from first principles* (manuscript, 2026). Dataset: Zenodo record 23235290 (2026-10-08), DOI 10.5281/zenodo.23235290, https://zenodo.org/records/23235290.
- **License:** CC BY 4.0 (whole record).
- **File used:** `mantel_data_files_v1.tar.gz` (892,080,187 bytes, md5 31b5f98dafc279ce627f7a1473d08d3d). It is streamed once; only the 459 members listed with their md5 in `benchmarks/construction/hydride-coulomb-mustar/members.md5` are written (about 100 MB). The 3.2 GB of text `W_ee.dat` files are not used: the archive README states, and the construction checked on three directories, that `W_ee.dat` equals `W_ee_raw.npy` divided by the outer product of `dos_raw.npy`.
- **Directories used:** `INPUTS/hydrides/{Mg2IrH6, CaH6, LaH10, YSbH6}`, `INPUTS/Mg2IrH6_pressure_dependence/Mg2IrH6-{5,10,15,20,25,30}GPa`, `INPUTS/H3S_pressure_dependence/H3S-{140,160,170,180,190,200,210,220,230,240}GPa`. Per directory: `IsoME/rescaled-alpha2F.dat`, `IsoME/<seed>.dos.dat` (or `DOS/<seed>.dos.dat` at H3S 170 GPa), `IsoME/*/Summary.dat` and `Info.txt`, `IsoME/*.jl`, `W_ee_raw.npy`, `dos_raw.npy`, `mu_results.dat`, `<seed>_mantel.log`, `<seed>.scf.in`, `relax/*.relax.in`. `FIGURES/pressure_dependence/*/Tc*.txt` were read for cross-checks only and are not staged.

## Transformations

The construction scripts are in `benchmarks/construction/hydride-coulomb-mustar/`:

- `fetch.sh` streams the tarball, checks its md5 and extracts and verifies the members.
- `extract.py` reads them into `pool.json.gz` (kept in the cache, not in the repository) without changing any number.
- `make_task.py` writes the catalog, the lab data and the private and ablation files:
  - α²F: 4 of the 50 phonon-smearing columns are staged (1, 10, 20 and the column used by the archived full-W run, `ind_smear` in Info.txt).
  - W(ε,ε′) = `W_ee_raw` / (dos ⊗ dos), zero where the Mantel DOS vanishes; the solver and `w_kernel` use the window |ε|, |ε′| ≤ 5.6 eV (the solver's energy window is ±5 eV). For CaH₆ and LaH₁₀ (2000-point grids) only the diagonal is staged.
  - Archived run files (Summary.dat, Info.txt) and `mu_results.dat` are served verbatim.
- `eliashberg.py` is a numpy re-implementation of the IsoME imaginary-axis equations (IsoME.jl, MIT licence, was read for its conventions; no code was copied), linearized at Tc. It reproduces the archived full-W Tc of 15 runs to within 1.6 % and the archived fixed-μ\* runs to within 0.5 %.

## What is hidden

- At Mg₂IrH₆ 15 and 25 GPa and H₃S 180 and 190 GPa: W, `mu_results.dat` and all archived runs. (The archive itself holds no W at H₃S 180 GPa.)
- Graded: the authors' archived full-W IsoME Tc at those four pressures; the constant-DOS μ\* (1 eV cutoff) that reproduces the archived full-W Tc of Mg₂IrH₆ at 20 GPa, computed with the validated lab solver; and N(E_F) of H₃S at 180 GPa per formula unit per spin from its DOS file.
- Not graded: the values quoted in the abstract (μ = 0.62, μ\* ≈ 0.39 and the 40 K suppression at 0 GPa).

## Scope

Retrospective reconstruction from released data; the Eliashberg runs are new solutions with an equivalent solver on archived inputs. The manuscript accompanying the record was released in October 2026 and may not yet be peer reviewed.
