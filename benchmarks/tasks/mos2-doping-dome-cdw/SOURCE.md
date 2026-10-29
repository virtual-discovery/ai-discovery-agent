# Source and transformations

- **Source:** N. Girotto Erhardt, J. Berges, S. Poncé and D. Novko, *Understanding the origin of superconducting dome in electron-doped MoS₂ monolayer*, npj 2D Mater. Appl. 9, 44 (2025); arXiv:2412.02822.
- **Archive:** Materials Cloud record 51yyy-ksf72, DOI 10.24435/materialscloud:z1-aw, https://archive.materialscloud.org/records/51yyy-ksf72.
- **License:** CC BY 4.0 (whole record).

## Files used

| File | Bytes | md5 |
|---|---|---|
| `fig2a-d.tar.gz` | 740,670 | 4cdb4ac5e83aaa4e3067d6304cb5e97e |
| `fig2e.tar.gz` | 8,777,795 | 0514e9cf3bd294ef9f49e5d79912aefb |
| `fig3.tar.gz` | 213,504 | dbf67323c8b26e984707a9ab512f979c |
| `fig4.tar.gz` | 517,963 | ad524c1f940d23e3280ac430a32fe286 |
| `fig5.tar.gz` | 2,127 | 25bf433a60498198d5aec2e36ba287bd |
| `input.tar.gz` | 202,049 | 6f541b537c66c7de1f3f92ef30242a17 |
| `README.md` | 1,076 | 3b22839753c99e088435274842d0449d |

`benchmarks/construction/mos2-doping-dome-cdw/members.md5` lists the 93 extracted members that are read, with their checksums.

### Deliberately not used

- `fig6.tar.gz` (97 MB): the elphmod tight-binding supercell model, the same "model III" used by the excluded Schobert et al. (2024) study.
- `figS11.tar.gz`: substrate calculations.
- The Wannier Hamiltonians (`mos2_hr.dat`, `mos2_wsvec.dat`), the phonon dispersions and the band structures.

No model III data and no quantity of the Schobert study is used or reproduced. This task uses only the DFPT/EPW spectral functions, the Wannier-derived DOS, the DFT potential-energy surfaces, soft modes and relaxations, and the published Tc tables of this paper.

## Construction

The scripts are in `benchmarks/construction/mos2-doping-dome-cdw/`:

- `fetch.sh` downloads the six archives and verifies them.
- `extract.py` copies every number used into `pool.json`, unmodified.
- `eliashberg.py` is an independent numpy implementation of the paper's isotropic Eliashberg equations (Eqs. 1–9), in the full-DOS and constant-DOS forms. It follows the numerical conventions of the authors' `fig2e/tc.py` workflow; the authors' GPL solver ebmb was read for those conventions but is not copied or used. `validate_solver.py` compares it with all 318 archived values in `fig2e/tc_*.dat` (McMillan, Allen–Dynes and four Eliashberg variants, μ* = 0–0.20): the largest deviation is 0.21 % and most are below 0.02 %. The same file runs in the lab.
- `make_task.py` writes the catalog, the lab data, the private targets and the ablation values.

## Transformations

- **Spectra:** served exactly as archived (ω in meV, α²F, cumulative λ). The solver rounds α²F and ω (in eV) to 6 decimals as `tc.py` does before solving.
- **DOS:** `dos.dat` and `dos2.dat` are merged into total, dz²-weighted (K) and remainder (Q) columns.
- **Archived Tc tables:** `fig5/*` served as archived, without the 1×1 rows at 0.15 and 0.17.
- **Units:** doping conversions use the archived lattice constant (celldm(1) = 6.019885 bohr, a = 3.1856 Å). With it, 1 electron per 1×1 cell equals 11.3787 × 10¹⁴ cm⁻², which reproduces the archive's charge-concentration axes.

## What is hidden

- **Withheld 1×1 H calculations:** the spectra, DOS, `tc_*.dat` rows and `fig5` rows at x = 0.15 and 0.17. The `lambda_2d.txt` and `omegalog_2d.txt` rows at these dopings are served.
- **`tc_1x1_K`:** archived `tc_013.dat` column DA (Eliashberg, full DOS, FA, μ* = 0.13): 0.15 and 0.17.
- **`critical_doping_e_per_cell`:** crossing of E(last PES point) − E(first PES point) between `pes_074` and `pes_076`, linearly interpolated and divided by 4 (2×2 cell).
- **`tc_1x1_soc_at_critical_doping_K`:** the archived SOC DA values at 0.18 and 0.19, interpolated linearly to the critical doping.
- **`tc_cdw_high_K`:** constant-DOS Eliashberg Tc of the archived FA spectra `2x2cdw_dop1` and `2x2cdw_dop3`. Their dopings (1.4 and 1.6 electrons per 2×2 cell) follow from matching the Allen–Dynes Tc of each spectrum to `fig5/tcs_CDW_800K.txt` within 0.2 %. The paper's Fig. 4 caption quotes n1 = 3.0 × 10¹⁴ cm⁻² for dop1, which is inconsistent with this match (1.4 per 2×2 cell is 3.98 × 10¹⁴ cm⁻²).

## Scope

This is a retrospective analysis of released data. The only new calculations are Eliashberg solutions on archived (or agent-modified) spectra and DOS, made with a solver validated against the authors' own tables.
