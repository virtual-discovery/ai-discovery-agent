# Construction of bab2-hole-doping-two-gap

This is a released-start task. It is built from Bercx, Poncé, Zhang, Trezza, Ghorbani Ghezeljehmeidan, Bastonero, Qiao, von Rohr, Pizzi, Chiavazzo and Marzari, *Charting the landscape of Bardeen–Cooper–Schrieffer superconductors in experimentally known compounds*, PRX Energy 4, 033012 (2025). It uses Materials Cloud record 0kqen-w8h90 (DOI 10.24435/materialscloud:x2-ap, CC BY 4.0), file `bab2-analysis.zip`. Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. Run `python3 -B fetch.py [DEST]`. It downloads `bab2-analysis.zip` (28.3 MB, md5 a63a59c6…), `README.md` and `files_description.md` through the archive's per-file links. It then extracts and verifies the 381 members listed in `members.md5`, leaving out the Wannier outputs, ph.out and the supercell force constants. The default DEST is `~/.cache/ai-theorist/bab2-hole-doping-two-gap`.
2. Run `bash run_in_image.sh extract.py /data`. It writes `pool.json` (1.4 MB) with the archived inputs, selected output lines, DFPT frequencies, α²F tables, logs and gap files, and asserts the folder dopings.
3. Run `bash run_in_image.sh make_task.py`. It writes:
   - `catalog.json`;
   - `lab/data.json`;
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json`.

   It asserts the archive facts listed below. `model.py` holds the readers and estimators. The isotropic references come from the lab's own `eliashberg.py`, loaded without writing bytecode.
4. Run `.venv/bin/python -m benchmarks.tools.check_task bab2-hole-doping-two-gap --cpuset <cpus>`. It writes `check.json`.
5. Run `python3 -B record_oracle.py`, then `python3 -B write_design.py`, which writes `design.json` with the calibration table.

`run_in_image.sh` runs a script in `ai-theorist-benchmark-runtime:20261007` on the builder's cpuset (CPUSET, default 4,5,20,21).

## What the task is

The agent buys archived results under the archive's folder names. The catalog quotes the README verbatim, including its doping list (0.0, 0.02, 0.04, 0.06, 0.065, 0.08, 0.1). The runs are:

- `sensitivity/0`–`6`: primitive cell, homogeneous hole doping;
- `331` and `331_cs`: 3×3×1 supercells.

The requests and costs are:

| Request | Cost | Returns |
|---|---|---|
| `dfpt_phonons` | 2 | frequencies at the 21 irreducible q of the 6×6×4 grid |
| `a2f` | 1 | the bab2.a2f file |
| `iso_log` | 1 | iso.in and the solver section of iso.out |
| `aniso_gap` | 2 | one archived gap-distribution file (imaginary axis or Padé real axis) with that temperature's log lines |
| `eliashberg_iso` | 1 | new numpy isotropic Migdal–Eliashberg solutions, up to 8 temperatures |
| `supercell_phonons` | 2 | archived supercell frequencies at Γ and A |

Every primitive-cell purchase carries the run's scf.in, selected scf.out lines and the relaxed cell. The budget is 34 units. The catalog is worth about 300 units plus any number of solver runs. The DFPT results of `sensitivity/3` and `sensitivity/5` are withheld.

The graded fields are:

- `critical_hole_doping`;
- `iso_tc_K` at x = 0.065 and 0.10 (μ* = 0.10, 500 meV cutoff);
- `aniso_tc_K` at x = 0.065;
- `gap_ratio_25K` at x = 0.065;
- `soft_mode_cm1`, for homogeneous 1/9 h/f.u. doping and for the Cs-substituted cell.

## What the archive shows (and the task is built on)

- **Folders are not in README order.** From each run's `tot_charge` (asserted against the Fermi energy and electron count), folders 0–6 hold x = 0.02, 0.04, 0.06, 0.08, 0.10, 0.065 and 0.00. Folders 2 and 5 both print "number of electrons = 15.94"; only scf.in and the Fermi energy separate them.
- **One soft mode at H.** For x ≤ 0.06 a single mode at H = (1/3, 1/3, 1/2) is imaginary: −32.12, −26.22, −19.11 and −7.85 cm⁻¹ at x = 0, 0.02, 0.04 and 0.06. It is stable at +3.45, +15.04 and +22.41 cm⁻¹ at x = 0.065, 0.08 and 0.10. The mode is mainly Ba out-of-plane motion coupled to B in-plane displacements. A 24×24×16 matdyn check confirms that the instability is confined to H, and that x = 0.065 is stable on the dense grid. ω²(x) is nearly linear; its slope falls from 17,200 to 13,800 cm⁻² per h/f.u.
- **λ is inflated by the soft branch.** λ is 2.10 at x = 0 and 1.27 at x = 0.10, while N(E_F) changes by only 5% (0.744 against 0.705 states/spin/eV).
- **The '331' supercell is undoped.** The README describes `331` as "homogeneously 1-hole doped". Its scf input has no tot_charge, and its output shows 144 electrons (9 × 16) and the undoped Fermi energy, 11.3085 eV. Its cell equals the x = 0 primitive cell ×3×3×1, and its folded-H mode (−32.04 cm⁻¹) reproduces the undoped primitive value (−32.12). The Cs cell has 143 electrons (1/9 h/f.u.) and stays unstable at its A point (folded H), at −39.54 and −25.47 cm⁻¹.
- **The archived isotropic runs** use μ* = 0.13, a 0.5 eV cutoff and 15–30 K, which is below Tc for most dopings. The x = 0.08 run diverges to NaN at 30 K, which lies above its Tc (29.6 K at μ* = 0.13), so the gap really has collapsed there. The lab solver reproduces every converged archived gap to within 0.021 meV (median 0.001 meV).
- **The archived anisotropic runs** extend to 60 K, where the gap is still open (Δ_max = 3.94 meV at x = 0.065). The x = 0.065 run starts at 25 K; the other runs start at 15 K. The Padé real-axis files contain negative gaps from 30 K upwards. The imaginary-axis distributions show two clean clusters: 11.63 and 2.34 meV at 25 K for x = 0.065.
- **α²F smearing columns are immaterial.** Near threshold the 10 columns change λ by less than 0.3%, and Tc by less than 0.01 K.

## Deviations from the design (candidates.json / DETAILS.md section 1)

- **What the hidden dopings expose.** Only the DFPT results of x = 0.065 and 0.08 are withheld; their α²F, logs and anisotropic files are served. The design exposed only scalar summaries for those dopings. That version is not uniquely recoverable:
  - Tc and gap ratios are non-monotonic in x near threshold (iso Tc 34.5, 36.7 and 32.4 K at x = 0.06, 0.065 and 0.08);
  - the anisotropic gap at 60 K is 2.75, 3.94 and 2.39 meV at the same three dopings.
- **Targets.**
  - **x_c** is kept. The reference is 0.06419 and the tolerance is ±0.0025, tighter than the design's ±0.004.
  - **Anisotropic Tc at 0.065** is kept. It must be extrapolated beyond 60 K, because no archived temperature closes the gap. The tolerance is ±4% (design ±10%).
  - **Isotropic Tc** is graded at x = 0.065 and 0.10, with μ* = 0.10 at the archived 500 meV cutoff. The design graded x = 0.08 at ω_c = 10 ω_max. Changing the cutoff at fixed μ* moves Tc by only about 1%, so it is not a discriminating fork, and 0.08 is no longer hidden.
  - **Gap ratio** is taken at 25 K (the lowest archived temperature of the 0.065 run), with a tolerance of ±10% (design ±15%).
  - **Added:** the two soft-mode frequencies. They grade the realisability question and the `331` mislabel found during exploration.
- **Decision points.**
  - The design's dp1 (λ inflated below x_c) is a rubric item only. A "highest stable Tc" target was considered and dropped as ambiguous: a global linear fit of ω² puts x_c at 0.0655, above 0.065.
  - The design's dp3 (NaN logs) became part of dp4 and dp5. The archived NaN is above Tc, so reading it as Δ = 0 is correct. The real convergence fork is in the agent's own solver runs (critical slowing down).
  - The design's dp4 (α²F columns) is a rubric item only, since the columns are immaterial.
  - The design's dp5 (μ*/ω_c) became dp4: archived μ* = 0.13 logs, or Allen–Dynes, against new solutions.
  - New dp1 (README order) and dp3 (supercell contents) come from pathologies verified in the archive. dp7 tests the Padé real-axis gaps.

## Decision points (7) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | Folder dopings from tot_charge, not the README order | metadata inconsistency, found in headers |
| dp2 | x_c from ω² (analytic) rather than ω, with the bracketing runs withheld | regime / model choice |
| dp3 | `331` is undoped despite the README (144 electrons); the folded H point is supercell A, not Γ | mislabelled run, found in raw outputs |
| dp4 | Archived μ* = 0.13 logs (15–30 K) or Allen–Dynes, against new Eliashberg solutions at the stated μ* | which calculation answers |
| dp5 | Unconverged gaps above Tc (critical slowing down) | data-quality judgement on own runs |
| dp6 | Anisotropic Tc beyond the archived 60 K: Δ² extrapolation, not Δ-linear, and not the isotropic value | regime / extrapolation |
| dp7 | Gap structure from the imaginary-axis, not the Padé, distribution | data quality |

dp1 cascades into x_c and both isotropic Tc values; dp2 cascades into the homogeneous soft mode.

## Calibration

The full table is in `design.json` under `calibration`.

| Field | Reference | Oracle | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|---|
| x_c (h/f.u.) | 0.06419 | 0.06418 | 0.06437 (ω² piecewise, released), 0.06402 (ω² quadratic), 0.06551 (ω² linear regression) | 0.07038 (ω linear), 0.04273 (README order) | ±0.0025 |
| iso Tc x = 0.065 (K) | 36.71 | 36.71 | 36.71 (converged Δ² extrapolation), 36.71 (column 10) | 34.68 (archived μ* = 0.13 log, −5.5%), 40 (unconverged scan, +9.0%), 32.09 (Allen–Dynes), 33.48 (README) | ±3% |
| iso Tc x = 0.10 (K) | 33.48 | 33.49 | 33.48, 33.48 | 30.58 (archived log, −8.7%), 36 (unconverged, +7.5%), 30.38 (AD), 41.79 (README) | ±3% |
| aniso Tc x = 0.065 (K) | 62.62 | 62.62 | 62.67 (large cluster²), 62.43 (√−F), 61.95 (3-point Δ² fit), 62.28 (Padé max²) | 67.09 (Δ linear, +7.1%), 36.71 (isotropic) | ±4% |
| gap ratio 25 K | 4.966 | 4.966 | 5.107 (medians), 4.574 (unweighted bins, −7.9%) | 5.910 (Padé, +19%) | ±10% |
| homogeneous 1/9 (cm⁻¹) | 25.60 | 25.46 | 25.66, 25.40, 26.40 | 30.81 (ω linear), −32.04 (`331` read as doped) | ±3.5 |
| Cs substituted (cm⁻¹) | −39.54 | −39.54 | none | −0.79 (Γ point) | ±1.5 |

Two values are not ablations:

- **Memorised paper values:** x_c = 0.065 and Tc_aniso = 61.6 K pass; Tc_iso = 32.0 K (μ* = 0.13) and the gap ratio 11.54/1.50 = 7.7 fail.
- **"60 K"** (the last archived temperature, a lower bound) fails by 0.12 K.

## Oracle

The oracle makes 21 evaluations for 30 of 34 units (88%):

- 7 `iso_log` purchases, to map the folders;
- 3 `dfpt_phonons` purchases (x = 0.04, 0.06 and 0.10);
- 5 solver runs: one validation against the archived μ* = 0.13 gaps, then a coarse and a fine scan at each doping;
- 4 imaginary-axis gap files (x = 0.065 at 25, 50, 55 and 60 K);
- both supercells.

It reads the doping from each header, fits ω² quadratically, and uses only converged solver gaps, extrapolated in Δ². It checks the electron count of `331` before using it.

Results:

- **Gate (2026-10-19, cpuset 4,5,20,21):** passed. All 10 ablations fail, the 11 robustness variants pass, and rounding to 4 significant digits passes. There are 7 decision points and 8 rubric items; the package is 1.8 MB with 7 hidden values. The oracle takes 4 s in the gate container.
- **Harbor oracle (`runs/build-bab2-hole-doping-two-gap/oracle-r1-bab2-hole-doping-two-gap`):** reward 1.0, runtime 47 s, 30 of 34 units charged.

## Open risks

- **Dev-easy risk.**
  - Every purchase prints `tot_charge`, so an agent that reads headers defuses dp1 at once.
  - dp2 (ω² near a soft-mode instability) is textbook.
  - The strongest discriminators should be dp3 (an archived run that contradicts its README description), dp5 (unconverged solver output) and dp4 (reusing the archived μ* = 0.13 runs).
- **Memorisation.** The paper's x_c = 0.065 and Tc_aniso = 61.6 K fall inside their tolerances. The anisotropic Tc check mainly rejects the isotropic value and a Δ-linear extrapolation.
- **Thin margins.**
  - The homogeneous soft mode read straight from x = 0.10 without extrapolation (22.41) passes, 0.31 cm⁻¹ inside the tolerance.
  - A global linear fit of ω (the wrong variable) gives x_c = 0.0669 and fails by only 0.0002; it is not listed as an ablation.
  - The archived-log ablation at x = 0.065 fails by 2.5 percentage points.
- **Solver-defined target.** The isotropic Tc is defined by the lab solver, including EPW's truncated Matsubara sum for Z. A different discretisation at the same μ* would differ by about 1%.
- **Reference uses withheld data.** The x_c reference uses the withheld x = 0.065 DFPT. From the purchasable runs, defensible fits give 0.0640–0.0655, all inside the tolerance.

## Disk

- **Peak:** about 325 MB: 236 MB in the cache (the 28 MB zip and the full extraction made during exploration) and 88 MB of matdyn checks in the scratchpad. All were deleted after the Harbor run, together with the gate's lab image.
- **Kept:** the refetch script with checksums and `pool.json` (1.4 MB).
- **Task package:** 1.8 MB.
