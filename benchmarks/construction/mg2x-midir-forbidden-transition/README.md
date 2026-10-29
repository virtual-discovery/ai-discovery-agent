# Construction of mg2x-midir-forbidden-transition

This is a released-start task. It is built from H. Sugawara, *A dipole-forbidden conduction-band transition and the mid-infrared absorption of n-type Mg2Si, Mg2Ge and Mg2Sn: all-electron quasiparticle self-consistent GW calculations* (to be submitted to PRB, 2026), using Zenodo record 23319151 (DOI 10.5281/zenodo.23319151, released 2026-10-13). Licenses are mixed per file; see the task's `SOURCE.md`. No AGPL ecalj or Questaal input is used.

## Rebuild

1. Run `python3 -B fetch.py [DEST]`. It downloads the 77.3 MB archive zip and its README, LICENSES and MD5SUMS through Zenodo's per-file links and checks their md5. It then extracts the 776 members of `members.md5` into `DEST/x` and checks each one; the members' checksums equal those in the archive's own `MD5SUMS.txt`. DEST defaults to `~/.cache/ai-theorist/mg2x-midir-forbidden-transition`, which `model.py` reads.
2. Run `python3 -B make_task.py`. It writes:
   - `environment/agent/catalog.json` and `environment/agent/experiment.json`;
   - `environment/lab/data.json` (0.6 MB);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json`.

   It asserts the archive facts the task rests on.
3. Run `.venv/bin/python -m benchmarks.tools.check_task mg2x-midir-forbidden-transition --cpuset 6,7,22,23`, which writes `check.json`.
4. Run `python3 -B write_design.py`, which writes `design.json` with the decision points and the calibration table.

`model.py` holds the readers and a numpy port of the archive's Julia estimators (`qbox.jl`, `ir_raman_ratio.jl`, `effective_masses.jl`). The port reproduces `analysis/ir_raman_ratio.json` to 1e-9 and the masses of `effective-masses.json`.

## What the task is

The agent buys archived results:

- ecalj LDA and QSGW bands, fine scans at X, and QSGW boxes around X (5 units each);
- TOTE files of the G0W0, QSGW, QSGW80 and 100 K-lattice runs;
- Questaal QSGW and QSGW+ladder eigenvalues;
- BerkeleyGW pseudopotential G0W0;
- the PBE momentum matrix elements of Mg2Si on the box grid;
- AHC and ZG electron–phonon renormalizations.

The budget is 32 units, about 28% of a catalogue worth about 114. The digitized spectra are staged free (CC0).

The 16 graded values are:

- the signed X3 − X1 splitting at the LDA, G0W0 and QSGW levels for the three compounds;
- the ladder correction for Mg2Si;
- the infrared/Raman ratio Q = ⟨q⊥²⟩(300 K)/⟨q⊥²⟩(77 K) at fixed n, for Mg2Si at 1e16 and 3e18 cm⁻³ and for Mg2Sn at 1e16 cm⁻³;
- the predicted splitting of Mg2Si at T → 0 and at 300 K, and of Mg2Ge at 300 K.

## What the archive shows (and the task is built on)

- **Mg2Sn is inverted at every level.** X3 lies below X1: X3 − X1 = −0.360 (LDA), −0.165 (G0W0) and −0.131 eV (QSGW). The two states are told apart by their dispersion: X1 has m_l ≈ 0.54–0.59 and the heavier m_t, X3 has m_l ≈ 1.03–1.09. The LDA of Mg2Sn is a semimetal.
- **The TOTE column label is garbled.** The README describes the columns as "eQP (Z = 1), eQP without Z, Z". But column 7 equals ε + Z(column 8 − ε) in every row, so it is the linearized, Z-renormalized solution; the paper's G0W0 uses it. Column 8 is the Z = 1 value: 0.531 / 0.660 / −0.101 eV against 0.439 / 0.579 / −0.165. In the QSGW iteration files, column 6 is the QSGW level, and columns 7–8 are not (834 meV for Mg2Si).
- **The ladder correction is about −2 meV.** Within Questaal (production basis, window 4/4), the ladder diagrams change the Mg2Si splitting by −2.0 meV, or −2.8, −4.5 and −0.3 meV for the larger windows and the larger basis, but lower the gap by 90–125 meV. QSGW80 lowers the splitting by 62.7 meV. Questaal and ecalj QSGW differ by 28.4 meV (basis).
- **Selection rule.** The PBE |p_13,14|² vanishes at X (about 1e-21) and along the valley axis, grows as q⊥², and does not depend on q_par at leading order.
- **Non-parabolicity.** The lower state flattens along the transverse [011] diagonal (Mg2Si: 155 meV at (0, 0.08, 0.08), where |q⊥| = 0.113, against 264 meV at (0, 0, 0.12) on the axis; Mg2Sn: 106 against 303 meV). The box Q is therefore 4.47 for Mg2Si against 3.90 for a parabolic valley or 3.96 for the X–W line. For Mg2Sn it is 6.50.
- **Temperature.** The production QSGW uses the room-temperature lattice. The 100 K lattice (6.3218 Å) raises the Mg2Si splitting by 23.0 meV. The AHC change of X3 − X1 (24³, η = 10 meV) is +11.9 meV at 0 K and +26.7 meV at 300 K for Mg2Si, and +20.0 meV at 300 K for Mg2Ge.

## Deviations from the design (candidates.json / DETAILS.md section 33)

- **Staging.** The design staged the LDA and QSGW path bands, which would show the QSGW splitting, a graded target, directly. All levels are now purchasable. The splittings are graded with their sign at three levels, so the Mg2Sn inversion (design dp1) and the TOTE columns are tested.
- **Targets.**
  - The design graded "QSGW splitting of Mg2Ge and Mg2Sn"; this is now the full signed 3 × 3 table.
  - The design's "predicted absorption peak of n-Mg2Ge" was dropped. On the box the band peak sits within 1–4 meV of the splitting at X (`lineshape.json`), so it duplicates the splitting, and the peak definition (raw or smoothed histogram) moves it by 3 meV.
  - The "300 K / 77 K strength ratio" (design tolerance ±15%) became Q at two densities for Mg2Si plus Mg2Sn, at ±4%. The design's ±15% would have accepted the parabolic 3.90.
  - The ladder correction and the temperature-consistent predicted splittings were added. They implement the design's dp3 ("level of theory, QSGW overestimation / 0.8 scaling") and dp6 ("zero-point renormalization"), and the coordinator's "temperature or strain shifts".
  - The X1/X3 assignment is graded through the signs.
- **Decision points.** The design's dp5 (experimental peak reading) is in the rubric only: the reading is estimator dependent (maximum 0.559 against half-maximum midpoint 0.564 eV for Mg2Ge), and its spread is not a property of the calculations. dp4 (non-parabolicity) and dp2 (selection rule) are kept, the latter in its "lowest non-vanishing order" form. Three were added:
  - degeneracy and box counting;
  - Mg2Sn transfer;
  - lattice-temperature consistency.
- **Mg2Sn Q is not an archive number.** It is the archive's own estimator (FD over the coarse box) applied to the archived Mg2Sn box. The Mg2Si Q values are the archive's (`ir_raman_ratio.json`).
- **Predicted splittings are a stated combination.** QSGW plus the AHC renormalization, with the 100 K-lattice QSGW for T → 0. They are a combination of archived numbers under conventions given in the instruction, not values printed in the archive.
- **Not staged:** the Julia programs' outputs that would leak targets (`ir_raman_ratio.json`, `lineshape.json`, `absorption_window.json`), the β overlaps (binary), the GaP runs, the Questaal debugging runs and the ecalj numerical-test runs (gw888, pw3–6). The last were left out so that "production settings" is unambiguous.

## Decision points (9) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | X1/X3 identity by dispersion versus energy order (Mg2Sn inverted) | evidence in the dispersions; a cascade into the G0W0 labels through the LDA band numbering |
| dp2 | TOTE columns: garbled README label for G0W0; the QSGW files' QP columns are not QSGW levels | metadata inconsistency found in the raw rows |
| dp3 | Ladder correction within one code and basis, versus across codes or QSGW80 | a discriminating purchase (2 units) |
| dp4 | Lowest-order matrix element q⊥² (forbidden at X), versus allowed, \|q\|² or full PBE | evidence in the purchased matrix elements |
| dp5 | Box over the whole transverse plane, versus a parabolic valley or the X–W line | 5-unit purchase; non-parabolicity visible only off the axes |
| dp6 | Fermi–Dirac at 3e18, with 2 × 3 × 8 counting and half weights on the symmetry planes | regime judgement plus integration bookkeeping |
| dp7 | Mg2Sn needs its own box (lower state X3, strong non-parabolicity) | transfer; a second 5-unit purchase |
| dp8 | T → 0 needs the 100 K-lattice QSGW run | metadata (lattice constants of the runs) |
| dp9 | Electron–phonon renormalization of X3 − X1: zero point, 300 K, sign | combining runs |

## Calibration

The full table is in `design.json` under `calibration` (source: `calibration.json`).

| Field | Reference | Defensible variants | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| X3 − X1, all 9 | see above | path/scan/box/TOTE agree < 0.5 meV | G0W0 from PBE start +20.1/+18.3 meV; QSGW80 −62.7/−53.2/−45.8 meV; Mg2Sn sign; G0W0 Z = 1 column +92/+80/+64 meV | ±0.010 eV |
| Ladder correction (meV) | −2.00 | −2.77 (Nc 8), −4.55 (Nc 12), −0.28 (basis 2) | +26.4 (cross-code), −62.7 (QSGW80), −90.5 (gap change) | ±8 meV |
| Q Mg2Si 1e16 | 4.469 | 4.496 (coarse box at both T), 4.479 (nested at both), 4.497 (unfolded sum) | 4.166 (\|q\|²), 3.964 (X–W line), 3.892 (parabolic), 3.385 (full PBE \|p\|²), 5.468 (no plane weights) | ±4% |
| Q Mg2Si 3e18 | 3.435 | 3.430, 3.442, 3.430 | 3.207 (\|q\|²), 2.990 (parabolic FD), 4.474 (Boltzmann), 2.873 (no spin), 1.785 (one octant) | ±4% |
| Q Mg2Sn 1e16 | 6.502 | 6.512 (Boltzmann), 6.521 (unfolded sum) | 4.469 (Mg2Si value), 4.408 (\|q\|²), 3.892 (parabolic), 3.980 (X–W line), 3.608 (X1 occupied) | ±4% |
| Predicted Mg2Si T0 | 0.4913 eV | 0.4932 (ZG 4³), 0.4929 (adiabatic), 0.4933 (16³) | 0.4794 (no ZPR, −11.9 meV), 0.4683 (room-temperature lattice), 0.4676 (sign reversed) | ±0.008 eV |
| Predicted Mg2Si 300 K | 0.4831 eV | 0.4891 (ZG 4³, +6.0 meV), 0.4841, 0.4845 | 0.4683 (zero point only, −14.9), 0.4564 (no e-ph), 0.5062 (100 K lattice) | ±0.008 eV |
| Predicted Mg2Ge 300 K | 0.5964 eV | 0.5964 (ZG), 0.5960 (adiabatic), 0.5934 (Ge 3d valence, η = 40 meV) | 0.5860 (zero point only, −10.3 meV), 0.5763 (no e-ph) | ±0.008 eV |

Rounding to 4 significant digits passes (checked by the gate).

## Oracle

The oracle makes 19 evaluations for 27 of 32 units (84%):

- six X–Γ scans (LDA, QSGW × 3 compounds), which give the levels at X and identify X1/X3 by the longitudinal mass;
- three G0W0 TOTE files, read in the column checked to equal ε + Z(c8 − ε);
- Questaal prod/qsgw and prod/qsgw_ladder;
- one momentum_matrix request (10 points), checking that p vanishes at X and along q_par and goes as q⊥²;
- the coarse QSGW boxes of Mg2Si and Mg2Sn (Fermi–Dirac, 2 × 3 × 8, trapezoid weights);
- AHC for Mg2Si at 0 and 300 K and Mg2Ge at 300 K;
- the production and 100 K-lattice QSGW TOTE files of Mg2Si (the first checks the column).

The cheapest route to all targets costs about 24 units. The estimated agent length is 45–75 tool calls.

## Results

- **Gate (2026-10-25, cpuset 6,7,22,23):** passed. All 24 ablations fail, the 8 robustness variants pass, rounding to 4 significant digits passes, and grading is deterministic. The oracle took 3 s in the gate container. There are 9 decision points and 8 rubric items; the package is 0.8 MB with 16 hidden values.
- **Harbor oracle (`runs/build-mg2x-midir-forbidden-transition/oracle-r1-mg2x-midir-forbidden-transition`):** reward 1.0, runtime 47 s, charge 27 of 32.

## Open risks

- **Thin margins.**
  - The no-ZPR Mg2Si T0 path fails by 3.9 meV.
  - The zero-point-only Mg2Ge 300 K path fails by 2.3 meV.
  - The ZG 4³ Mg2Si 300 K variant passes by 2.0 meV.
  - The \|q\|² paths fail by 2.7–2.9% beyond the edge of the ±4% window.
- **Physics caveat for Mg2Sn Q.** With a splitting of only 131 meV, the lowest-order matrix element is a poor approximation at 300 K (in PBE, \|p\|² saturates where the states mix). The graded value follows the stated convention (lowest non-vanishing order), and the rubric asks for the caveat. The Mg2Sn box at 300 K holds 0.3% of the weight on its outer faces, and truncating it to q_perp ≤ 0.14 moves Q by 1.1%.
- **Interpretation of the predicted splittings.** They add an LDA electron–phonon renormalization to the QSGW splitting and attribute the lattice effect to the 100 K-lattice run, as the archive does ("remove the temperature inconsistency of using the room-temperature lattice constant in a T = 0 calculation"). The paper may combine them differently. The QSGW + electron–phonon value overestimates the measured Mg2Si band by about 80 meV; this is a fact of the archive, not a grading issue.
- **Dev-easy risk.**
  - dp3, dp4 and dp9 are discoverable by careful reading.
  - The hard parts are dp5–dp7 (two box purchases and a correct Fermi–Dirac integration over an octant), dp1 (Mg2Sn labels), dp2 (TOTE columns) and dp8 (the lattice).
  - The instruction's lattice convention ("lattice constant below 100 K equals its 100 K value") hints that the lattice matters, but leaves finding and using the run to the agent.
- **Not peer reviewed.** The paper had not been reviewed as of October 2026. The graded values are archive outputs or stated combinations of them, so they do not depend on the paper's final text.

## Disk

- **Peak:** about 330 MB in `~/.cache/ai-theorist/mg2x-midir-forbidden-transition` (the zip, 77 MB, and the full extraction, 249 MB). It was deleted after the gate and the Harbor run.
- **Kept:** `fetch.py` and `members.md5`.
- **Task package:** 0.8 MB.
