# Construction of sngeo2-coherent-miscibility

Released-start task built on Müller, Kurbay, Zhang, Kioupakis & Natarajan, arXiv:2601.12184 (2026). The data are Materials Cloud record 1qt4h-tr585 (DOI 10.24435/materialscloud:e2-xp, CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. `python3 -B fetch.py` downloads `data.tar.xz` (8.6 MB) and `README.txt` to `~/.cache/ai-theorist/sngeo2-coherent-miscibility/`. It verifies md5 and sha256, then unpacks to `x/` (163 MB, mostly HSE band-structure outputs that are not used).
2. `python3 -B extract.py` writes `pool.json.gz` (0.8 MB). It contains:
   - the orderings, SQS and strain scans of both functionals, with structures as lattice, species and fractional coordinates;
   - the HSE-setting PBE runs.

   It also asserts the archive relations the task relies on (see `SOURCE.md`).
3. `python3 -B make_task.py` writes:
   - `environment/agent/catalog.json` and the 80 PBE ordering POSCARs in `environment/agent/orderings_PBE/`;
   - `environment/lab/data.json`;
   - `private/hidden.json` and `private/hidden_values.json`;
   - `solution/ablation_values.json`.

   It prints the reference, ablation and robustness payloads.
4. `python3 -B calibrate.py` prints the spread of the correct method over purchase plans (`calibration.json`).
5. `ce_mc/run_all.sh` runs the cluster expansion and Monte Carlo used for the incoherent-Tc window. Run it in the runtime image after step 1; the logs are in `ce_mc/logs/`.
6. `.venv/bin/python -m benchmarks.tools.check_task sngeo2-coherent-miscibility --cpuset 4,5,20,21` writes `check.json`.

`model.py` holds the thermodynamics:

- SQS screening and subregular fits;
- strain-scan minima and Δe_cs(x) on the Vegard line;
- spinodal, common-tangent binodal and critical temperatures;
- `analysis(**choices)`, which produces the reference and every wrong path from one code path.

`solution/solve.py` re-implements the oracle stand-alone from purchases only.

## What the released data show

### Clustering, not ordering

All 80 small-cell orderings and all 21 SQS per functional have positive formation energies. The PBE mixing enthalpy is about 110–115 meV/cation at x = 0.5.

### SQS units (dp1)

`data_SQS_*.json` holds total energies (eV per 6-atom cell) for x = 0 and 1. For 0 < x < 1 it holds the formation energy per rutile primitive cell (2 cations, the paper's Eq. 1); the README does not say so.

- The subregular fit of these numbers reproduces the paper's Table 1 exactly (PBE α = 0.915, β = −0.117 eV; r2SCAN 1.071, −0.219 eV).
- The archive's separate PBE relaxations of the same SQS (`Alloy_PBE_Values`, hybrid-functional settings) give 105.4–110.6 meV/cation at x = 0.5. That is half the released 0.212–0.223 eV, within 1 meV, so the purchase that settles the units exists.

### Unrelaxed r2SCAN SQS (dp2)

Three r2SCAN SQS did not reach their equilibrium volume: x = 0.625 index 1, and x = 0.75 indices 0 and 1.

- Their decorations are identical to the PBE SQS (checked atom by atom).
- For the other 18 SQS the r2SCAN/PBE volume ratio is 0.960; for these three it is 0.936. They lie about 0.7 Å³/cation (2.4 %) below the end-member Vegard line.
- Their energies are about 20 meV/cation above the trend: 0.281 / 0.263 / 0.264 eV per cell, against 0.242 and 0.225 for the valid SQS at those compositions.

The authors' r2SCAN Table 1 includes them. Excluding them moves the r2SCAN coherent spinodal from 904 K at x = 0.67 to 665 K at x = 0.42. The graded truth excludes them; the paper's 890 K is the wrong path.

### Strain datasets

- **Same settings as the SQS.** The scan minima equal the SQS end-member energies to 1e-5 eV.
- **Grids.** The in-plane grids sit on the end-member Vegard line, in steps of x = 0.05 for (001) and 0.1 for the 2D grids.
- **The archive's "101" dataset.** It fixes c and |a2 − a3| and steps |a2 + a3|: a (110)-type plane in the standard setting. The label is kept, and the point is rubric item 5.
- **Reproducing Table 2.** Parabolic scan minima give α_cs = 0.619 / 1.200 / 1.861 eV (PBE) against the paper's 0.601 / 1.176 / 1.833, and 0.811 / 1.361 / 2.133 (r2SCAN) against 0.801 / 1.331 / 2.108. The values are 2–3 % higher; the graded strain tolerance of ±5 % admits the paper's values.

### Incoherent Tc is not unique

These runs are the coordinator's Fe lesson again.

| Estimate | Tc |
|---|---|
| Ideal-entropy (mean-field) estimate on the SQS enthalpy | 2746 K (the paper's ≈2750 K) |
| Lasso CE with 24 clusters (RMSE 4.4, 5-fold CV 7.2 meV/cation) + SGC Monte Carlo, Binder crossings L = 6/8 | about 2550 K |
| 17 clusters (CV 7.6) | about 2350–2400 K |
| 12 clusters (CV 11.4) | about 2250 K |
| 8 clusters (CV 20) | about 1650 K |
| The paper's CE | ≈2350 K |

The gap between mean field and validated Monte Carlo (about 7–15 %) is no larger than the spread between validated CEs. So the incoherent Tc is graded only as a window, [2000, 3000] K. The window admits mean field and validated CE + MC. It rejects the unit error (5490 K), a coherent value (≈870 K) and an over-truncated CE (≈1650 K, dp9).

The CE's random-alloy enthalpy is its empty-cluster ECI, 121 meV/cation for every λ tried. The paper's CE gives about 106. This number is not graded either.

## Deviations from the candidate design

- **Incoherent Tc** is graded as a window, not ±15 % of a CE + MC value (see above).
- **Ge solubility at 1000 K** was replaced by the 700 K coherent binodal and spinodal (PBE). At 1000 K every method gives incoherent solubility below 1 %, and the coherent alloy is fully miscible, so the original target did not discriminate.
- **Mixing-enthalpy tolerance** is the window [100, 124] meV/cation instead of ±5 meV. Defensible estimates span 106–121: single SQS, SQS mean, subregular fit and CE random alloy.
- **r2SCAN coherent Tc** is graded with the unrelaxed SQS excluded (see dp2), and the ±20 % tolerance reflects purchase-plan scatter.
- **Orderings.** The 80 PBE orderings are staged free as the CE training set, so the budget binds on SQS, strain scans and the units check. They are staged as POSCAR files (the archive's format), with the total energy in the comment line.
  - The gate's leak scan had flagged coincidental 4-digit substrings in `catalog.json` when the orderings were inlined there. A total energy of −114.4398 eV matched the 114.4 meV/cation target, and a coordinate 0.1280967 matched the 0.128 binodal composition.
  - Neither is derivable from the other. The matches were checked by hand.
- **Coherent model.** The coherent model (Cahn approximation with ideal entropy) is defined in the prompt, as the design asked, because otherwise the coherent numbers are not unique.

## Calibration

The full table is in `design.json` → `calibration`. Purchase-plan percentiles come from `calibration.json`, with one random valid SQS per composition for 3, 4, 5 or 7 compositions.

| Check | Reference | Robustness / plans (5–95 %) | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| mixing enthalpy x = 0.5 (meV/cation) | 114.4 | 106–121 | 129 (r2SCAN), 229 (per-cation reading) | [100, 124] |
| Δe_cs(0.5) 001 / 100 / 101 (meV/cation) | 78.0 / 149.7 / 232.5 | paper α/8: 75.1 / 147.0 / 229.1 | no relaxation 104.9 / 223.6 / 246.5; c unconstrained 78 / 21.8 / 54.0 | ±5 % |
| habit plane | 001 | 001 | 100 (c unconstrained) | exact |
| incoherent Tc (K) | 2746 (mean field) | CE + MC 2250–2550 | 1650 (over-truncated CE), 3064 (r2SCAN), 5491 (units) | [2000, 3000] |
| coherent Tc PBE (K) | 866 | 816–950 | 614 (no relaxation), 311 (functional mixing), 2299 (c unconstrained) | ±12 % |
| coherent Tc r2SCAN (K) | 665 | 647–788 | 904 (unrelaxed SQS kept), 960 (index 0 unscreened) | ±20 % |
| 700 K binodal (PBE) | [0.128, 0.818] | lo 0.089–0.182, hi 0.80–0.898 | swapped [0.246, 0.677] | ±0.10 each |
| 700 K spinodal (PBE) | [0.246, 0.677] | lo 0.215–0.322, hi 0.658–0.773 | swapped [0.128, 0.818] | ±0.10 each |

Tightest margins:

- **Spinodal upper composition.** About 5 % of 3-composition plans exceed ±0.10 (the 95th percentile is 0.773 against 0.677).
- **Binodal/spinodal swap.** It misses by 0.118 against the 0.10 tolerance.
- **"101" with no relaxation.** It sits 6.0 % from the reference, outside ±5 %, and that path also fails on 001 and 100.
- **r2SCAN coherent Tc.** The window upper edge is 798 K against the 788 K 95th percentile of 3-composition plans; the nearest wrong path is 904 K.

Rounding every number to 4 significant digits passes.

## Oracle and length

The oracle uses 37 evaluations and 52 of 58 units (90 %):

- **Units check:** 3 HSE-setting runs (both end members and `Sn_50_1`).
- **PBE SQS:** 7 (index 0 at every composition) plus 2 more at x = 0.5.
- **r2SCAN SQS:** 9. Index 0 at every composition is screened against the Vegard line; x = 0.75 needs indices 1 and 2.
- **PBE strain scans:** 10. All three planes at x = 0.5, plus (001) at x = 0.25 and 0.75.
- **r2SCAN strain scans:** 6, on (001) at x = 0.25, 0.5 and 0.75.

The oracle reports the mean-field incoherent Tc and states the Monte Carlo correction in the report. It takes 6 s in the gate container; the Harbor run took 48 s.

The catalog totals about 1840 units, so the budget covers 3 %. The expected agent length is 50–80 tool calls: structure parsing, the units cross-check, SQS screening, strain-grid mapping, free-energy numerics and, optionally, CE + MC.

## Gate and Harbor

- **Gate (2026-10-26):** passed.
  - 9 decision points;
  - all 11 ablations fail;
  - all 7 robustness variants pass;
  - rounding passes;
  - 11 hidden values, no leaks.
- **Harbor oracle:** reward 1.0, `runs/build-sngeo2-coherent-miscibility/oracle-r1-sngeo2-coherent-miscibility`.

## Open risks

- **Dev-easy risk.** The hard parts are evidence-based: the SQS units, the three unrelaxed r2SCAN SQS, and the in-plane constraint on the 2D grids. A careful agent that cross-checks units with the HSE-setting runs and screens volumes will pass. The coherent model itself is defined in the prompt.
- **Memorisation.** The paper's PBE coherent Tc (930 K) lies inside the PBE window, and its "above 2300 K" lies inside the incoherent window. Its r2SCAN 890 K fails, and the strain energies, the 700 K boundaries and the units are not in the paper in graded form.
- **Exclusion of the r2SCAN SQS.** The graded r2SCAN truth departs from the authors' fit. The evidence is strong (identical decorations, volume ratio and Vegard deviation, energy outliers), and it is stated in `SOURCE.md` and rubric item 3.
- **Incoherent window.** A CE with CV around 11–15 meV/cation sits near 2000–2250 K, close to the lower edge.

## Disk

- **Peak:** about 175 MB in `~/.cache/ai-theorist/sngeo2-coherent-miscibility/` (the archive, its unpacked tree and the arXiv HTML read for the paper's tables). Everything was deleted after the Harbor run.
- **Kept:** `pool.json.gz` (0.8 MB) and the fetch script with checksums.
- **Task package:** 5.6 MB, of which `environment/lab/data.json` is 5.5 MB.
