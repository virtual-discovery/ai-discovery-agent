# Construction of cnt-kohn-anomaly-nonadiabatic

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** Marazzi, Poncé, Charlier and Rignanese, Phys. Rev. Lett. 135, 256402 (2025), with Materials Cloud record vs5fd-dbv98 (DOI 10.24435/materialscloud:yx-x1, CC BY 4.0).
- **Provenance, checksums, transformations and what is hidden:** the task's `SOURCE.md`.

## Rebuild

From this folder, with compute in `ai-theorist-benchmark-runtime:20261007` via `run_in_image.sh` (cpuset `CPUSET`, default 4,5,20,21; cache mounted at `/cache`):

1. `python3 -B fetch.py [DEST]` downloads the four record files (27 MB) through the per-file API links and checks their md5. It unpacks them into `DEST/raw` and checks the 75 members listed in `members.md5`. The default DEST is `~/.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic`.
2. `bash run_in_image.sh extract.py /cache` writes `DEST/pool/` (28 MB): the coupling tables, all 37 self-energy tables (as arrays and as verbatim lines for ν = 1, 20 and 23), the bands, the phonons, the DFPT scan and the small text files.
3. `bash run_in_image.sh make_task.py` takes about 2 minutes. It:
   - checks the lab evaluator against the archived tables (maximum deviation 0.17 meV);
   - computes the references, robustness variants and ablations;
   - rebuilds `coupling_calculator.c` with gcc and a parametrised `main()`, and asserts that it reproduces `output.csv`;
   - writes the catalog, the lab data (`selfenergy.py` is copied from here), `private/hidden.json`, `hidden_values.json`, `solution/ablation_values.json`, the budget in both `public.json` files, and `calibration.json` here.
4. From the repository root: `.venv/bin/python -m benchmarks.tools.check_task cnt-kohn-anomaly-nonadiabatic --cpuset <set>` writes `check.json`.
5. `python3 -B write_design.py` writes `design.json`, including the calibration table and the oracle length taken from `check.json`.

## What the data showed (design inputs)

- **The evaluator can be validated exactly.** The archived g±q tables hold, on a 2000-point grid, the energy-ordered Wannier bands at k and k±Q and |g| for the four band pairs. Two steps reproduce every archived self-energy table (modes 1, 20 and 23, 50–1578.87 K, static, Re and Im) to ≤ 0.17 meV:
  - follow the two bands through their crossings and interpolate bands and |g|² onto 10⁶ points;
  - use spin factor 2, a single process (k → k+Q), all four pairs, η = k_BT/100 (from `specfun/input.py`) and μ = −1.0380 eV.

  The archived series A and B temperatures are mutually consistent. Without interpolation, the raw 2000-point grid gets the dynamic part badly wrong (Re Π(ω_Q) at 300 K is −19.4 instead of −26.6 meV). With 10⁵ points, Im Π at 50–100 K is wrong (η is not resolved), although high-temperature tables still agree.
- **The k−q tables give the same self-energy as k+q.** Summing both doubles it. The authors' coupling-plot script sums |g+q| + |g−q| "because they are equivalent in 1D", but their self-energies contain one process.
- **The Dirac crossing sits exactly at the EPW Fermi level (−1.0380 eV).** The lower Wannier band is therefore fully occupied at T = 0, and the 2k_F backscattering appears as the (1,2)/(2,1) pairs in the energy-ordered basis. The intraband pairs alone give about 3% of the static self-energy. The band-plot script uses −1.0426 eV, which changes nothing at the target temperatures.
- **The self-energy table header is wrong.** It names 10 columns in another order than the 8 data columns. Column 7 is ω-independent (static), and column 8 is Im Π. Reading the header puts Im Π(ω ≈ 0) ≈ 0 in place of the static self-energy.
- **Smearing schemes differ.** ω_Q in the self-energy and coupling data (35.28, 130.11 and 143.42 meV) is the cold-smearing (`mv`, 0.01 Ry) DFPT value. The DFPT smearing scan is Fermi–Dirac: at 10 mRy it gives 48.93 and 137.83 meV. The authors subtract the Fermi–Dirac static self-energy at 1578.87 K (0.01 Ry) from the cold-smearing ω_Q. Between 5 and 0.5 mRy the scan agrees with the EPW static self-energy; the 10 mRy point does not.
- **For the RBM, nonadiabaticity barely moves T\*.**
  - Pole (Ω → 0 limit with the linewidth term): 527.4 K.
  - Static (adiabatic) criterion: 538.9 K.
  - The authors' iteration: 528 K; the paper quotes about 525 K.

  The paper's contrast between about 220 K (DFPT smearing proxy) and 525 K comes from the reference frequency and smearing scheme, not from dynamics. Anchoring on the Fermi–Dirac 10 mRy frequency gives 329 K.
- **ν = 20.** The nonadiabatic Ω is flat and non-monotonic (125.2 meV at 100 K, a minimum of about 123.4 meV near 340 K). The adiabatic Ω is 106.3 meV at 100 K. The adiabatic instability of ν = 20 extrapolates to about 0.3 K, which these data cannot locate.
- **Threshold program.** The parameter called `mass` in `coupling_calculator.c` is a = 1/m* (function documentation; `fig1.py` uses the reversed labels). At m* = 2 it gives 60.33 (ν = 1), 151.17 (ν = 20) and 163.40 meV (ν = 23). Eq. (9) of the paper gives values about 3% higher. The paper text's 54.5, 129.0 and 137.5 meV match neither.
- **Spurious low-temperature roots.** At 50–70 K, Im Π of the RBM has a spike near 4–5 meV that creates spurious pole roots. T\* is therefore defined as the onset on cooling.

## Deviations from the design (candidates.json / DETAILS.md section 10)

- **Archive size.** The archive has 37 self-energy temperatures, not 49.
- **Graded fields.**
  - **Target 1** (nonadiabatic Ω of the A1 mode at 100 and 300 K, ±1.5 meV) is kept, with ±1.2 meV, for mode ν = 20 by archive index. The archive's A1(L)/A1(T) labels contradict each other between the paper text and the scripts. The linewidths at 100 and 300 K are added (±6%); they carry the mesh fork.
  - **Target 2** (adiabatic T\*) is replaced by the adiabatic Ω of ν = 20 at 100 K (±1.5 meV). For the RBM the static T\* (539 K) lies within 2% of the nonadiabatic one, so it does not discriminate. For ν = 20 it lies far below 50 K and is not recoverable.
  - **Target 3** (whether a nonadiabatic instability occurs above 50 K, and T\*) is graded as the RBM T\* (±8%, `null` allowed and failing). Mode 20 is stable under nearly every path, so its verdict carries no information.
  - **Target 4** (λ_c, ±10%) becomes g₀ of the archived program for ν = 1 and ν = 20 (±6%). ±10% would accept the paper-text 54.5 meV for the RBM.
- **The threshold model is not a lab request.** A lab model run at the CNT parameters would return the hidden g₀ directly. The archived C source and its output are purchasable instead, and the agent image has gcc.
- **Decision points.**
  - Kept from the design: dp1 (double counting) as task dp1/dp2, dp2 (static versus dynamic) as dp3, dp3 (k-mesh) as dp7, and dp4 (bookkeeping) as dp5/dp6.
  - Not a decision point: the design's dp5 (Dyson versus first order). The complex-pole convention is stated for uniqueness, so first order would simply violate it.
  - Dropped: the design's dp6 (mode or q selection). The archive has coupling only at Q, and the modes are named because the labels are inconsistent.
  - Added, from verified archive evidence: dp4 (header), dp8 (smearing scheme and DFPT proxy), dp9 (mass parameter) and dp10 (non-monotonic Ω, interpolation).
- **Released tables.** 50, 700, 800, 900, 1000, 1250 and 1578.87 K; the 70–600 K tables are withheld. The plotting scripts are not served because they contain the protocol; their data arrays (the DFPT scan, `FERMI_SHIFT`, the k±q comment) are served with attribution.

## Decision points (10)

| dp | Fork | Kind |
|---|---|---|
| dp1 | ω_Q already contains an electronic self-energy (double counting) | physics; DFPT scan evidence |
| dp2 | Which one: static at the DFPT smearing 0.01 Ry = 1578.87 K, not the same T, not 0.01 eV | metadata (degauss, units) |
| dp3 | Static versus dynamic for ν = 20 | model class |
| dp4 | Table header mislabels the columns | archive defect in the raw data |
| dp5 | One process (k+Q); the authors' comment sums k±q | validation against an archived table |
| dp6 | Energy-ordered band pairs: backscattering is (1,2)/(2,1) | coupling data |
| dp7 | Low-T convergence (η = k_BT/100): 10⁶ points; high-T validation does not transfer | regime; needs the 50 K table |
| dp8 | Cold-smearing ω_Q versus the Fermi–Dirac DFPT scan; scan as temperature proxy | inconsistent archived calculations |
| dp9 | Threshold program parameter a = 1/m* | archive naming inconsistency |
| dp10 | Non-monotonic Ω(T): compute withheld temperatures, do not interpolate | regime judgement, binding budget |

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Defensible spread (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| Ω ν = 20, 100 K (meV) | 125.18 | −0.71 (no γ² term) to +0.03 (authors' iteration) | coarse mesh +0.34 (fails via γ), interpolation +0.25 (fails at 300 K), double processes −4.2, intraband only +5.0, raw grid +5.8, FD-scan ω_Q +8.2, adiabatic −18.9, same-T reference +21.8, no subtraction / header −25.2 | ±1.2 meV |
| Ω ν = 20, 300 K (meV) | 123.49 | −0.51 to +0.03 | interpolation +1.93, double −6.9, adiabatic −7.2, raw grid +7.1 | ±1.2 meV |
| γ ν = 20, 100 K (meV) | 12.866 | −3.8% (on-shell −Im Π(ω_Q)) to +0.3% | coarse mesh −9.8%, same-T −14.8%, no subtraction +24.6%, double +107%, raw grid −97% | ±6% |
| γ ν = 20, 300 K (meV) | 10.853 | −2.9% to +0.2% | interpolation −5.3% (passes; not the failing field), same-T −6.1%, smearing in eV −9.1%, no subtraction +12% | ±6% |
| Ω_ad ν = 20, 100 K (meV) | 106.27 | −0.47 (dynamic reference) to +0.34 (log fit of released tables) | FD-scan ω_Q +7.9, NA +18.9, intraband +22.8, same-T +23.8, no subtraction −28.6, double −31.0 | ±1.5 meV |
| RBM T\* (K) | 527.4 | 527.7 (lab 10⁶) to 538.9 (static criterion), 531.4 (log fit), 528.0 (authors) | FD-scan ω_Q 328.9 (−38%), interpolation 367.0 (−30%), DFPT proxy 203.7, double 910, no subtraction ≥ 1579, same-T / eV / intraband: no instability | ±8% |
| g₀ ν = 1 (meV) | 60.33 | Eq. (9) +3.1%, output.csv interpolation −3.7% | a = m* +101%; paper text 54.5 (−9.7%, not an ablation) | ±6% |
| g₀ ν = 20 (meV) | 151.17 | Eq. (9) +3.1%, output.csv +0.1% | a = m* +102%; paper text 129.0 (−14.7%) | ±6% |

Rounding to 4 significant digits passes.

## Oracle

The oracle makes 18 evaluations for 21 of 26 units (81%):

- `scf.in` and `specfun/input.py`: 2 units;
- the DFPT scan: 1 unit;
- tables of ν = 20 and ν = 1 at 1578.87 K and of ν = 20 at 50 K: 3 units;
- 50 K validation at 10⁶ and 10⁵ points: 3 units;
- ν = 20 at 100 and 300 K (10⁶ points, 40 frequencies plus the static point): 4 units;
- RBM static check at 1578.87 K plus 500, 520, 540 and 560 K (10⁵ points): 5 units;
- `coupling_calculator.c`, compiled with gcc: 1 unit;
- coupling tables of ν = 1 and ν = 20 (|g| at k_F): 2 units.

The oracle identifies the static column from the data (the ω-independent one) and asserts the 50 K validation. It solves the complex pole by interpolation and takes T\* from the Ω → 0 limit.

Graded errors: Ω −0.045 and −0.005 meV; γ +0.1% and −0.1%; Ω_ad +0.012 meV; T\* −0.3%; g₀ 0.

A worst-case evaluator call (4×10⁶ points, 2000 K, 41 frequencies, both processes) takes 33 s and 0.8 GB. The catalog is unbounded because the evaluator can be bought repeatedly; the 21 tables of the three modes alone cost 21 units. An agent needs roughly 45–75 tool calls.

## Results

- **Gate (2026-10-20, cpuset 4,5,20,21):** passed.
  - 10 decision points; all 14 ablations fail.
  - 9 robustness variants and 4-significant-digit rounding pass.
  - 8 rubric items and 8 hidden values, with no leaks.
  - The package is 4.2 MB; the oracle took 35 s in the gate container.
- **Harbor oracle (`runs/build-cnt-kohn-anomaly-nonadiabatic/oracle-r1-cnt-kohn-anomaly-nonadiabatic`):** reward 1.0, solver charge 21, runtime 62 s.

## Open risks

- **Dev-easy risk.** The pole and bare-frequency conventions are stated for uniqueness. An agent that knows the double-counting correction and validates the evaluator against the 50 K and 1578.87 K tables defuses dp1–dp3 and dp5–dp7. The discriminating traps are:
  - dp4: a header that contradicts the data;
  - dp5: the authors' own k±q comment;
  - dp7: cheap meshes, which validate at high T;
  - dp8: the Fermi–Dirac scan, which looks like the consistent reference;
  - dp10: the non-monotonic Ω.
- **Memorization.** The paper's RBM T\* (about 525 K) passes the T\* check. The paper's thresholds (54.5 and 129 meV) fail. The ν = 20 frequencies appear only in a figure.
- **Thin margins.**
  - Interpolation at 300 K fails by 0.73 meV beyond tolerance.
  - The coarse mesh fails the 100 K linewidth by 3.8 points.
  - The no-γ² variant passes with 0.49 meV to spare.
  - The on-shell linewidth passes with 2.2 points to spare.
- **Authors' protocol as truth.** The reference follows the authors' protocol, which combines a cold-smearing ω_Q with a Fermi–Dirac subtraction. The instruction pins ω_Q, so an agent that "repairs" the scheme mismatch fails T\* (329 K). The rubric asks the report to discuss the mismatch.
- **Evaluator fidelity.** The lab evaluator interpolates bands that are rounded to 0.1 meV. Its residual 0.03–0.17 meV deviation is far below every tolerance. The spike in Im Π of the RBM near 4–5 meV at 50–70 K is an archived feature, not reproduced exactly.

## Disk

- **Peak:** about 265 MB:
  - in `~/.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic`, the 27 MB of archives, the 138 MB unpacked tree and the 28 MB pool;
  - in the scratch subfolder, 70 MB (a pickle of the tables and the arXiv PDF).
- **After the Harbor run:** all were deleted, together with the gate's lab image. `fetch.py`, `members.md5` and `extract.py` regenerate them.
