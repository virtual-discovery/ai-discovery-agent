# Construction of mn5si3-noncollinear-spinwaves

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). Sources: Biniskos, dos Santos et al., PRB 105, 104404 (2022) and APL Mater. 11, 081103 (2023); Materials Cloud records paxx3-a1k09 and 4tgt6-1p378 (CC BY 4.0).

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g)`), with the repository and the cache mounted at their host paths. CACHE = `~/.cache/ai-theorist/mn5si3-noncollinear-spinwaves`.

1. `python3 -B fetch.py CACHE` downloads both `Mn5Si3.zip` files (126 MB) through the per-file API, checks the record md5s and unpacks only the spin-wave, neutron and `DFT/inputcard` members.
2. `python3 -B extract.py CACHE` writes `CACHE/pool.json`: structure, exchange tables with bond vectors in Å, dispersion tables, cuts, energy scans, gap table, Q scans, the phase diagram and the archived ground state.
3. `python3 -B make_task.py CACHE` writes `data.json`, `catalog.json`, `private/hidden.json`, `hidden_values.json` and both `public.json` files.
4. `python3 -B calibrate.py` (~15 min on 4 cores) runs the reference route, the defensible variants and the wrong paths from purchasable data only. It writes `calibration_table.json` and `solution/ablation_values.json`. `calibrate.py rescore` re-applies changed thresholds.
5. `.venv/bin/python -m benchmarks.tools.check_task mn5si3-noncollinear-spinwaves --cpuset 4,5,20,21` writes `check.json`.

`solution/mnsi_lswt.py` is the LSWT library: Landau–Lifshitz linearisation in local frames, symplectic (Colpa-type) normalisation with degenerate clusters, ground-state relaxation, and the scan/peak-extraction model. It reproduces the authors' SWIS eigenvalues exactly (archived `dispersion_MnSiDisper.dat`: same table, state and anisotropy, S = 1.22). The intensities were checked against the analytic two-sublattice chain.

## What the archive showed (design inputs)

- **Convention.** SWIS energies follow E = −Σ_listed J e_i·e_j − Σ K (n·e)², with S = μ_s/γ = 1.22, where the table lists every bond from both ends.
- **Bond-offset bug in the authors' run.** The `_temp` 4.5 Å table used for the paper's spin-wave run has cell offsets that do not match its own archived Cartesian bond vectors (2024 `pair_MnSi_rad4.5noncol.txt`, columns 10–13). With the offsets as archived, several "2.9/4.4 Å" bonds actually connect sites 5.8–17 Å apart. The task sells bond vectors, re-derived from the Cartesian columns and checked; offsets are not sold. The archived 16-site ground state is stationary, but with the archived offsets it is unstable at (½,0,0).
- **Mn1 couplings** (−2.17 meV Mn1–Mn2a, +11.99 meV Mn1–Mn1 along c) are the authors' choices, not KKR results. The catalog says so.
- **DFT model quality.** No exchange table fits all directions within the archived errors: h prefers an energy scale of about 0.8–0.9, k prefers 1.0–1.1, and the k plateau near 2.75 r.l.u. (9–17 meV) is not reproduced. A free 5-parameter fit to the cut line shapes reached χ²r ≈ 6 on h/k but predicted l badly (rms 0.10–0.14): overfitting.
- **Isotropic cone.** In Å⁻¹ the cone is nearly isotropic, so a phenomenological cone fitted to h and k predicts the l branches with rms 0.024. The l prediction alone therefore cannot discriminate microscopic models (see Deviations).
- **(2,1,0).** The withheld 5–20 meV energy scan at (2,1,0) has a peak at 9.84 ± 0.97 meV (archived). In the Mn1 + Mn2a models this is the strong Mn1–Mn2a optical mode (9.3–14.5 meV across tables and scales). The Mn2-only model has nothing between 5 and 20 meV; its next mode is at 32 meV.
- **AFM2 (80 K) gaps** of 0.18 and 0.41 meV are far below any model with the anisotropy that fits the 10 K gap of 0.71 meV. They stay in the catalog as context.

## Deviations from the design (candidates.json / DETAILS.md section 23)

- **Lab type.** It is a lookup lab. It sells archived exchange tables (not "pair sets for requested orderings") and archived neutron data. The agent relaxes ground states and runs LSWT itself; the catalogue is worth 81 units against a budget of 25.
- **Graded targets.**
  - (1) The third direction is kept: all 26 (0,3,Ql) branch positions at 13 energies, with rms ≤ 0.065 r.l.u.
  - (2) "α/β mode energies within 0.5 meV" was dropped. Those are 80 K AFM2 fits, and no honest route reproduces them. In their place, the withheld (2,1,0) zone-centre excitation is graded in a window of 8–15.5 meV.
  - (3) The canting is graded in a 9–19° window against the authors' archived ground state (13.93°). The archived state is a model output, but any no-recall relaxation of any sold table gives 12–15.6°.
  - (4) A categorical Mn1 answer was added.
- **Reachability without recall.** The design's "rescale DFT J uniformly" and "first two shells only" forks did not survive as discriminating ablations. Overall scale is calibrated from purchased data, and table range or reference changes all pass. They are robustness variants instead.
- **Difficulty.** The l-branch target is weakly discriminating, because an isotropic cone passes. The model-based content is carried by the optical mode, the Mn1 decision and the canting.

## Decision points (6)

| dp | Wrong path | Fails on |
|---|---|---|
| dp1 | Mn2a-only sites (AFM2-like) | optical (32.1 meV), Mn1, canting 0 |
| dp2 | Bonds double counted | l rms 0.067, optical 6.9 |
| dp3 | Free-ion S = 5/2 | l rms 0.128, optical 6.2 |
| dp4 | Lowest branch / equal branch weights | l rms 0.091 / 0.148 (+ optical 6.4) |
| dp5 | Collinear state assumed, not relaxed | canting 0 |
| dp6 | Angle to the Mn1 axis instead of the tilt | canting 76 |

## Calibration

| Field | Reference | Defensible spread | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| l rms (r.l.u.) | 0.032 | 0.027–0.041 (cone 0.024; no anisotropy 0.048) | double counting 0.067 | ≤ 0.065 |
| (2,1,0) mode (meV) | 12.76 (truth 9.84) | 9.30–14.54 | double counting 6.90; Mn2-only 32.1 | [8.0, 15.5] |
| canting (deg) | 13.90 (archived 13.93) | 12.05–15.57 | collinear 0; angle-to-axis 76.1 | [9, 19] |

Tightest margins:

- **l rms:** the double-counting wrong path fails by only 0.0016. The no-anisotropy variant (0.048) and the collinear-reference table (0.041) pass by 0.017–0.024.
- **Optical window:** the collinear-reference table (14.54) passes by 0.96. The scale-1.1 variant (14.04) passes by 1.5.

## Oracle

- **Purchases:** 4 exchange sets, the h and k tables, 6 cuts (2 units each), 3 energy scans and the gap table. That is 16 evaluations and 22 of 25 units (88%).
- **Run time:** container time 75–78 s.
- **Gate (2026-10-26):** passed. All 7 ablations fail; the 8 robustness variants and 4-significant-digit rounding pass.
- **Harbor oracle:** reward 1.0 (`runs/build-mn5si3-noncollinear-spinwaves/oracle-r1-*/*/verifier/reward.txt`, 1 min 58 s).

## Open risks

- **Easy l target:** a cone fit passes the l branches, so the target is easy on its own.
- **Literature recall:** an agent that recalls "AFM1 with Mn1 moments" passes the Mn1 and canting checks without the analysis. The optical mode still requires the Mn1 model with a calibrated scale.
- **Wide optical window:** the window is broad (8–15.5 meV) because the table and scale choices spread the predicted mode.
- **Shared double-counting margin:** the l tolerance sits close to the double-counting path (margin 0.0016), and that path also fails on the optical mode.
- **Peak disk:** about 0.45 GB of cache during the build. Raw downloads are deleted; `pool.json` (0.4 MB) and the record JSONs are kept.
