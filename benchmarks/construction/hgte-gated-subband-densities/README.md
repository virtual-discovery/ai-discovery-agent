# Construction of hgte-gated-subband-densities

This is a released-start task built from M. Hofer et al., *Self-consistent electrostatic modeling of gated narrow-gap topological insulators*, Nano Lett. 26, 3042 (2026). It uses two Zenodo records, both CC BY 4.0:

- **18709991:** kdotpy runs and processed SdH charts.
- **18429111:** raw mK transport.

Provenance and checksums are in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

1. **Fetch.** `python3 -B fetch.py [DEST]` reads the 397 needed members through HTTP range requests on the two remote zips. It skips the 0.8 GB and 0.5 GB dispersion CSVs and the Origin projects. It checks every member against `members.md5`. The download is about 800 MB unpacked, and the default DEST is `~/.cache/ai-theorist/hgte-gated-subband-densities`.
2. **Extract.** `bash run_in_image.sh extract.py` writes `/scratch/pool.npz` (36 MB, in DEST/scratch). It holds:
   - the transport grids, B ≤ 2 T;
   - all 80 FigureS3 kdotpy runs (`.kpout`, potential, band-resolved IDOS).
3. **Make the task.** `bash run_in_image.sh make_task.py > make_task.log` writes the lab pool, catalog, private targets, ablation and robustness values, and `calibration.json`. On the way it asserts four things:
   - the authors' 70 nm gate line is reproduced exactly;
   - the archived 70 nm SdH chart is reproduced from the raw record to 1e-10;
   - the wedge runs satisfy band −1 = 2 × (majority + minority);
   - every robustness variant passes and every ablation fails.
4. **Gate and design.** Run in this order:
   1. `python3 -B write_design.py`
   2. `.venv/bin/python -m benchmarks.tools.check_task hgte-gated-subband-densities --cpuset 4,5,20,21`
   3. `python3 -B record_oracle.py`
   4. `python3 -B write_design.py` again, to record the oracle length.

Other files:

- **`sdh.py`:** the SdH Fourier analysis, an exact port of `FFT_experiment.ipynb` plus variants.
- **`model.py`:** readers for the kdotpy runs, following `FigureS3.ipynb`.
- **`explore_calib.py`, `explore_ridge.py`:** the exploration scripts used to find the ridge spreads and the hole-ridge line search.
- **`targets_70nm.png`:** the archived 70 nm chart with the graded points. It is private and not staged.

## What the task is

The lab serves archived data only:

| Request | Cost | Returns |
|---|---|---|
| `hall_scan` | 1 | Low-field (|B| ≤ 0.2 T) slope of Vxy/VI and Vxx/VI at B ≈ 0, at all gate voltages of one sample. |
| `transport` | 1 | Raw lock-in rows (VG, B ≤ 2 T, VI, Vxx, Vxy) for all gate voltages within ±0.025 V of a 0.05 V-grid centre. |
| `calculation` | 1 | One main self-consistent kdotpy run: `.kpout`, V(z), IDOS of bands −3…6 and the total, −60…80 meV. |
| `pocket_calculation` | 1 | One kphi-wedge run (105–165° or 165–285°) in the same potential. |

There are three samples: 26, 70 and 107 nm. The catalog totals 162 units and the budget is 25.

The graded values, all for 70 nm:

- **Six SdH densities:**
  - top surface at +0.20 V;
  - bottom surface at +0.05 V;
  - first bulk subband at +0.30 V;
  - top and bottom surface at −0.45 V;
  - majority hole pocket at −0.70 V.
- **The hole-onset gate voltage.**
- **Three model densities:**
  - bottom surface at +0.20 V;
  - top surface at −0.45 V;
  - majority pocket at −0.60 V.
- **The majority-pocket multiplicity, 2.**

## What the archive shows (and the task is built on)

- **Exact SdH reproduction.** `sdh.chart` with the notebook's parameters reproduces `70nm_FT.csv` from `70nm.csv` to 1e-14: σxx, B ≤ 0.8 T, Hamming window, linear detrend, n = e f/h on a 0.05e11 grid. The chart axis is already e f/h, so each peak is one spin-resolved pocket.
- **Gate line.** The authors' 70 nm line (9.8197 Vg + 3.2060) is the single-carrier Hall density, |B| ≤ 0.2 T, fitted over 0.05 < Vg ≤ 0.29 V. The window is sample-specific: 26 and 107 nm use 0.2–0.38 V. With other windows the 70 nm single-carrier line moves a lot (CNP from −0.33 to −0.44 V), because a low-mobility bulk subband fills above ~0 V. The convention is therefore stated in the instruction.
- **Surface labels in kdotpy.** The bands are numbered by energy. On the electron side, band 1 is the gate-responsive top surface (1.9→3.9) and band 2 the screened bottom surface (0.94→1.07). Their labels cross near n ≈ 0.001 nm⁻². On the hole side, band 2 is the top surface (0.57→0.33), and the authors swap the labels below Vg = −0.3 V.
- **Hole pockets.** The wedge runs give single-pocket densities as IDOS(−1) at the main run's E_F divided by 360°/width (6 or 3). kdotpy extrapolates a wedge as if isotropic, so each wedge run's own `[ef]` is wrong for an anisotropic band. The main run's band −1 = 2 × (p[11] + p[1̄1]) to ≤ 0.5 % for n ≤ −0.002 nm⁻² (2–5 % near the onset, where the wedge sampling is coarse), which means two inequivalent pairs; bulk inversion asymmetry splits the four side maxima. The [11] pocket dominates: 2.19 vs 0.32 at n = −0.004 nm⁻².
- **Experiment vs model (70 nm).**
  - Top surface: SdH 2.81 vs model 2.76 at +0.20 V.
  - First bulk subband: 1.22 vs 1.20 at +0.30 V.
  - Bottom surface: 0.75–0.83 vs model 0.94–1.0, so the model is too high.
  - Majority hole pocket: SdH 2.49 vs model 2.04 at −0.70 V, so the model is ~20 % too low.
  - Charge balance at −0.70 V: holes = 0.28 + 0.64 + 3.66 = 4.6 ≈ 1.85 × 2.49.
- **Hole onset.** The SdH ridge of the majority pocket (slope ≈ 5.1–5.2 × 10^11 cm⁻² per V, about half the gate lever arm) extrapolates to zero at −0.21 V, on the electron side of the CNP (−0.33 V). The model onset is −0.21 V. At single gate voltages the ridge is weaker than harmonics and sum frequencies of the surface states: the strongest peak above 1e11 at −0.70 V is at 1.17. Between −0.50 and −0.40 V it is lost among them, so the onset fit uses −0.75…−0.55 V.

## Deviations from the design (candidates.json / DETAILS.md section 51)

- **No kdotpy on request.** The design wanted a lab that runs kdotpy, or a parameter grid, over ε, VBO, E_P, boundary condition, self-consistent on/off and effective-mass vs full-band neutrality. The archive holds one parameter set per sample (E_P = 21 eV, VBO = −620 meV, ε = 20.8, sample-specific dv1). The authors' runs need GPUs (jax_eigh, 4–6 GPUs) and CPU-days for a grid. A numpy 8-band Schrödinger–Poisson evaluator validated to kdotpy's material database, strain and BIA terms is not feasible within a few CPU-hours of build. So there is no calibration fork for ε/VBO and no fork for effective-mass vs full-band neutrality. These are rubric items (4, 8). The self-consistency fork survives as the rigid-band ablation of an archived run (dp7).
- **Targets.**
  - The 70 nm SdH densities are still graded as designed, but the agent extracts them from the raw 70 nm record rather than predicting them from a calibrated model. The archived model misses some of them by 15–20 % (bottom surface, hole pocket), so a model-only prediction is not recoverable within any fair tolerance.
  - The model side is graded separately: three densities of the released 70 nm runs under a stated gate convention.
  - The CNP voltage was dropped. It is a convention of the gate fit window (the spread across windows is 0.1 V).
  - The hole-onset target is graded as a gate voltage, not a density: the onset density is close to 0 and has no stable relative tolerance.
- **Multiplicity graded explicitly.** The majority-pocket multiplicity (design dp7) is graded as an exact integer.
- **Samples.** 45 nm (CompMat-Bench) is not used at all. 26 and 107 nm are sold for calibration only.

## Decision points (8)

| dp | Fork | Kind |
|---|---|---|
| dp1 | n = e f/h for spin-split bulk subbands and pockets vs 2e f/h | conversion, checkable by charge balance |
| dp2 | Hole ridge vs harmonics and sum frequencies (strongest peak 1.17, nearest-to-model 2.15) | evidence across gate voltages; costs budget |
| dp3 | Top vs bottom surface on the hole side; kdotpy label crossing | discovery in the model runs and data |
| dp4 | Majority-pocket multiplicity: 2 (BIA) vs 4 (cubic) vs 1 | charge balance plus wedge runs |
| dp5 | Wedge run read at its own E_F vs the main run's E_F | metadata trap in the `.kpout` |
| dp6 | Wedge normalisation 360°/width; band −1 is not 2 or 4 equal pockets | inference from the command lines |
| dp7 | Self-consistent runs vs rigid-band filling of one run | model class; budget temptation |
| dp8 | Hole onset vs CNP / Hall sign change | regime judgement |

Cascades:

- dp2 feeds the hole density, the onset and the multiplicity.
- dp3 feeds both experimental values at −0.45 V and the model top surface.
- dp5 and dp6 together determine the model pocket.

## Calibration

The full table is in `design.json` → `calibration` and `calibration.json`.

| Field | Reference | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| top +0.20 | 2.807 | 2.803–2.819 | — | ±0.168 |
| bottom +0.05 | 0.827 | 0.824–0.835 | — | ±0.07 |
| first bulk +0.30 | 1.224 | 1.196–1.232 | 2.448 (2e f/h) | ±0.073 |
| top −0.45 | 0.446 | 0.433–0.452 | 0.631 (swapped) | ±0.07 |
| bottom −0.45 | 0.631 | 0.574–0.643 | 0.446 (swapped) | ±0.07 |
| hole −0.70 | 2.491 | 2.411–2.511 | 2.147 (nearest model), 1.172 (strongest), 4.98 (2e f/h) | ±0.150 |
| hole onset (V) | −0.214 | −0.272 to −0.187; model −0.211 | −0.327 (CNP), −0.350 (Hall sign change) | ±0.06 |
| model bottom +0.20 | 0.971 | 0.970, 0.971 | 1.592 (rigid band from n = 0) | ±0.058 |
| model top −0.45 | 0.485 | 0.483, 0.486 | 0.675 (band 1) | ±0.04 |
| model pocket −0.60 | 1.604 | 1.595, 1.601 | 1.858 (band −1 / 2), 1.255 (wedge E_F), 0.929, 4.81 | ±0.096 |
| multiplicity | 2 | — | 4, 1 | exact |

The SdH variants are:

- single spectrum vs ±0.01/0.02/0.024 V averages;
- B_max 0.6/0.8/1.0/1.5 T, with zero-padding;
- σxx vs ρxx;
- Hamming vs Hann window.

The model variants are the gate-fit window (0.06–0.28 V) and quadratic interpolation.

Tightest margins:

- The hole onset with B_max = 0.6 T is −0.272, 0.058 V from the reference against a 0.06 V tolerance.
- The bottom surface at −0.45 V with B_max = 0.6 T is 0.574, 0.057 below the reference against 0.07.
- The model top surface at −0.45 V: the experimental top value 0.446 would pass the model check (0.039 < 0.04).

## Oracle

The oracle uses no literature value. It buys:

- the 70 nm Hall scan;
- 9 main runs: 0.003…0.007 and −0.001…−0.004 nm⁻²;
- 3 wedge runs: 165–285° at −0.002 and −0.003, and 105–165° at −0.003;
- 9 transport blocks: +0.30, +0.20, +0.05, −0.45, and −0.55 to −0.75.

That is 22 evaluations, 22 of 25 units (88 %). The oracle's method:

- **Assignment.** Each experimental peak is the strongest local maximum within ±25 % of the model value. The top surface is chosen as the band that follows the gate.
- **Hole ridge.** A line search on the peak prominence (log power minus a ±0.25e11 running mean) finds the ridge; an absolute-power search locks onto harmonics.
- **Multiplicity.** From charge balance.

Largest oracle errors: hole pocket +0.019 (tolerance 0.150), bottom +0.05 V +0.006, onset +0.006 V.

- **Gate (2026-10-27, cpuset 4,5,20,21):** passed. 13 ablations fail, 9 robustness variants pass, 4-significant-digit rounding passes and grading is deterministic. 8 decision points, 8 rubric items, 14.7 MB, 10 hidden values. The oracle container ran for 22 s. After the oracle length was written into `design.json` the gate was re-graded with `--reuse-state` (passed).
- **Harbor oracle (`runs/build-hgte-gated-subband-densities/oracle-r1-hgte-gated-subband-densities`):** reward 1.0, 58 s.

The estimated agent length is 50–80 tool calls: an SdH pipeline from raw data, kdotpy output parsing, ridge tracking and assignments.

## Open risks

- **Experimental truth from a hand-tracked chart.** The ridge guesses (2.81, 0.83, 1.22, 0.45, 0.63, 2.49) are mine, from the 2D chart. Each is the dominant local maximum in its window, and the windows exclude harmonics and sum frequencies. The hole ridge is weak (log power −12.3 against −12.0 for neighbouring harmonics), so an agent without ridge tracking fails dp2 by design.
- **The model check at −0.45 V is not independent of the experiment.** The SdH top surface (0.446) lies inside the model window (0.485 ± 0.04); this is a margin to note.
- **Hole onset rests on an extrapolation.** It extrapolates from −0.55 V to −0.21 V. Variants span 0.085 V; the tolerance is 0.06 V around −0.214.
- **Memorisation.** Fig. S3 of the paper shows the 70 nm chart with model points (low precision). The multiplicity 2 and the BIA-split pockets are in the authors' notebook labels ("doubly-degenerate"), not in the staged files.
- **Result size.** A 26 nm transport block is about 0.45 MB of JSON; `state.json` reached 5.4 MB for the oracle's 22 calls.

## Disk

- **Peak:** about 0.83 GB in `~/.cache/ai-theorist/hgte-gated-subband-densities` (members plus the 36 MB pool), plus the gate's lab image layer (15 MB over the shared runtime image).
- **After the build:** raw downloads, the scratch pool and the gate image were removed. No Harbor `hgte-gated-subband-densities__*` image remained. `fetch.py` and `members.md5` are kept for a refetch.

## Pilot follow-up (2026-10-27): positive-gate assignment (dp9)

In the dev pilots (`runs/dev-hgte-pilot1`), Sonnet 5 and GPT-5.6 Terra missed the same two fields:

| Field | Sonnet | Terra | Reference |
|---|---|---|---|
| top surface at +0.20 V | 1.54 | 1.76 | 2.81 |
| first bulk subband at +0.30 V | 3.84 | 3.17 | 1.22 |

Both took the steep ridge as the bulk subband. For the top surface they took a flat feature near 1.5 (twice the bottom surface, 0.76) or a sum frequency near 1.73. Sonnet justified this with the Hall sum rule: 0.82 + 1.54 + 2.81 = 5.17 = n_tot(+0.20 V). The rule holds for the correct assignment as well, so it cannot tell the two apart.

The raw data rule the alternative out. From prominence-filtered spectra of the archived chart, averaged over ±0.02 V:

- **The steep ridge does not start from zero.** It runs 1.63 (−0.05 V) → 1.97 (0.00) → 2.20 (+0.05) → 2.40 → 2.64 → 2.80 (+0.20) → 2.97 → 3.17 (+0.30) → 3.58 (+0.40). At Vg ≤ 0 it already holds 1.6–2.0, so it cannot be a subband that begins filling there.
- **The newly filled subband does start from zero.** The ridge 0.34 (+0.05) → 0.48 → 0.60 → 1.10 (+0.25) → 1.23 (+0.30) → 1.50 (+0.40) → 1.74 (+0.50) extrapolates to zero near −0.07 V. That is the first bulk subband, and it agrees with the model onset of bands 3/4 (−0.04 to −0.05 V).
- **The ~1.5 feature is flat** at 1.47 / 1.52 / 1.49 for +0.20 / +0.25 / +0.40 V, i.e. twice the bottom surface. It is a harmonic.
- **The model labels agree.** Band 1, which is continuous across neutrality, is the steep branch (2.76 at +0.20 V). Bands 3/4 appear only above n ≈ 0.0027 nm⁻².

Resolution: case (a). The misassignment is legitimate and the evidence rules it out, so the fields stay as they are. Two ablations were added under dp9:

- `steep_ridge_as_bulk_harmonic_as_top`: top 1.420, E1 3.248;
- `steep_ridge_as_bulk_sum_as_top`: top 1.720, E1 3.248.

Both fail by more than 1.3 against tolerances of 0.168 and 0.073. Rubric item 2 now names the onset-from-zero argument.

- **Gate (cpuset 14,15,30,31):** passed. 9 decision points, all 15 ablations fail, all 9 robustness variants pass. The oracle is unchanged (22 of 25 units).
- **Harbor oracle (`runs/build-hgte-gated-subband-densities-fix`):** reward 1.0, 57 s.

The raw members were refetched for this check (checksums verified) and deleted afterwards.
