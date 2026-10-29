# Construction of mnwo4-halfmoon-frustration

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). Source: N. Natsch, T. N. Tošić, J.-R. Soh, N. A. Spaldin, Materials Cloud record p8297-g9n39 (DOI 10.24435/materialscloud:s5-a8, CC BY 4.0), paper arXiv:2503.23704. The record holds the authors' UppASD inputs and outputs only; the experimental neutron data are not archived, so the authors' model is the system under study ("teacher") and its archived spin-dynamics outputs are the data.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) ...`), with this folder mounted at its repository path so that `make_task.py` finds `../../tasks/<id>`.

1. `python3 -B fetch.py CACHE/raw` downloads `Fig_2.zip`, `Fig_4.zip`, `Fig_7.zip` and `Fig_5_and_Fig_10.zip` (113 MB) through the record's per-file API, checks the record md5s and unpacks only the members used (126 MB).
2. `python3 -B extract.py CACHE` writes `CACHE/extracted.npz` (21 MB) and `CACHE/model.json`. It asserts identical Hamiltonians across runs, that the archived q files reproduce the outputs' Cartesian q, and the map plane q = (−l/2, k, l). Re-running fetch + extract reproduced the arrays bit for bit.
3. `python3 -B calibrate.py CACHE` (~2 min) writes `calibration_fits.json`: the oracle path, defensible variants and the wrong paths, each with every graded quantity.
4. `python3 -B make_task.py CACHE` writes `data.npz`, `data_meta.json`, `catalog.json`, `private/hidden.json`, `hidden_values.json`, `solution/ablation_values.json`, both `public.json` files and `calibration_table.json`, and prints the table.
5. Optional: `OMP_NUM_THREADS=1 python3 -B nonrecall.py CACHE SEED` (~13 min per seed) reruns the no-literature route; its outputs `nonrecall_seed{1..4}.json` feed `calibrate.py`.
6. Optional: `python3 -B mc_check.py CACHE WORKDIR` (~20 min, compiles `mc.c`) reproduces the Monte Carlo cross-check quoted in SOURCE.md.
7. `.venv/bin/python -m benchmarks.tools.check_task mnwo4-halfmoon-frustration --cpuset 6,7,22,23` writes `check.json`.

Supporting files: `model.py` (J(q), Luttinger–Tisza, AF1 energies, 8-sublattice AF1 linear spin waves with anisotropy, SCGA with anisotropy, block grids), `fitlib.py` (spectral variable-projection fit, null-direction scan, predictions). `solution/solve.py` embeds the same code.

## What the archive showed (design inputs)

- **Conventions (verified from raw outputs).** The archived 1.5 K exchange energy −0.13974 mRy/atom equals −(1/N) Σ_i Σ_j J_ij ê_i·ê_j over all `jfile` rows for the archived final spin configuration (ordered pairs, unit vectors); the anisotropy energy equals K1⟨(ê·n)²⟩ with K1 = −1.41×10⁻² mRy (easy axis). The `jfile` (mRy) differs slightly from the paper's rounded Table I (J5 −0.0148 vs −0.0176, J8 −0.0127 vs −0.016, J12 −0.0215 vs −0.0256 meV).
- **Spectra.** AF1 linear spin waves with S = μ/(gμB) = 5/2 (precession of 5 μB classical moments) reproduce the archived 1.5 K S(q,ω) peaks (gap 0.574 meV at k_c, top 2.58 meV); the archived energy axis is correct.
- **dp1, exact null vector.** The 808 branch energies along the archived [102] + [010] path have rank 12 of 13 parameters; the null vector is J3 − J11 (the two bonds swap parallel/antiparallel roles in AF1 and their phases are complex conjugates on both path segments). Along that line the spectra are unchanged but the ground state is not: the min-norm split puts the flat spiral 0.155 meV/Mn below AF1 and moves k_LT from (−0.228, ½, 0.441) to (−0.179, ½, 0.416). The 20 K map χ² has its minimum near the archived split (J3 ≈ −0.28 vs −0.270).
- **SCGA vs Monte Carlo.** Equilibrium MC of the archived Hamiltonian matches the archived maps no better than SCGA (7×8-block Pearson 0.923/0.869/0.701 vs 0.933/0.889/0.725 at 14/20/30 K; MC vs SCGA 0.99). A map-only SCGA fit does not fix the absolute exchange scale (it drifted to J's 1.7× too large, or ran away), and a joint fit that weights the maps heavily lands on a wrong point (spectral cost 536 vs 342 for the archived J on 51 spectra). The oracle therefore lets the spectra fix everything except the null direction, which the maps fix.
- **Box quantisation (dp5).** The archived 11 K full-zone maxima are (−8/36, ½, 16/36) and (−11/48, ½, 11/24): commensurate approximants of the two boxes; the thermodynamic-limit Luttinger–Tisza vector is (−0.2277, ½, 0.4413).
- **Mean field.** T_MF = (2/3)λ_max/k_B = 16.1 K (easy-axis component 17.6 K); SCGA T_c ≈ 11.2 K, close to the paper's heat-capacity ~11 K. Reduced model 5e has T_MF = 19.5 K (easy axis ~21 K), so RPA is singular or near-singular at 20 K.
- **Noise.** The archived maps are single runs; mirror pairs (k → −k) give block noise σ ≈ 7×10⁻⁶, 3.3×10⁻⁶, 2.8×10⁻⁶, 2.2×10⁻⁶ at 14/20/30/60 K, so even a perfect model reaches only r ≈ 0.97 (14 K) on 7×8 blocks; the 60 K map is noise-dominated.

## Deviations from the design (candidates.json / DETAILS.md section 14)

- **Lab serves archived outputs, not a calculator.** LSWT, SCGA and ground-state energies are a few hundred lines of numpy, so a calculator lab would not bind the budget. The agent buys archived spectra, maps, reduced-model maps and zone slices (catalogue 169 units, budget 33).
- **No synthetic data.** The design proposed LSWT-generated dispersions with added noise; instead the archived spin-dynamics S(q,ω) (with its own thermal broadening and statistical noise) is sold as is.
- **Graded targets.** (1) The withheld 14 K map is kept (Pearson ≥ 0.85 on 27×29 block averages instead of > 0.8 on raw pixels). (2) "Identity of the 3 dominant AFM bonds" is replaced by predictions of the withheld reduced-model maps 5d (no J9) and 5e (no J3/J7/J10/J11) plus a categorical origin; a |J| ranking is unstable (J4 −0.334 vs J9 −0.317). (3) Model T_N is not graded: no temperature sweep is archived and the paper's ~11 K comes from unarchived heat-capacity runs; it is a rubric item. (4) The withheld zone-boundary magnon energy is dropped (on-path points are interpolation; off-path points would be builder forward modelling); J, D and the spiral–AF1 energy balance (the design's second question) are requested and reviewed but not graded (see the coordinator review below). Added and graded: the Luttinger–Tisza propagation vector, computed from the archived `jfile`.
- **Polarisation channel (design dp6) became a form-factor decision.** The archived |S| has neither polarisation nor form factor, and in the nearly isotropic paramagnet the perpendicular projection changes nothing (r 0.929 vs 0.928), so it cannot be an ablation; applying the form factor can.

## Coordinator review: grading without literature recall (2026-10-20)

The first version graded J1, J6, J9, D and the spiral − AF1 energy, and the oracle started its spectral fit from the paper's rounded Table I. Question: can an agent without those numbers reach the windows with compute that fits the agent container? Tested routes (single-threaded, numpy/scipy):

- **Generic multistart, spectra only:** 20 random starts with sigma continuation 0.12 → 0.03, and 12 more with 0.3 → 0.03. All landed in wrong minima (spectral cost 118–192 vs 121.5 at the archived parameters).
- **Map-informed start:** an SCGA fit of the 20/30 K maps from a flat start gives sensible ratios for the strong couplings but no absolute scale. A scale scan plus spectral refinement then failed (J1 −0.29).
- **DE on spectra only (10 spectra, ~1.5 min per seed):** seed 2 found a point that fits the 20 spectra *better* than the archived parameters (cost 109.5 vs 121.5; J1 = +0.02, D = −0.37). Its maps are poor (map χ² 1246 vs 421). The archived noisy spectra alone do not identify the Hamiltonian.
- **Joint DE (spectra + maps with a free SCGA temperature scale + the 11 K ordering-peak constraint; `nonrecall.py`, `fitlib.Joint`; ~13 min per seed):**

  | Seed | J1 | D | Spiral − AF1 | k distance | r14 | r5d | r5e | Outcome |
  |---|---|---|---|---|---|---|---|---|
  | 1 | −0.62 | −0.11 | −0.020 | 0.0159 | 0.935 | 0.735 | 0.959 | passes |
  | 2 (lowest objective) | −0.37 | −0.18 | +0.048 | 0.0156 | 0.917 | 0.721 | 0.910 | passes |
  | 3 | −0.62 | −0.17 | +0.033 | 0.0076 | 0.938 | 0.600 | 0.602 | fails 5d, 5e |
  | 4 | −0.87 | −0.10 | −0.076 | 0.0214 | 0.922 | 0.705 | 0.954 | fails k |

**Conclusion (option b):** a feasible no-literature route reproduces the forward-modelling targets: the propagation vector, the three maps and the origin. It does not pin the individual J values, D or the sign of the spiral − AF1 balance. Those fields are still required in the submission and reviewed with the report, but they are no longer graded.

The k tolerance was widened from 0.013 to 0.020 to accept the no-literature fits. At that tolerance the box-peak / experimental-k decision point (old dp5) no longer fails (k_ic 0.0157, 48-box peak 0.017), so it was dropped from grading and kept as rubric item 7. The 14 K threshold was set to 0.85: defensible fits score ≥ 0.915, and the nearest wrong path (temperature-insensitive) scores 0.825.

The oracle still starts from Table I, as a fast deterministic reference. The best-of-4-seeds no-literature solution is listed in ROBUSTNESS and passes every graded check.

## Decision points (7)

| dp | Fork | Ablation result |
|---|---|---|
| dp1 | Dispersion only, min-norm split of J3 − J11 | k off 0.049 |
| dp2 | Spin length 1 in the dynamics | r14 0.72, r5e 0.81, k off 0.031 |
| dp3 | Mean-field RPA | r14 0.50 (0.81 at 17 K), r5e 0.31 |
| dp4 | Isotropic model | k off 0.034 (maps pass) |
| dp5 | Magnetic form factor applied | r14 0.71, r5d 0.29, r5e 0.75 |
| dp6 | Temperature-insensitive (spin-liquid) half-moons | r14 0.83, r5d 0.61, r5e 0.79, origin wrong |
| dp7 | Quasi-1D chain (J1, J4) | all maps 0.35–0.51, k wrong, origin wrong |

## Calibration

Full table: `design.json` → `calibration` and `calibration_table.json`.

| Field | Reference | Oracle (calibrate / in lab) | Defensible spread | Nearest wrong path | Tolerance |
|---|---|---|---|---|---|
| k distance (r.l.u.) | 0 | 0.0109 / 0.0045 | 0.004–0.0159 (no-literature best seeds 0.0156, 0.0159) | S = 1 0.031; isotropic 0.034; min-norm 0.049 | ≤ 0.020 |
| r 14 K map | 0.943 | 0.916 / 0.936 | 0.915–0.938 (MC ≈ 0.92) | T-insensitive 0.825; RPA 17 K 0.815 | ≥ 0.85 |
| r 5d | 0.744 | 0.736 / 0.740 | 0.721–0.746 | T-insensitive 0.612 | ≥ 0.68 |
| r 5e | 0.941 | 0.937 / 0.942 | 0.910–0.959 | S = 1 0.811; T-insensitive 0.788 | ≥ 0.88 |

Reported but not graded: J1 −0.587 (oracle −0.598; no-literature seeds −0.37 to −0.87), D −0.192 (−0.10 to −0.18), spiral − AF1 +0.031 meV/Mn (−0.076 to +0.048).

Tightest margins:
- **k:** the no-literature best seed (0.0156) passes by 0.0044; spin length 1 fails by 0.011.
- **5d:** the no-literature best seed (0.721) passes by 0.041; temperature-insensitive fails by 0.068.
- **5e:** the no-literature best seed (0.910) passes by 0.030.
- **14 K map:** temperature-insensitive fails by 0.025.

## Oracle

Buys 20 spectra (indices 0–100, step ~5), the 20 K and 30 K maps and the 36×24×36 zone slice at k = ½: 23 evaluations, 29 of 33 units (88 %). Fits the spectra from the paper's rounded Table I (sigma continuation 0.06 → 0.03 meV), scans the J3 − J11 null direction against the block-averaged maps, then predicts with SCGA, Luttinger–Tisza and the AF1/flat-spiral energies. Container time 25 s.

- **Gate (2026-10-20 11:15Z, after the review changes):** passed. 7 decision points; all 8 ablations fail; 5 robustness variants (including the no-literature route) and 4-significant-digit rounding pass.
- **Harbor oracle:** reward 1.0 (`runs/build-mnwo4-halfmoon-frustration/oracle-r1-*/*/verifier/reward.txt`, 1 min 4 s).

## Open risks

- **Memorisation** now affects only the reported, ungraded J and D. An agent that recalls Table I still has to forward-model k, the three maps and the origin, and recall helps it pass them.
- **Without literature values the global search is the hard part.** 2 of 4 joint-DE seeds pass. Choosing the lowest-objective seed passes, but only if the agent weights the spectra and the maps jointly. Spectra-only fits can beat the archived parameters on the spectra and still fail the maps.
- **5d map discriminates weakly:** the archived model itself scores 0.744.
- **SCGA vs spin-dynamics differences.** All calibration predictions use SCGA. MC was checked only with the archived Hamiltonian, where it scores about 0.01 lower on the 14 K map.
- **Peak disk use:** about 0.71 GB in the cache during the build, plus 0.25 GB for the review refetch; raw downloads deleted afterwards.
