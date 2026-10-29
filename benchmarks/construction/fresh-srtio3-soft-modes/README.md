# Construction of fresh-srtio3-soft-modes

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7 `pw.x`
(and `ph.x` electric-field) calculations for perovskite SrTiO3 that the agent requests; the hidden truth comes from
reference calculations made through the same `task_lab.run` code path.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --network none --cpuset-cpus
10,11,26,27 --user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none` and `mpirun -np 2` (done by `qe_batch.py`).

1. `python3 -B fetch_pseudos.py --install`: downloads SSSP 1.3.0 PBE precision (record rcyfm-68h65) into
   `~/.cache/ai-theorist/fresh-srtio3-soft-modes/`, checks the record MD5s, extracts `Sr_pbe_v1.uspp.F.UPF`,
   `ti_pbe_v1.4.uspp.F.UPF` and `O.pbe-n-kjpaw_psl.0.1.UPF`, checks them against the SSSP metadata MD5 and the SHA256
   in `pseudo_checksums.json`, and copies them plus `SSSP_LICENSE.txt` into `environment/lab/`.
2. `python3 -B run_reference.py --stage A`, then `--stage B` (needs stage A). Requests are built in `reference_plan.py`
   with the oracle's own request builders; results append to `reference_runs.jsonl` (89 runs, resumable; 2.7 h of
   pw.x/ph.x on 2 cores). Job folders (pw.in/pw.out, ph.out) are in `runs/build-fresh-srtio3-soft-modes/reference_jobs/`.
3. `.venv/bin/python -B make_task.py`: applies the estimators of `solution/solve.py` and writes `private/hidden.json`,
   `private/hidden_values.json`, `solution/ablation_values.json`, `calibration.json` and `design.json` -> `calibration`.
4. `.venv/bin/python -B replay.py`: offline replay of the oracle against `reference_runs.jsonl` (seconds, no QE),
   grading ablations, robustness variants and rounding as check_task does.
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-srtio3-soft-modes --cpuset 10,11,26,27` writes `check.json`.

Exploration runs (early lab versions incl. ph.x phonons at Gamma and R) are in
`runs/build-fresh-srtio3-soft-modes/jobs/`; no graded number comes from them.

## Lab

- **Request (no defaults):** kind (`scf` | `dielectric`), cell (`perovskite_5` | `rotated_10`, the R-point doubled
  fcc-type cell with two formula units), a, c (|c/a - 1| <= 5%), one Cartesian displacement per atom (<= 0.4 A),
  ecutwfc (45/60/75/90 Ry, dual 8), Gamma-centred mesh (4/6/8 for 5 atoms, 2/3/4/6 for 10 atoms).
- **Fixed:** PBE, SSSP precision Sr/Ti GBRV ultrasoft (semicore), O PSlibrary PAW; fixed occupations;
  conv_thr 2e-12 Ry per atom; ph.x tr2_ph 1e-16, epsil only (trans false).
- **Returned:** energy, forces, stress, positions, masses, pw.x k count / symmetry; for `dielectric` the eps_inf tensor
  and raw Born charges (no sum rule). pw.in/pw.out/ph.out kept; scratch deleted.
- **Charge = canonical size:** scf units = ceil(N_k x weight(cell) x cutoff factor) with weights 0.03 (5 atoms) and
  0.15 (10 atoms) per symmetry-distinct k point and cutoff factors 0.7/1/1.4/2; N_k is counted by `task_lab` from the
  space group of the displaced cell (matches pw.x's count on all 89 reference runs); dielectric 4/5/9 units at
  k4/k6/k8 (scaled by cutoff). Max 10 units per request; budget 80. About 1 unit per 40-60 s on 2 loaded cores.
- **Measured wall times (2 cores, host load 15-20):** cubic 5-atom scf 30-60 s; displaced 5-atom k6 70-120 s, k8
  120-240 s; dielectric k4 3.3 min, k6 5.1 min, k8 8.5 min (the slowest offered request); 10-atom k4 2.1-4.6 min,
  k6 3.6-4.2 min.

## Physics found in construction (PBE, 60 Ry unless stated)

- **a0:** 3.9359 A (cubic fit of E(a), 60 and 75 Ry agree to 0.0001 A); zero pressure 3.9338 (60 Ry) / 3.9346 (75 Ry).
- **Zone-centre polar block (0.01 A frozen phonons, ASR, z block 4x4):** soft TO -192.0 (k4), -143.4 (k6), -138.4 (k8)
  cm^-1; other TO 144/509 cm^-1. 0.04 A displacements: -133.4 (k6). At a = 3.905 A: -92.0 (k6). Ti-only on-site
  force constant gives +53 (k6), +60 (k8), -73 (k4) cm^-1. Soft eigendisplacement (Ti = 1): O_z -2.23, O_x/O_y -1.99.
  The DFPT phonon run in exploration (k4) gave -197 cm^-1 and Z*/eps identical to the k4 dielectric run.
- **Dielectric:** eps_inf 6.849 / 6.349 / 6.294 and Z*(Ti) 7.811 / 7.386 / 7.346 (k4/6/8); raw charge-neutrality
  violation 0.265 / 0.036 / 0.017 e. Highest LO (nonanalytic term on the full block): 731.5 cm^-1 (k8/k8).
- **Polar well along the soft eigendisplacement (fixed cubic a0):** depth 3.33 (k8), 3.67 (k6), 8.74 (k4 pattern and
  mesh) meV/f.u. at Ti ~0.028 A; amplitudes 0.04/0.08/0.12 alone give 0.98 meV/f.u.
- **a0a0c- tilt (rotated_10, k4, 8 angles, same-cell reference):** 4.67 deg, depth 4.03 meV/f.u., harmonic R mode
  -73.0 cm^-1. 6x6x6 changes the 4-deg energy by 0.11 meV/f.u. The undistorted 10-atom k4 energy is 0.70 meV/f.u. above
  the 5-atom k6 energy (16 meV/f.u. below the 5-atom k4 energy). Exploration DFPT at R (k4): -74.7 cm^-1.
- **Competition:** tilt 4.03 > polar 3.33 meV/f.u. (k6 for both: ~3.9 vs 3.67); k4 force constants flip it.
- **Biaxial +1%:** c = 3.9122 A (75 Ry energy minimum; 60 Ry 3.9121; stress zero 3.9084 at 60 Ry, 3.9097 at 75 Ry,
  a Pulay offset inside tolerance); in-plane polar mode -187.2 (k8), -191.0 (k6); with c = a0, -192.4 (k6).

## Calibration (full table: `calibration.json`, `design.json` -> `calibration`)

| field | reference | tolerance | max robust error | nearest failing ablation (error) |
|---|---|---|---|---|
| a0_A | 3.9360 | 0.004 | 0.0014 (zero pressure 75 Ry) | experimental 3.905 (-0.031) |
| gamma_soft_mode | -138.4 | 10 | 5.0 (k6; 0.04 A) | experimental lattice (+46); k4 (-54); Ti-only (+192) |
| gamma_LO_highest | 731.5 | 4 | 0.66 | k4 force constants+tensors (+3.4, also fails soft mode/eps/Z); k4 tensors (+12); analytic (-222) |
| eps_inf | 6.294 | 0.12 | 0.054 (k6) | k4 (+0.555) |
| born_charge_Ti | 7.346 | 0.10 | 0.041 (k6) | k4 (+0.465) |
| polar_well_depth | 3.33 | 0.5 | 0.34 (oracle k6) | coarse amplitudes (-2.34); k4 (+5.4); Ti-only (-3.3) |
| R_mode | -73.0 | 8 | 0.21 | per 10-atom cell (-30) |
| tilt_angle_deg | 4.674 | 0.3 | 0.006 | (no ablation targets it) |
| tilt_well_depth | 4.03 | 0.5 | 0.004 | cross-cell reference (-0.70); per cell (+4.0) |
| strained_c_A | 3.9122 | 0.005 | 0.0038 (stress zero 60 Ry) | c = a0 (+0.024) |
| strained_inplane | -187.2 | 10 | 3.8 (k6) | cubic soft mode (+49); c = a0 (-5.2, passes; fails via c) |
| larger_energy_gain | tilt | exact | - | k4 path: polar |

Tightest margins: tilt depth cross-cell ablation fails by 0.20 meV/f.u.; strained c stress-zero robustness at 60 Ry
uses 0.0038 of 0.005 A; oracle polar depth (k6) uses 0.34 of 0.5 meV/f.u.; k4 LO alone (3.4 of 4) passes LO but the same
path fails four other fields.

## Oracle

`solution/solve.py`: EOS (4 scf), z-polar frozen phonons at k6 and k8 (6), dielectric k6 and k8 (2), cubic
reference + 4-point polar scan (5), tilt scan 0/2/3/4/5 deg in rotated_10 k4 (5), strained c scan (3), strained
in-plane frozen phonons (4): 29 evaluations, 66 of 80 units (82%), 62 min of summed reference wall time. Replay
reproduces the references exactly.

## Deviations from the candidate design

- **No DFPT phonons.** ph.x phonons took 13 min (Gamma) and 10 min (R) at k4 on 2 cores, and k4 is unconverged for the
  soft mode; the lab offers pw.x frozen phonons and a `dielectric` ph.x kind (eps_inf, Z*) only. The design's
  "frozen-phonon amplitude vs DFPT" fork became dp8 (amplitude sampling of the anharmonic polar well).
- **Cutoff fork dropped:** 60 vs 75 Ry changes a0 by 0.0001 A and strained c by 0.0001 A; the Pulay stress shifts the
  stress-zero c by 0.004 A (inside tolerance). Cutoffs are still a budget choice.
- **ASR fork not graded separately:** the frozen-phonon Sr column from the sum rule and from a direct Sr
  displacement agree to 1e-4; the ASR remains part of the rubric. Pseudopotential semicore not offered as a choice
  (SSSP files include semicore states).
- **No hidden SCFs:** graded quantities are the agent's own calculations at specified configurations; the polar well
  is defined along the harmonic soft-mode eigendisplacement so that it is unique.
- **Strain:** the in-plane frequency is nearly independent of the Poisson relaxation, so dp7 is graded via strained_c_A.
- **Competition** graded as the verdict 'larger_energy_gain' (tilt); the biquadratic coupling is a rubric item.

## Open risks

- **Recall:** SrTiO3 soft modes are famous, but all graded numbers are protocol-specific (PBE a0, k-converged
  frozen phonons); the experimental lattice gives -92 cm^-1 and fails.
- **Near-tie:** the polar (3.3-3.7) and tilt (3.9-4.0) depths differ by ~0.3-0.7 meV/f.u. on converged meshes; the
  verdict is unique on 6x6x6 and 8x8x8 for both, and flips only with 4x4x4 polar data.
- **Mesh convergence:** 8x8x8 is the densest offered 5-atom mesh; soft-mode and polar-well differences 6->8 (5 cm^-1,
  0.34 meV/f.u.) are inside tolerance, but the truly converged values may differ slightly from k8.
- **Timing:** the dielectric k8 request takes 8.5 min under load; a slower host could exceed 10 min (timeout 40 min).

## Gate and Harbor

- **Gate** (`check_task --cpuset 10,11,26,27`, 2026-10-27 ~03:40Z, final package): passed, no issues. 8 decision
  points, 11 ablations all failing, 11 robustness variants and 4-significant-digit rounding passing, deterministic,
  11 hidden values, task 2.1 MB. Live oracle: 29 evaluations, 66 of 80 units, 0 failed, 3629 s container time.
- **Harbor oracle** (`runs/build-fresh-srtio3-soft-modes/oracle-r1-fresh-srtio3-soft-modes`): reward 1.0, 59.4 min
  total runtime with the coordinator's CPU-pinning daemon. Per-trial images removed.
- **Disk:** peak ~70 MB (63 MB SSSP archive in `~/.cache/ai-theorist/fresh-srtio3-soft-modes/`, deleted after the
  gate; `fetch_pseudos.py` refetches it); 5 MB of job folders under `runs/build-fresh-srtio3-soft-modes/`.
