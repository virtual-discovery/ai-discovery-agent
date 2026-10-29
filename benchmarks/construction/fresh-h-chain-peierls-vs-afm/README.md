# Construction of fresh-h-chain-peierls-vs-afm

First fresh-start (live pw.x) task of the benchmark. See `benchmarks/tools/BUILDING.md`.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --cpuset-cpus 10,11,26,27 --user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none` (done by `qe_batch.py`).

1. `python3 benchmarks/construction/fresh-h-chain-peierls-vs-afm/fetch_pseudo.py` — downloads SSSP 1.3.0 PBE precision (record rcyfm-68h65) into `~/.cache/ai-theorist/fresh-h-chain-peierls-vs-afm/`, checks the record MD5s, extracts `H_ONCV_PBE-1.0.oncvpsp.upf` (md5 1790becc920ee074925cf490c71280fe, sha256 69d63178…5c09) and copies it into `environment/lab/` only.
2. `python3 benchmarks/construction/fresh-h-chain-peierls-vs-afm/run_reference.py` — runs the 90 requests of `reference_plan.py` through the task's own `task_lab.run` (2 MPI ranks) and appends them to `reference_runs.jsonl` (resumable). The committed file was built with `--seed` from an exploration run: 48 + 2 rows whose canonical request was identical were copied (tags in `reference_seeded.json`), 40 were computed directly; the two 4-atom rows come from the worst-case timing runs. Total compute ≈ 2.6 h on 2 cores.
3. `python3 -B benchmarks/construction/fresh-h-chain-peierls-vs-afm/make_task.py` — applies the estimators of `solution/solve.py` to the reference runs and writes `private/hidden.json`, `private/hidden_values.json`, `solution/ablation_values.json` and `design.json` → `calibration`.
4. `.venv/bin/python -m benchmarks.tools.check_task fresh-h-chain-peierls-vs-afm --cpuset 10,11,26,27` writes `check.json`.

`qe_batch.py` is the generic batch runner (list of `{tag, request}` → JSONL of lab results).

## Physics found in construction (2-atom cell, 48 k, MV 0.005 Ry, 10 Å box, 80 Ry)

| R (Å) | NM uniform → AFM gain (meV/atom) | AFM site moment (μB) | FM start | δ_eq | dimerization energy (meV/atom) | δ_c |
|---|---|---|---|---|---|---|
| 1.0 | 0.41 | 0.013 | — | ≈0.26 | ≈200 | — |
| 1.25 | 0.96 (0.31 at 192 k / 0.002) | 0.028 (0.021 converged) | → NM | 0.402 | 646 | ≈0 |
| 1.5 | 5.4 | 0.098 | → NM | — | — | ≈0.045 |
| 1.75 | 35.4 | 0.264 | → NM | 0.572 | 1612 | 0.115 |
| 2.0 | 128.2 | 0.465 | → NM | ≈0.625 | 1929 (vs AFM) | 0.201 |

- A ferromagnetic start relaxes to the nonmagnetic state at every R up to 2.0 Å (no Stoner instability at the DOS minimum of a half-filled 1D band), while the antiparallel start finds the nesting-driven SDW at every R. This is the main discovered-evidence fork.
- The weak SDW at R = 1.25 Å is destroyed by MV 0.02 Ry and by Fermi–Dirac 0.005 Ry, survives Gaussian/MV/FD at ≤0.005 / 0.002 Ry, and its gap needs dense k: 0.288 (48 k), 0.229 (96 k), 0.213 eV (192 k, the maximum the lab allows at this R). Everything gapped (dimers, magnetism at R ≥ 1.75 Å, δ_c) is converged at 48 k / 0.005 Ry to <1% (e.g. δ = 0.11 at R = 1.75: site moment 0.0678 vs 0.0665 at 96 k / 0.002 Ry).
- Cutoff 60 → 80 → 100 Ry changes the 1.75 Å dimerization energy by 1.3 / 0.1 meV/atom; box 8 → 10 → 12 Å by <0.5 meV/atom. Box size is therefore not a decision point (deviation from the candidate design).
- The moment vanishes continuously, with m² linear in δ near δ_c (R = 1.75: m² = 0.00931, 0.00460, 0.00014 at δ = 0.105, 0.110, 0.115; zero at 0.12). The reference δ_c is the m² = A − Bδ² interpolation through the two magnetic points nearest the boundary; starting moments ±0.5 and ±1.0 give the same solution.
- The global minimum is the nonmagnetic H₂-like dimerized chain at every R (AFM starts at δ_eq give zero moments), lower than the best uniform chain by 0.2–1.9 eV/atom.

## Deviations from the candidate design

- **No withheld R and no hidden SCF.** Every graded number is derived from several lab calculations (force-zero minimum, energy differences, boundary interpolation, converged-limit gap); nothing graded is returned by the lab as such.
- **R_c replaced.** "R where the AFM-uniform state undercuts the dimerized NM minimum" does not exist (the dimerized chain wins by ≥0.2 eV/atom everywhere). It is replaced by the coexistence boundary δ_c(R) at 1.75 and 2.0 Å and the ground-state labels at 1.0 and 2.0 Å.
- **"AFM moment at δ = 0"** replaced by the magnetic energy gain (1.75, 2.0 Å) and the converged SDW gap at 1.25 Å (the 1.25 Å magnetic energy gain itself is not k-converged even at 192 k, 0.96 → 0.53 → 0.31 meV, so it is not graded).
- **Box convergence dropped** as a decision point (see above); the 4-atom cell is kept only as an allowed alternative (↑↑↓↓ or period-4 patterns).
- Energies are returned per cell in Ry, as pw.x prints them; normalization is a (weak, convention) decision point.

## Calibration (full table: `design.json` → `calibration`)

| field | reference | tolerance | robustness (must pass) | nearest wrong path (must fail) |
|---|---|---|---|---|
| δ_eq 1.25 / 1.75 | 0.4021 / 0.5717 | ±0.01 | energy-parabola vertex 0.4020 / 0.5715 | 0.1-step grid minimum 0.40 / 0.60 |
| dimerization energy 1.25 / 1.75 | 646.0 / 1611.7 | ±10 meV | 60 Ry: 1610.4; converged uniform reference 645.5 | NM uniform reference 1647.3 (+35.6); grid minimum 1588.3 (−23.4) |
| magnetic gain 1.75 / 2.0 | 35.44 / 128.2 (96 k, 0.002 Ry) | 10% | standard settings 35.65 / 128.3 | FM start 0.03 / 0.04; per cell ×2 |
| δ_c 1.75 / 2.0 | 0.1151 / 0.2011 | ±0.015 | bisection midpoint 0.1175 / 0.2025; converged settings 0.1144 / 0.2011 | 0.1-grid midpoint 0.15 / 0.25; linear moment extrapolation 0.180 / 0.213; FM start 0 |
| gap 1.25 | 0.2131 eV (192 k, MV 0.002) | ±0.04 eV | 96 k MV 0.002: 0.2294; Gaussian 96 k: 0.2265 | 48 k: 0.288; MV 0.02: 0.0008; FD 0.005: 0.002; FM start: 0 |
| ground state 1.0 / 2.0 | dimerized_nonmagnetic ×2 | exact | — | AFM wins at 2.0 |

Tightest margins: gap (ablation 48 k is 0.075 eV from the reference against a 0.04 tolerance; robustness 96 k is 0.016 inside), δ_c at 2.0 (linear extrapolation 0.012 off, passes on its own; the ablation fails through 1.75), dimerization energy (60 Ry is 1.3 meV off; nearest wrong path 23 meV).

## Oracle

`solution/solve.py`: 25 evaluations, 29 of 36 units (81%): uniform NM/AFM at 1.25, 1.75, 2.0 Å; force-zero δ_eq scans (3 points) at 1.25 and 1.75 Å; the 1.25 Å SDW gap at 96 and 192 k (MV 0.002); boundary searches on a 0.01 grid at 1.75 (3 runs) and 2.0 Å (4 runs); AFM starts at the dimerized minima of 1.0 and 2.0 Å. Sum of reference SCF wall times for these requests: 44 min. An offline replay against `reference_runs.jsonl` reproduces the references (errors ≤ 0.6% / 0.0005).

## Gate and Harbor

- **Gate** (full run 2026-10-19T21:13Z, 2719 s container time; `check.json` was then regraded with `--reuse-state` on 2026-10-20 after a `design.json`-only edit, which is why it shows 0 container seconds): passed. 8 decision points, 11 ablations all failing, 9 robustness variations and 4-significant-digit rounding passing, 9 hidden values, task 0.17 MB. The live oracle (25 evaluations, 29/36 units, no failures, 45 min in the container) reproduced the references exactly where it repeated a reference request (run-to-run differences 0) and to ≤0.6% / 0.0005 elsewhere.
- **Harbor oracle** (`runs/build-fresh-h-chain-peierls-vs-afm/oracle-r1-fresh-h-chain-peierls-vs-afm`, 2026-10-20): reward 1.0; 50 min wall, 25 evaluations, 29 units, 47.8 min of SCF time on the cores pinned by the coordinator's daemon.
- Peak disk: 61 MB of downloads (deleted after the gate; `fetch_pseudo.py` refetches them) + ~4 MB of job folders under `runs/`.

## Notes and open risks

- Reference rows were computed before two cosmetic parser changes: `total/absolute_magnetization` are now read from the XML at full precision (they were the 2-decimal stdout values) and the request ranges in `public.json` were narrowed afterwards (box ≤ 10 Å, ecutwfc ≤ 80 Ry; the box-12 and 100-Ry rows serve only as convergence evidence). No estimator uses those fields; energies, forces, site moments and gaps are parsed identically.
- Timing: the heaviest allowed request (k-points × cell length = 480 Å, e.g. 2 atoms at 1.25 Å with 192 k or 4 atoms with 96 k) takes ≈6 min; close to the magnetic boundary the SCF needs 20–30 iterations (up to ~4 min at 48 k, ~5.5 min at 96 k for R = 2.0 Å). The hard pw.x timeout is 20 min; a near-critical run at the largest k could exceed 10 minutes.
- The H-chain physics is textbook (Peierls, Hubbard), so the main difficulty is in the discovered forks: FM start → NM, smearing/k-sensitivity of the weak SDW, and the budgeted boundary search. δ_eq ≈ 1 − 0.75/R can be guessed from the H₂ bond length.
- Run-to-run noise: identical requests reproduce energies to 1e-8 Ry and site moments to 1e-4 μB (seeded vs recomputed rows, e.g. R1.25_d0_afm); all tolerances are orders of magnitude larger.
