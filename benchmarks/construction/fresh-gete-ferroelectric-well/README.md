# Construction of fresh-gete-ferroelectric-well

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7 `pw.x`
SCF calculations of the two-atom rhombohedral GeTe cell that the agent requests; the hidden truth comes from
reference calculations made through the same `task_lab.run` code path.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --network none --cpuset-cpus 8,9,24,25 --user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none`, `mpirun -np 2` (done by `qe_batch.py`).

1. `python3 -B fetch_pseudos.py --install`: downloads SSSP 1.3.0 PBE precision (Materials Cloud rcyfm-68h65; 63 MB)
   and four PSlibrary 1.0.0 files from the QE library into `~/.cache/ai-theorist/fresh-gete-ferroelectric-well/`,
   checks every MD5/SHA256 against `pseudo_checksums.json` and the SSSP metadata, and installs the six UPF files and
   `SSSP_LICENSE.txt` into `environment/lab/`.
2. `python3 -B run_reference.py --stage A` (67 runs: rocksalt EOS on k12/16/20, cutoff/k/smearing ladders, a
   10-point R3m stencil on k8/12/16, the experimental cell), then `python3 -B make_task.py --stage-a` (writes
   `reference_stage_a.json`: rocksalt a0 and the oracle's rounded R3m model minimum 4.375 Å / 57.80° / 0.4695), then
   `python3 -B run_reference.py --stage B` (52 runs). Requests are in `reference_plan.py`; results go to
   `reference_runs.jsonl` (resumable). Total ≈ 6 h on 2 cores (host load 5–10).
3. `python3 -B make_task.py`: applies the estimators of `solution/solve.py` and writes `private/hidden.json`,
   `private/hidden_values.json`, `solution/ablation_values.json`, `calibration.json` and `design.json` → `calibration`.
4. `.venv/bin/python -B replay.py`: offline replay (seconds, no QE) of `solve.py` against `reference_runs.jsonl` with a
   stand-in lab client (requests matched on the canonical request, charged by `task_lab.validate`), then every
   ablation (must fail), robustness variant and 4-digit rounding (must pass).
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-gete-ferroelectric-well --cpuset 8,9,24,25`.

Exploration runs (an earlier draft of the same lab, including an ibrav=5 cell orientation and relax/vc-relax kinds
that were later removed) are in `runs/build-fresh-gete-ferroelectric-well/explore_runs.jsonl`; no graded number
comes from them.

## Lab

- **One request schema, no defaults:** `kind` (`scf` only), `pseudopotential_set` (`sssp` = protocol; `pslibrary`;
  `pslibrary_fullrel` = spin–orbit), `a_rh_A`, `alpha_deg`, `x_Te`, `ecutwfc_Ry` (40/50/60/75, dual 8), `kmesh`
  (6–24), `kshift`, `occupations` (fixed | MV smearing) and `degauss_Ry`.
- **Cell:** ibrav = 0 with a_i = a (p, q, q) cyclic, threefold axis along Cartesian [111]; at α = 60° the cube axes are
  x, y, z, so pw.x finds all 48 operations of rocksalt. (With ibrav = 5 pw.x found only 12 and used 2.6× more k points.)
- **Returned:** energy, stress (pw.x sign), pressure, forces, final structure echo, band edges on the mesh with insulator
  filling (valence top, conduction bottom, indirect and smallest direct gap with k positions), Fermi energy and −TS
  when smeared, k-point and symmetry counts, wall time. `pw.in`/`pw.out` kept; scratch deleted, also on failure.
- **Charge = canonical size:** units = ceil(N_k × set factor (1 / 0.8 / 3.5) × ecut factor (0.75/1/1.35/2) / 160),
  N_k = pw.x's symmetry-distinct k points (reproduced from kpoint_grid + irreducible_BZ; matched pw.x on all 119
  reference runs). ≈ 1 unit per minute. Cap 8 units per request (≤ ~9 min: polar k20 = 5 units, 290–520 s;
  spin–orbit R3m k12 = 5 units, 370 s). Budget 70.

## Physics found in construction (SSSP, PBE, no SOC, 50 Ry)

- **Rocksalt:** a0 = 6.0093 Å (k20 six-point P(a)); direct gap at L 0.35 eV. Its energy converges slowly with k:
  relative to k24, +1.4 (k12), +0.4 (k16), +0.1 meV (k20) at a0; shifted k12 is close to k24 (contains no L).
  Cutoff: pressures at 40/50/60/75 Ry agree within 0.04 GPa (no Pulay fork; not a decision point).
- **R3m (P = 0):** a_rh = 4.3775 Å, α = 57.743°, x = 0.4691, V = 3.6 % above rocksalt; gap 0.56 eV. Energy
  0.55 meV too low at k12, converged to 0.01 meV at k16. Well depth 31.86 meV/f.u.
- **Clamped rocksalt cell:** depth 18.43 meV/f.u. at x ≈ 0.481 → shear + volume coupling supply 42 % of the depth.
- **Under compression:** rocksalt gap at L 0.15 eV at 3 GPa, closed (overlap) by 6.6 GPa; the harmonic polar force
  constant changes sign at 8.58 GPa (k24, MV 0.002; k20 MV 0.005: 8.66; k20 MV 0.002: 8.87). Fixed occupations
  (k16) give 9.89, k12/MV 0.01 gives 9.65. At x = 0.48 the force includes the quartic term and the apparent force
  constant is positive from 7.4 GPa up (linear extrapolation to zero: 1.5 GPa).
- **Spin–orbit (PSlibrary pair, fixed structures):** +1.46 meV/f.u. at k16 (+1.86 at k12). The SSSP and PSlibrary
  scalar depths differ by about the same amount, so comparing the SOC result with SSSP gives +0.06.

## Calibration (full table: `design.json` → `calibration`, `calibration.json`)

| field | reference | tolerance | max robust error | nearest failing wrong path |
|---|---|---|---|---|
| cubic_a0_A | 6.0093 | ±0.008 | 0.0007 (oracle 3-point k16) | (no ablation; graded for the volume reference) |
| r3m_a_rh_A | 4.3775 | ±0.005 | 0.0025 (k16 10-point model) | one-sided model +0.0076; experimental cell −0.098 |
| r3m_alpha_deg | 57.743 | ±0.1 | 0.049 | one-sided model −0.152; experimental cell +0.62 |
| r3m_x_Te | 0.4691 | ±0.0008 | 0.0004 | one-sided model −0.0010; experimental cell +0.0068 |
| well_depth (meV/f.u.) | 31.86 | 4 % | 0.4 % (k20 both phases) | k12 +6.1 %; cubic at R3m volume +31 %; clamped −42 %; per atom −50 % |
| clamped_well_depth | 18.43 | 5 % | 2.5 % (oracle: 2 points + cubic) | k12 +8.6 %; per atom −50 % |
| soft_mode_critical_pressure_GPa | 8.578 | ±0.6 | 0.29 (k20 MV 0.002) | k12/MV 0.01 +1.07; fixed occupations +1.31; x = 0.48 −7.1 |
| soc_well_depth_change | 1.460 | ±0.6 | 0.40 (k12, the oracle route) | mixed with SSSP −1.40 |

Tightest margins: α for the one-sided model (fails by 0.05°; it also fails a_rh and x), SOC at k12 (passes by 0.2 meV),
soft-mode k12/MV 0.01 (fails by 0.47 GPa), depth at k12 (fails by 2.1 %).

## Oracle

`solution/solve.py`: 3 rocksalt pressures (k16) → a0; centred 7-point R3m stencil (k12, energies + full gradients) →
quadratic model minimum (rounded) → one k20 SCF there + Newton step; rocksalt at a0 on k20 and k24; clamped scan
x = 0.480, 0.485 (k20) with E0 + A d² + B d⁴; soft mode at a_cub = 5.80/5.75/5.70 Å, x = 0.495 (k20, MV 0.005);
PSlibrary scalar/fully-relativistic pair at both structures (k12). 22 evaluations, 61 of 70 units (87 %), ≈ 75 min
of pw.x time. No literature value is used: the starting stencil centre (4.36 Å, 58.3°, 0.470) is only a starting
guess, and the Newton step from the k20 point lands within 0.009° of the reference.

## Gate and Harbor

- **Gate** (`check_task --cpuset 8,9,24,25`, 2026-10-26T19:05Z): passed, no issues. 9 decision points, 10 ablations all
  failing, 7 robustness variants and 4-significant-digit rounding passing, grading deterministic, 8 hidden values,
  task 9.6 MB (six UPF files). Live oracle: 22 evaluations, 61 of 70 units, 0 failed, 4647 s container time; errors
  identical to the offline replay (identical requests reproduce the construction energies exactly).
- **Harbor oracle** (`runs/build-fresh-gete-ferroelectric-well/oracle-r1-fresh-gete-ferroelectric-well`,
  2026-10-26 22:42 – 10-27 00:30Z): reward 1.0, all 8 numerical checks true, charge 61; lab pinned by the daemon to
  cores 2,18,3,19; 108 min wall. Under that host load job times were 32–681 s (two k20 polar SCFs took 654 and 681 s,
  slightly above the 10-minute guideline; 290–520 s in construction).
- **Disk:** peak ≈ 85 MB (70 MB downloads in `~/.cache/ai-theorist/fresh-gete-ferroelectric-well/`, deleted after the
  gate; `fetch_pseudos.py` refetches them) + 4.5 MB of job folders under `runs/`. No per-trial Harbor images remained.

## Deviations from the candidate design

- **SCF-only lab.** The design offered relax / fixed-volume shape relax / vc-relax. A vc-relax from the experimental
  cell took > 25 min on 2 cores and hit the time limit; relaxations were dropped. Every SCF returns the full stress and
  forces, so the structure is found from gradients (stencil + Newton), which is itself a decision point (dp4).
- **Pressure question:** graded as the harmonic soft-mode instability pressure of rocksalt (where ∂²E/∂x² at x = 0.5
  changes sign) instead of the pressure where the relaxed distortion vanishes, which would need relaxations near a flat
  surface. In a Landau picture with η·P² strain coupling both coincide for a continuous transition; continuity is a
  rubric item.
- **No hidden validation SCFs** (the lab could compute any point directly).
- **Spin–orbit added as a graded fork** using the licensed PSlibrary 1.0.0 fully relativistic pair; the protocol itself
  is scalar-relativistic SSSP.
- **Polarization (Born charges / Berry phase) not offered** (would need ph.x and is not needed for the question).
- **Tolerances:** depth 4 % (design 10 %), structure ±0.005 Å / ±0.1° / ±0.0008 (design ±0.002 in x), critical
  pressure ±0.6 GPa (design ±2 % in volume), set between the robustness spread and the nearest wrong paths.
- **Cutoff fork** dropped (pressures converged within 0.04 GPa from 40 Ry).

## Open risks

- **SOC margin:** the reachable k12 route passes by 0.2 meV; the k16 reference R3m spin–orbit run (10 units) exceeds the
  per-request cap. The k12→k16 change (0.4 meV) shows the quantity is not fully k-converged in the lab's range.
- **Timing under load:** 5-unit requests (k20 polar, k12 spin–orbit) take 5–11 min when the host is busy.
- **Reference runs above the cap** (k24 polar, k16 spin–orbit R3m, 715–820 s on loaded cores) are construction-only.
- **Dev-easy risk:** the per-atom and clamped/relaxed forks are recall-level; the discriminating ones are the slow
  k convergence of rocksalt (dp3), the gap closing under pressure (dp6–dp7), the harmonic window (dp8), the
  anharmonic stencil (dp4) and the mixed-family SOC confounder (dp9).
- **Recall:** experimental α = 58.36°, x = 0.475 and literature PBE values fail the structure windows; the depth,
  clamped depth, soft-mode pressure and SOC change are protocol-specific.
