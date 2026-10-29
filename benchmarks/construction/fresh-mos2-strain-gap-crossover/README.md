# Construction of fresh-mos2-strain-gap-crossover

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7 `pw.x`/`pp.x`
calculations for monolayer 1H-MoS2 that the agent requests; the hidden truth comes from reference calculations made
through the same `task_lab.run` code path.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --network none --cpuset-cpus 10,11,26,27
--user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none` (done by `qe_batch.py`).

1. `python3 -B fetch_pseudos.py --install`: downloads the four PSlibrary 1.0.0 UPF files (scalar and fully relativistic Mo,
   S) and the PSlibrary GPL v2 license into `~/.cache/ai-theorist/fresh-mos2-strain-gap-crossover/`, checks size, MD5 and
   SHA256 (`SOURCE.md`) and installs them into `environment/lab/`.
2. `python3 -B run_reference.py --stage A` (relaxations at 3.13-3.20 A, 35/45/55/70 Ry; a 15 A cell; a 9x9 spin-orbit run),
   then `.venv/bin/python -B make_task.py --stage-a` (writes `reference_lattice.json`: A0 = 3.1866 A from the oracle's
   three-point energy minimum, the S-height line and the Pulay offset of the stress), then `--stage B` (relaxed and clamped
   biaxial ladders, oracle spin-orbit ladder, uniaxial relaxations; 27 runs) and `--stage C` (spin-orbit runs at the
   relaxed geometries and the uniaxial states; 17 runs). Requests are in `reference_plan.py` (built with the oracle's own
   request builder); results append to `reference_runs.jsonl` (resumable). Total reference compute about 5.5 h on 2 cores.
3. `.venv/bin/python -B make_task.py`: applies the estimators of `solution/solve.py` and writes `private/hidden.json`,
   `private/hidden_values.json`, `solution/ablation_values.json` and `design.json` -> `calibration` (every robustness and
   ablation variation graded by `task_lab.grade`).
4. `.venv/bin/python -B replay.py`: offline replay of the oracle against `reference_runs.jsonl` (seconds, no QE).
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-mos2-strain-gap-crossover --cpuset 10,11,26,27` writes `check.json`.

`exploration` runs (107 pw.x runs incl. reference, early lab versions with 8 empty bands and other costs) live in
`runs/build-fresh-mos2-strain-gap-crossover/explore.jsonl` and `diag/` (construction-only 70/90 Ry and dual-12 stress
check); no graded number comes from them.

## Lab

- One request schema (`public.json` -> `request`), every field required: kind (scf | relax, relax scalar-relativistic
  only), relativistic (scalar | full), lattice parameter, 2x2 in-plane strain tensor, S height, in-plane S shift, cell
  height (15/20/25/30 A), ecutwfc (35/45/55/70 Ry, dual 8), Gamma-centred n x n x 1 scf mesh (6-18), up to 60 band k points.
- `run`: pw.x scf/relax (mpirun -np 2, -npool 2, one thread per rank), pp.x plot_num 11 -> planar average in numpy,
  pw.x `bands` at the requested points; returns raw eigenvalues, the potential profile, stress, forces, positions, mesh
  HOMO/LUMO. Inputs/outputs kept; wavefunctions and cube files deleted.
- **Charge = canonical size:** symmetry-distinct scf k points (lab reproduces pw.x's count from the strained, displaced
  cell; matched pw.x on all 107 runs) x3 for relax, + 0.25 per band k point, x4 spin-orbit, x cutoff factor, /25,
  min 1, max 7 units (~10 min on 2 loaded cores; slowest offered requests measured 6-8 min).

## Physics found in construction

- **a0.** Energy minimum of relaxed-ion cells: 3.1864 (45 Ry, 5 points), 3.1865 (55), 3.1863 (70). The analytic stress
  zero is 3.1815 at 45 and 55 Ry (they agree, so it looks converged) but 3.1852 at 70 Ry: a Pulay stress of ~0.26 N/m
  that persists to 55 Ry (35 Ry: zero at 3.148, -1.2 %). Diagnostic at 90 Ry and dual 12 confirmed it is a wavefunction-
  cutoff effect. Band-energy differences are converged at 45 Ry (< 1 meV vs 55; 35 Ry: 4-6 meV).
- **Scalar-relativistic** at a0: VB(Gamma) - VB(K) = +0.018 eV (tensile end ~ -0.1 %); **with SOC** -0.053 eV: the K valence
  band splits by 0.149 eV, Gamma is unsplit, Q conduction lower branch drops ~35 meV, K conduction ~2 meV.
- **Direct range with SOC (relaxed ions):** [-1.436, +0.432] %; CBM moves K -> Q (t = 0.515 of Gamma-K) under compression,
  VBM K -> Gamma under tension. Clamped S height: VB/CB ladders shift ~0.1 eV per %, lower end beyond -2 %.
- **Deformation potentials (vacuum-referenced):** VB(K) -5.13, CB(K) -16.35, Eg(K) -11.22 eV; raw eigenvalues give
  -13.7 (VB) and -7.7 (CB) because the cell-averaged potential moves by ~-8.6 eV per unit strain vs vacuum.
- **Vacuum:** edges vs vacuum identical to 0.1 meV for c = 15-30 A (not a decision point). **k:** 9x9 = 12x12 = 18x18 to 0.1 meV.
- **Uniaxial zigzag -4 %:** transverse strain +0.927 % (energy minimum in eyy; offset-corrected stress gives +0.934 %, raw
  stress +0.736 %); S shift -0.0049 A along armchair; gap 1.7315 eV, VBM at K, CBM at the Q pair along k_x
  (t = 0.516); the +-60 deg Q valleys are 26 meV higher.

## Calibration (full table: `design.json` -> `calibration`)

| check | tolerance | max robust error | nearest failing wrong path |
|---|---|---|---|
| a0 | +-0.003 A | 0.0012 (70 Ry stress zero) | 45 Ry analytic stress zero -0.0049; 35 Ry -0.039 |
| lower end (as a0(1+e)) | +-0.0035 A | 0.0009 (scalar ladder + SOC shifts of all edges) | K-only SOC / scalar -0.0053; 45 Ry stress a0 -0.0049; 9x9 mesh -0.045; clamped -0.026 |
| upper end | +-0.0035 A | 0.0012 | 45 Ry stress a0 -0.0050; clamped +0.0089; scalar -0.018 |
| DP valence_K | 8 % | 2.6 % (scalar-relativistic) | clamped -39 %; raw eigenvalues +168 % |
| DP conduction_K | 8 % | 0.2 % | clamped -23 %; raw +53 % |
| DP gap_K | 6 % | 0.9 % | clamped -15 % |
| uniaxial eyy | +-0.1 % | 0.007 | uncorrected stress -0.19; clamped ions +0.36; clamped transverse -0.93 |
| uniaxial gap | +-0.015 eV | 0.0003 | one Gamma-K line +0.026; clamped transverse -0.017; clamped ions -0.020 |
| edge / valley k | 0.1 1/A | grid vs parabola 0.02 | Q60 instead of Q along k_x (0.7); 9x9 point t=1/3 |

Tightest margins: K-only SOC correction fails the lower end by 0.0018 A (~0.06 % strain); 45 Ry analytic-stress a0 fails
a0 by 0.0019 A and both ends by ~0.0014 A; clamped-transverse fails the gap by 0.002 eV (also fails eyy by far).

## Oracle

`solution/solve.py`: 10 evaluations, 38 of 48 units (79 %): scalar relaxations at 3.17/3.18/3.19 A (energy minimum, S-height
line, stress offset), spin-orbit 9x9 scf + Gamma-K band points at -2/-1/0/+1 % biaxial with relaxed heights, two uniaxial
relaxations (eyy 0.6/1.0 %), one spin-orbit uniaxial run sampling both families of Gamma-K lines. Sum of reference wall
times 42.7 min. Replay reproduces the references exactly.

## Deviations from the candidate design

- PSlibrary scalar/fully-relativistic pairs instead of SSSP; band energies graded with SOC (the scalar-relativistic
  tensile end sits at ~0 % and is decided by ~10 meV). Geometry stays scalar-relativistic.
- a0 from the energy minimum (Pulay stress discovered); strains graded as lattice parameters with the submitted a0.
- No hidden validation SCFs; uniaxial state is zigzag -4 % (CBM at Q along k_x) rather than armchair +3 %.
- Vacuum size and coarse interpolation dropped as decision points (measured negligible); design dp "strain reference"
  kept (dp5) though it is a convention point.
- Absolute deformation potentials graded for the K valence and conduction states only.

## Open risks

- Recall: the qualitative picture (K->Gamma VBM under tension, K->Q CBM under compression, SOC 150 meV at K) is
  well known; graded numbers are protocol-specific.
- dp5 (experimental lattice) and dp6 (vacuum alignment) are convention points strong models know.
- The Pulay trap (dp4/dp9) is subtle: an agent who trusts converged-looking 45/55 Ry stresses fails a0 by 0.002 A;
  it is discoverable by comparing E(a) and stress, or the 70 Ry stress.
- Spin-orbit runs take 4-7 min; low-symmetry spin-orbit 12x12 runs are refused by the 7-unit cap (agents must use 9x9,
  which is converged).

## Gate and Harbor

- **Gate** (`check_task --cpuset 10,11,26,27`, 2026-10-26): passed, no issues. 10 decision points, 12 ablations all failing,
  14 robustness variations and 4-significant-digit rounding passing, deterministic, 10 hidden values, task 9.8 MB.
  Live oracle: 10 evaluations, 38 of 48 units, 0 failed, 2868 s container time.
- **Harbor oracle** (`runs/build-fresh-mos2-strain-gap-crossover/oracle-r1-fresh-mos2-strain-gap-crossover`, 2026-10-26):
  reward 1.0, 47.7 min wall (08:24-09:12 UTC) with the CPU-pinning daemon. Per-trial images removed.
- **Disk:** peak ~35 MB (9.7 MB pseudopotential downloads, deleted after the gate; `fetch_pseudos.py` refetches them;
  13 MB job folders under `runs/`; 1.5 MB gate workdir).
