# Construction of fresh-ge-internal-strain

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7 `pw.x`
calculations for diamond Ge that the agent requests; the hidden truth comes from reference calculations made
through the same `task_lab.run` code path.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --network none --cpuset-cpus 10,11,26,27 --user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none` (done by `qe_batch.py`).

1. `python3 benchmarks/construction/fresh-ge-internal-strain/fetch_pseudos.py --install`: downloads SSSP 1.3.0 PBE
   precision (Materials Cloud rcyfm-68h65) and the PSlibrary file from the QE pseudopotential library into
   `~/.cache/ai-theorist/fresh-ge-internal-strain/`, checks every MD5/SHA256 (see `SOURCE.md`) and installs
   `ge_pbe_v1.4.uspp.F.UPF`, `Ge.pbe-n-rrkjus_psl.1.0.0.UPF` and `SSSP_LICENSE.txt` into `environment/lab/`.
2. `python3 run_reference.py --stage A` (30 equation-of-state runs), then
   `.venv/bin/python -B make_task.py --stage-a` (writes `reference_lattice.json`: A0 = 5.7585, A6 = 5.6012 Å), then
   `python3 run_reference.py --stage B` (57 runs at A0 and A6). Requests are listed in `reference_plan.py`; results
   are appended to `reference_runs.jsonl` (resumable). Total compute 1.3 h on 2 cores.
3. `.venv/bin/python -B make_task.py`: applies the estimators of `solution/solve.py` to the reference runs and writes
   `private/hidden.json`, `private/hidden_values.json`, `solution/ablation_values.json` and `design.json` → `calibration`.
4. `.venv/bin/python -B replay.py`: offline replay (seconds, no QE). A stand-in `lab_client.py` serves the oracle's
   requests from `reference_runs.jsonl` (matched on the canonical request, charged by `task_lab.validate`), grades
   with `task_lab.grade`, and checks every ablation (must fail), robustness variation and 4-digit rounding (must pass).
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-ge-internal-strain --cpuset 10,11,26,27` writes `check.json`.

## Lab

- One request schema (`public.json` → `request`), every field required: kind (`scf` | `relax`), pseudopotential
  (GBRV 3d-in-valence | PSlibrary frozen-3d), lattice parameter, symmetric strain tensor, relative sublattice
  displacement (Cartesian), ecutwfc (30/35/45/60 Ry, dual 8), k mesh (4–16) and shift, occupations
  (fixed | MV smearing) with degauss, relaxation force threshold.
- `run` writes `pw.in`, runs `mpirun -np 2 pw.x -npool 2` with one thread per rank, keeps `pw.in`/`pw.out`,
  deletes wavefunctions, and returns energy, stress (pw.x sign), pressure, forces, band edges or Fermi energy, all
  band energies at the mesh point nearest Γ, the cell, k-point and symmetry counts, and for relaxations the
  ionic steps and final crystal coordinates.
- **Charge = canonical request size.** `validate` counts the symmetry-distinct k points pw.x will use (the lab
  reproduces `kpoint_grid` + `irreducible_BZ` from the strain/displacement symmetry; it matched pw.x on all 87
  reference runs and 13 exploration cells), × a pseudopotential/cutoff factor, ×3 for relaxations, 100 k points
  per unit at GBRV/45 Ry (≈ 20–40 s per unit on 2 cores), minimum 1, maximum 15 per request.

## Physics found in construction (GBRV, PBE, no SOC)

- **Zero pressure.** a0 = 5.7585 Å from P = 0 at 60 Ry (45 Ry: 5.7539; 35 Ry: 5.7174 — Pulay stress −1.4 GPa;
  frozen-3d PSlibrary: 5.7676). On Γ-containing meshes E(Γ2′) − E(Γ25′) = −0.443 eV: PBE Ge is band-inverted.
  A shifted mesh reports a positive "gap" of 0.15 eV; a Γ-centred mesh with fixed occupations fills two thirds
  of the Γ25′ triplet (HOMO = LUMO).
- **Why 6 GPa (deviation from the design, accepted by the coordinator).** At a0 the internal-strain force Λ does
  not converge on any mesh the lab offers (shifted k6/k8/k10/k12: 0.180/0.225/0.250/0.263 Ry/bohr; Γ-centred
  0.29–0.35 depending on occupations; the Maxwell relation Ω ∂σ_yz/∂u_x = Λ fails by 7–55 %). Under 6 GPa
  (a = 5.6012 Å) the direct gap at Γ is open (+0.31 eV) and Γ-centred meshes converge: Λ = 0.3978/0.3997/0.4004/0.4008
  at k8/k10/k12/k16, Maxwell 0.4007. Shifted meshes still undersample the gap region (0.364 at k8, 0.392 at k12).
- **6 GPa values:** C44⁰ = 95.5 GPa, Λ = 0.3994 Ry/bohr per unit ε4, Φ = 0.2565 Ry/bohr², closure C44 = 64.6 GPa
  (32 % softening), ζ = 0.588 (scan 0.588, BFGS relaxation 0.588), ω = 305.5 cm⁻¹, C11 = 135.6, C12 = 62.7 GPa.
- **Grüneisen regime.** E(Γ2′) − E(Γ25′) falls from 0.59 eV (−1 % in a) through 0.31 (a6) to 0.03 eV (+1 %): the optical
  force constant softens anomalously on the expanded side, so ±1 % gives γ = 1.43 (2.43 on a Γ-centred k8 mesh),
  while ±0.5 % and one-sided compression-side steps give 1.026–1.033.
- **Smearing.** MV 0.01–0.02 Ry softens Φ by 3.5 % (ω −5.4 cm⁻¹) and raises C44⁰ by 2 %; MV 0.002 Ry on Γ k12 agrees with fixed occupations.
- **Anharmonicity.** A one-sided displacement along [111] (cubic xyz term) shifts ω by −7.6 cm⁻¹ at 0.02(1,1,1) Å; along x the force is linear to 0.1 cm⁻¹ up to 0.03 Å.
- **Relaxation threshold.** At ε4 = 0.002 the clamped force (8.0e-4 Ry/bohr) is below the 1e-3 default threshold: the relaxation "converges" after one step with zero displacement.
- **Cutoff.** Forces and shear stress are converged at 45 Ry (60 Ry: identical C44⁰, Λ to 0.1 %, ω to 0.01 cm⁻¹); only the pressure needs 60 Ry.

## Calibration (full table: `design.json` → `calibration`)

| field | reference | tolerance | max robust error | nearest failing wrong path |
|---|---|---|---|---|
| a0_A | 5.7585 | ±0.006 Å | 0.0046 (45 Ry pressures) | frozen 3d core +0.0091; 35 Ry −0.041; experiment −0.10 |
| gamma_gap_at_a0_eV | −0.4425 | ±0.05 | 0.0002 | shifted-mesh band edges +0.150; Γ-mesh edges 0 |
| a_6GPa_A | 5.6012 | ±0.006 Å | 0.0034 (45 Ry) | frozen 3d core +0.0100; 35 Ry −0.031 |
| C11_GPa | 135.6 | 2 % | 0.3 % | one-sided 2 % strain −2.7 % |
| C12_GPa | 62.65 | 3 % | 0.2 % | one-sided 2 % strain −7.7 % |
| C44_clamped_ion_GPa | 95.5 | 1.5 % | 0.4 % (energy route with −PΔV) | energy route without pre-stress work +2.6 %; shifted k8 +2.5 %; MV 0.02 +2.2 % |
| C44_relaxed_ion_GPa | 64.64 | 2 % | 1.25 % (MV 0.002) | MV 0.02 +5.1 %; shifted k8 relaxation +8.0 %; shifted k8 closure +11.9 %; clamped +48 % |
| kleinman_zeta | 0.5884 | ±0.012 | 0.0068 (MV 0.002) | shifted k8 −0.052; clamped/loose relaxation 0 |
| raman_frequency_cm-1 | 305.51 | ±1.5 | 0.26 (Γ k16) | smearing MV 0.01 −5.6; MV 0.02 −5.4; frozen 3d +5.7; one-sided [111] 0.02 −7.6 |
| raman_gruneisen | 1.0316 | ±0.05 | 0.005 | MV 0.02 (±1 %) +0.069; ±1 % steps +0.40; ±1 % on Γ k8 +1.40 |

Informational (neither robustness nor ablation): a shifted 12×12×12 mesh gives C44 relaxed 65.80 (+1.8 %) and
ζ 0.577 (−0.011), inside the tolerances; it is the costliest shifted mesh affordable (10 units per displaced cell)
and is close to converged. The frozen-3d-core file passes C44⁰, C44 and ζ individually (+1.2 %, +0.5 %, −0.004)
but fails a0, a(6 GPa) and ω.

Tightest margins: a0 (45 Ry pressures pass by 0.0014 Å; frozen-3d core fails by 0.003 Å), C11 one-sided strain
(fails by 0.7 %), C44 clamped (frozen 3d passes at 1.2 % < 1.5 %), gruneisen MV 0.02 (fails by 0.019).

## Oracle

`solution/solve.py`: 11 evaluations, 27 of 32 units (84 %): pressures at 5.75/5.77 and 5.59/5.61 Å (Γ k12, 60 Ry),
clamped shear e4 = 0.01, displacement u_x = 0.01 Å, shear + u_x = −0.012 Å scan, ±1 % normal strain, and the optical
displacement at a6 × (1 ± 0.005), all Γ-centred 12×12×12, 45 Ry, fixed occupations. Sum of reference wall times
for these requests: 11.2 min. The offline replay reproduces every reference exactly.

## Exploration record

`exploration_runs.jsonl` keeps the parsed results of the 111 exploratory pw.x runs (`explore_driver.py`, a
scratch driver that writes the same pw.x inputs as `task_lab.run` from an earlier request format with Voigt
strains; run with the same image and cpuset). They document the zero-pressure non-convergence (batches b4, b7
at 5.76 and 5.658 Å) and the 6 GPa convergence ladders (b8, b9 at 5.60 Å) that motivated the redesign; no
graded number comes from them.

## Gate and Harbor

- **Gate** (`check_task --cpuset 10,11,26,27`, 2026-10-20): passed, no issues. 10 decision points, 19 ablations all
  failing, 19 robustness variations and 4-significant-digit rounding passing, grading deterministic, 12 hidden
  values, task 1.6 MB. Live oracle: 11 evaluations, 27 of 32 units, 0 failed, 852 s container time. `check.json` was then regraded with `--reuse-state` after a `design.json`-only edit (recorded wall times), which is why it shows ~1 container second.
- **Harbor oracle** (`runs/build-fresh-ge-internal-strain/oracle-r1-fresh-ge-internal-strain`, 2026-10-20): reward
  1.0, all 10 numerical checks true; lab pinned by the coordinator's daemon to cores 6,22,7,23; 23.6 min wall
  (21.5 min of pw.x time under host load 14–20; the same requests took 11.2 min during construction).
- **Disk:** peak ≈ 75 MB (67 MB of downloads in `~/.cache/ai-theorist/fresh-ge-internal-strain/`, deleted after the
  gate; `fetch_pseudos.py` refetches them) + ≈ 7 MB of job folders under `runs/`. No per-trial Harbor images remained.

## Deviations from the candidate design

- **Shear and Raman analysis at 6 GPa instead of zero pressure** (see Physics). a0 and the signed Γ gap at a0 stay graded.
- **No hidden validation SCFs and no withheld strain/volume:** every graded number is derived from several
  agent-bought calculations (pressure roots, finite differences, closure relations); nothing graded is returned by the lab.
- **Strain is a symmetric tensor** (natural for a deformation) and the lab returns pw.x-convention stresses; the
  Voigt factor and stress sign are weak convention points (dp9 only).
- **Semicore treatment**: offered as a second pseudopotential (frozen 3d core, PSlibrary 1.0.0 from the QE library)
  rather than as an ecut/valence switch; the graded protocol is the SSSP GBRV file.
- **Relaxation is offered** with an agent-chosen force threshold (charged ×3); the displacement-scan and closure
  routes are equivalent and cheaper.
- **Experimental C_ij/ζ/ω** are rubric context only.

## Open risks

- The frozen-3d-core file passes the three shear quantities individually; it fails only through a0, a(6 GPa)
  (by 0.003–0.004 Å beyond the tolerance) and ω. An agent who uses it only for convergence tests is not penalized (intended).
- A shifted 12×12×12 mesh (10 units for a displaced cell) lands inside the C44/ζ tolerances (+1.8 %, −0.011);
  only coarser shifted meshes fail. The Maxwell relation exposes the shifted-mesh deficit even at k12.
- The Grüneisen parameter is graded at a regime boundary: it is local (±0.5 % steps or one-sided compression
  steps agree to 0.005) but ±1 % steps already cross into the gap-closing region. Tolerance ±0.05 sits between.
- Γ-centred k16 displaced/sheared cells take ~3 min (7 units); the slowest offered requests (relax, k12, 60 Ry,
  low symmetry) are capped at 15 units ≈ 6–10 min on 2 loaded cores; the pw.x timeout is 25 min.
- Host load roughly doubled pw.x wall times during the Harbor run; agents polling every 20 s will spend more calls.
