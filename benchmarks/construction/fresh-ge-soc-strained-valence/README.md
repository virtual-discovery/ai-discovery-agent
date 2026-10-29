# Construction of fresh-ge-soc-strained-valence

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7
calculations (`pw.x` SCF, optional `bands` run at agent-chosen k points, optional `bands.x` double-group labels) for
cubic or tetragonally strained diamond Ge with spin–orbit coupling; the hidden truth comes from reference
calculations made through the same `task_lab.run` code path.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --network none --cpuset-cpus 12,13,28,29 --user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none` (done by `qe_batch.py`).

1. `python3 -B fetch_pseudo.py --install`: downloads the SG15 fully relativistic and scalar-relativistic Ge files
   from quantum-simulation.org into `~/.cache/ai-theorist/fresh-ge-soc-strained-valence/`, checks size/MD5/SHA256
   and the CC BY-SA line, installs them into `environment/lab/`.
2. `python3 -B run_reference.py --stage A` (11 cubic pressure runs; a0 = 5.76851 Å from the spin–orbit 80 Ry Γ-k10
   pressures; `reference_lattice.json` holds A0 = 5.7685 and the c-scan grid), then
   `python3 -B run_reference.py --stage B` (34 runs: labelled Γ levels at A0, c scans at ε∥ = −3 … +0.5 %
   (spin–orbit and scalar, two c each), −1 % fixed-c / cubic / k8–k10 80 Ry checks), then
   `python3 -B oracle_chain.py` (stage C: replays `solution/solve.py` offline and runs each request the oracle asks
   for that is not yet in `reference_runs.jsonl`; 7 runs). Total compute ≈ 3.5 h on 2 cores.
3. `bash runs/build-fresh-ge-soc-strained-valence/mk.sh` equivalent: `python3 -B make_task.py` **inside the runtime
   image** (needs scipy): writes `private/hidden.json`, `private/hidden_values.json`, `solution/ablation_values.json`,
   `calibration.json` and `design.json` → `calibration`.
4. `python3 -B replay.py`: offline replay of the oracle and all ablations/robustness variations (seconds, no QE).
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-ge-soc-strained-valence --cpuset 12,13,28,29`.

`explore_plan.py` / `exploration_runs.jsonl` record the 13 first exploratory runs (convergence at a = 5.77 Å; one
invalid run with the fully relativistic file and spin–orbit off is kept and tagged `_INVALID`). `levels.py` is a
small helper for reading them.

## Lab

- One request schema, every field required (`kind` "pw", `a_par_A`, `c_A`, `spin_orbit`, `ecutwfc_Ry` 40/50/60/80,
  `kmesh` 4–12, `kshift`, `occupations` fixed/smearing, `degauss_Ry`, `band_kpoints_inverse_A` (≤ 30), `symmetry_labels`).
- **Charge ≈ 1 unit per minute**, computed in `validate` from the canonical size: the symmetry-distinct k points
  pw.x uses (reproduced by the lab for cubic/tetragonal cells; matched pw.x on every run) × seconds-per-k factors for
  spin–orbit/scalar and cutoff, plus a band-run overhead and per band k point, ÷ 60 s; maximum 10 units per request.
- Silent traps (nothing in the lab warns): the s-like level is inverted below the p manifold; band energies "nearest
  Γ" on a shifted mesh are off Γ; fixed occupations run without complaint in a semimetal; bands.x names the heavy hole
  G_6+ (LH and SO both G_7+); 40 Ry gives converged levels but an unconverged stress.

## Physics found in construction (PBE + SOC, SG15 FR)

- a0 = 5.7685 Å (60 Ry: 5.7677; scalar-relativistic 5.7676; 40 Ry: 5.738 — stress not converged).
- At a0, Γ: Γ8+ (fourfold), Γ7+ split-off at −0.2848 eV (Δ0), odd-parity Γ7− s-like level at −0.566 eV, i.e. below
  the split-off level; the first level above Γ8+ is the p-like Γ6− at +2.30 eV. E_F cuts Γ8+ (semimetal). Levels are
  converged to ≤ 2 meV for k6–k12, 40–80 Ry, MV 0.005–0.02 or fixed occupations on Γ meshes. Away from Γ the s band
  disperses by 0.12 eV within 0.02 Å⁻¹ and one Γ8 branch curves upward (electron-like light band).
- Relaxed biaxial path: −ε⊥/ε∥ ≈ 0.68–0.76 depending on cutoff/k (σ_zz noise ≈ 0.05 GPa in the semimetal);
  reference c/a(−1 %) = 1.0172 (Γ-k10, 80 Ry).
- At −1 % relaxed: HH on top, E(HH)−E(LH) = 0.0722, E(SO)−E(HH) = −0.345 eV. b (symmetric ±0.5 %) = −2.476 eV; linear
  at −1 %: −2.12; 6×6 Bir–Pikus fits −2.44 … −2.55 and reproduce SO−HH to 2 meV.
- E(s)−E(SO): −0.280 (0), −0.149 (−1 %), +0.009 (−2 %), +0.098 (−2.5 %), +0.191 (−3 %): crossing at ε∥ = −1.946 %.

## Calibration (full table: `design.json` → `calibration`)

| field | reference | tolerance | max robust error | nearest failing wrong path |
|---|---|---|---|---|
| a0_A | 5.76851 | ±0.005 | 0.0009 (60 Ry / scalar) | 40 Ry −0.031; experiment −0.11 |
| delta0_at_a0_eV | 0.28476 | ±0.005 | 0.0001 | shifted-mesh point −0.090; band counting +0.281 |
| s_level_at_a0_eV | −0.5656 | ±0.035 | 0.0013 | spin–orbit off +0.089; first level above Γ8 +2.86 |
| c_over_a_at_minus1pct | 1.01720 | ±0.0008 | 0.00044 (scalar 60 Ry) | c fixed −0.0071; exp. reference +0.033 |
| hh_minus_lh_at_minus1pct_eV | 0.07221 | ±0.005 | 0.0015 | c fixed −0.026; scalar +0.055; Koster labels −0.144 |
| so_minus_hh_at_minus1pct_eV | −0.34498 | ±0.008 | 0.0019 | c fixed +0.027; Koster labels +0.072 |
| b_shear_deformation_potential_eV | −2.4756 | ±0.12 | 0.078 (6×6 from −1 %) | linear at −0.5 % +0.155; linear at −1 % +0.36 |
| s_so_crossing_strain_percent | −1.946 | ±0.1 | 0.004 | linear extrapolation from 0/−1 % −0.19; exp. reference +1.9 |

Tightest margins: b (robust 0.078 vs tolerance 0.12 vs nearest wrong path 0.155), c/a (robust 0.00044 vs 0.0008).

## Oracle

`solution/solve.py`: 9 evaluations, 29 of 36 units (81 %): spin–orbit pressures at 5.75/5.79 Å (Γ-k8, 60 Ry), the
labelled cubic run at its a0, two scalar σ_zz runs at −1 % (D from σ_zz = 0), labelled spin–orbit runs at −1 %,
±0.5 % and −2 %. Sum of reference wall times ≈ 20 min. Replay errors: all ≤ 0.2 of the tolerances
(b +0.014, crossing +0.03 %).

## Deviations from the candidate design

- **No projwfc projections:** the SG15 files carry no atomic pseudo-wavefunctions; bands.x double-group labels
  replace them (parity for the s level; the shared LH/SO label for the heavy hole).
- **Withheld-strain validation SCFs dropped:** replaced by the graded s/split-off crossing strain and SO−HH at −1 %,
  which test the same nonlinear LH–SO physics from agent-bought runs.
- **Scalar runs use the scalar SG15 file** (QE 6.7 does not j-average the FR file with lspinorb off: P = −16.5 GPa).
- **Poisson ratio graded as c/a at −1 %** with a tolerance set by the semimetal's stress noise.
- **Effective masses not graded:** in the inverted protocol the light Γ8 branch is electron-like; rubric point.
- Budget 36 units; oracle uses 9 evaluations (below the 15-evaluation guideline; an agent's convergence checks and
  identification runs add 5–10 more).

## Open risks

- An agent that knows "HH on top under compression in Ge" bypasses the label trap (dp7); an agent that requests labels
  finds the s level easily (dp1). Difficulty rests mainly on dp3 (strain reference), dp5 (b nonlinearity) and dp8.
- c/a tolerance (±0.0008) is about twice the k/cutoff spread; Γ-k6 stress may fall outside.
- Host load roughly doubles pw.x wall time; a tetragonal labelled run is 3.5–7 min.

## Gate and Harbor

- **Gate** (`check_task --cpuset 12,13,28,29`, 2026-10-26): passed, no issues. 9 decision points, 15 ablations all
  failing, 13 robustness variations and 4-significant-digit rounding passing, grading deterministic, 8 hidden values,
  task 0.35 MB. Live oracle: 9 evaluations, 29 of 36 units, 0 failed, 1442 s container time.
- **Harbor oracle** (`runs/build-fresh-ge-soc-strained-valence/oracle-r1-fresh-ge-soc-strained-valence`): reward 1.0,
  all 8 numerical checks true, 33.8 min total runtime.
- **Disk:** peak ≈ 10 MB (0.3 MB of pseudopotential downloads, deleted after the gate; `fetch_pseudo.py` refetches
  them) + ≈ 4 MB of job folders under `runs/`. No per-trial Harbor images remained.
