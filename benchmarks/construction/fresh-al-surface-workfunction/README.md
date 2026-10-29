# Construction of fresh-al-surface-workfunction

A fresh-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). The lab runs new Quantum ESPRESSO 6.7 calculations of Al slabs and bulk Al: pw.x on two MPI ranks, with pp.x for the planar-averaged potentials.

## Rebuild

1. `python3 -B fetch_pseudo.py --install` downloads the SSSP 1.3.0 PBE-precision archive, its metadata and the license from Materials Cloud record rcyfm-68h65. It checks the md5s and the UPF SHA256, and copies `Al.pbe-n-kjpaw_psl.1.0.0.UPF` and `SSSP_LICENSE.txt` into `environment/lab/`. The 63 MB archive can be deleted afterwards.
2. `python3 -B reference_batch.py` writes `reference_batch.json`, the 127 construction requests.
3. `CPUSET=8,9,24,25 python3 -B qe_batch.py reference_batch.json` runs each request in the runtime image on the builder's cpuset, through `environment/lab/slab_qe.calculate`, the code path that `task_lab.run` uses. Results go to `~/.cache/ai-theorist/fresh-al-surface-workfunction/runs/`. The full set takes about 6 h on two cores; the relaxations dominate.
4. `python3 -B make_task.py` writes:
   - `environment/lab/private/hidden.json` and `private/hidden_values.json`;
   - `solution/ablation_values.json`;
   - `calibration.json`, and the `calibration` block of `design.json`;
   - `reference_runs.json.gz`, the parsed results without the potential profiles, from which `make_task.py` reruns when the run cache is gone.
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-al-surface-workfunction --cpuset 8,9,24,25` writes `check.json`.

`run_reference.sh` chains steps 1–4.

## Lab

- **Protocol:** PBE; the SSSP Al PAW at 30/240 Ry; a = 4.04 Å, the construction's own E(a) minimum of 4.0388 Å on a 24³ mesh, rounded; Marzari–Vanderbilt smearing; local-TF mixing (β = 0.5); conv_thr 2·10⁻¹⁰·N Ry; Γ-centred meshes.
- **Slabs:** symmetric and centred in the cell. The planar averages are computed in numpy from the pp.x cube files: `plot_num=11` gives the electrostatic potential, `plot_num=1` the total local potential.
- **Request limits:**
  - layers: up to 14 at k12, 12 at k18 and 10 at k24;
  - vacuum 8–20 Å; smearing 0.01, 0.02 or 0.04 Ry;
  - optional `nbnd`, from 1.5N + 1 to 4N + 8;
  - symmetric relaxation of the two outer planes, at k18 only, for 5–10 layers;
  - bulk cells: primitive (n³ mesh) or oriented (the 1×1 surface cell with 3 or 2 planes, n×n×kz).
- **Charges:** a slab costs ceil(N·f/4) units with f = 0.5, 1 or 1.6 for k12, k18 or k24; a relaxation costs 3× its SCF; a bulk cell costs 1. The budget is 58.
- **Measured wall times on two cores, with enough bands:**
  - k18: 80 s (N = 7), 135 s (N = 9), 210 s (N = 12);
  - k24: 310 s (N = 12, (111));
  - default bands at k18, N = 12: 440 s (44 SCF iterations);
  - symmetric relaxations at k18: 450–550 s (N = 7–10), 1100–1200 s for N = 11–12, which is why the lab caps relaxations at 10 layers;
  - 18-layer reference slabs: 440–470 s.
- **Failed jobs:** a non-converged SCF or BFGS fails the job, and the job stays charged. Wavefunctions and cube files are always deleted.

## Calibration (eV, J/m², %; full table in `calibration.json` and in `design.json` → `calibration`)

**Reference values.** The reference averages the bulk-terminated slabs with N = 7–14 and 18 at k18 / MV 0.02 / 15 Å, with nbnd = ⌈2.5N⌉ + 4. σ uses the oriented bulk with the same in-plane mesh: kz 24 for (111), 27 for (100).

| quantity | reference | spread of σ_N, φ_N (N = 7–18) | 18-layer slab |
|---|---|---|---|
| σ(111) | 0.8228 | 0.810–0.838 | 0.8243 |
| σ(100) | 0.9504 | 0.941–0.959 | 0.9539 |
| φ(111) | 4.0608 | 4.030–4.095 | 4.0642 |
| φ(100) | 4.2687 | 4.232–4.302 (low at N = 7–9) | 4.2688 |

- **Δd12(111):** +1.061 %, the mean of N = 9–12, which give 0.99, 0.96, 1.12 and 1.18 %. N = 7 gives 0.53 %.
- **Ordering:** φ(100) − φ(111) = +0.208 eV.
- **18-layer energies:** reported as offsets F18 + 18 × 537.46 eV, which are 0.7047 eV for (111) and 0.9442 eV for (100). Their residuals against 18 ε_b + 2Aσ_ref are +1.3 and +3.6 meV.

**Tolerances and margins:**

| field | tolerance | largest robustness error | nearest failing ablation |
|---|---|---|---|
| σ | ±0.008 | 0.0038 (window N = 8–12 or 9–12) | primitive bulk 40³: +0.0122 (111) / +0.0096 (100); k12 ladder +0.0115 (100); slope fit N = 7–12: +0.0146 / −0.0183 |
| φ | ±0.03 | 0.0115 (k24 ladder, (100)) | thin ladder N = 3–6: −0.050 / +0.074; total local potential −0.24 / −0.26 |
| d12 | ±0.3 | 0.10 (single N = 10) | 7-layer relaxation −0.53 |
| 18-layer energy | ±6 meV | 4.2 meV ((111), window N = 9–12) | slope fit N = 8–12 +8.2 (111) / +6.8 (100); primitive bulk −11.2 / −9.4; single N = 7 +11.7 / −12.0 |
| 18-layer φ | ±0.04 | 0.012 | thin ladder −0.053 / +0.074; total potential −0.25 |

The oracle plan errors are: σ −0.0002 / +0.0019, φ +0.002 / −0.006, 18-layer energies −1.4 / −2.3 meV and d12 −0.09. Every ablation fails at least one check; the per-check results are in `check.json`. Repeated runs of identical inputs gave energies equal to 10⁻⁸ Ry, so run-to-run noise is far below every tolerance.

## Decision points and evidence available to the agent

1. **Band count (dp1).** With the pw.x default nbnd, the highest band falls below E_F for (111) slabs from about N = 9. The result reports `highest_band_minimum_eV` and the band count, and the SCF needs 29–44 iterations instead of 15. E_N comes out too high, and σ_N jumps to 0.86–0.93 J/m², which looks like a Boettger divergence. (100) is affected only mildly up to N = 12. The ablation uses real default-band runs.
2. **Bulk reference (dp2).** A separately converged primitive bulk (40³) lies 1.15 meV/atom below the slab-consistent ε_b. The primitive mesh energies themselves scatter by ±3 meV up to 24³, and 24³ happens to lie close to the k18 oriented value.
3. **Slope fit (dp3).** Quantum-size oscillations make a straight-line fit of E_N over the purchasable window depend on the window (σ(111) 0.837 for N = 7–12, 0.802 for N = 8–12). The oriented bulk value is stable to 0.02 meV/atom in kz.
4. **Vacuum level (dp4).** The total local potential still contains v_xc of the spill-out density and plateaus 0.23–0.26 eV below the electrostatic potential. At 15 and 20 Å vacuum it gives 3.79 and 3.80 eV, against 4.03 eV from the electrostatic potential.
5. **Quantum-size effects (dp5).** σ_N and φ_N oscillate. φ(100) is about 0.05 eV lower for N = 7–9 than for N ≥ 10, and thin slabs (N = 3–6) are off by 0.05–0.07 eV. The ordering is resolved (0.21 eV).
6. **Relaxation (dp6).** Unrelaxed slabs show outward forces on the surface plane (+0.04–0.09 eV/Å). The thin-slab relaxation gives half the expansion.
7. **k-mesh (dp7).** On the 12×12 mesh, σ(111) is +0.026 J/m² even with a matched 12×12 oriented bulk. 18×18 agrees with 24×24 to 0.001 J/m² and 0.003 eV, averaged over the ladder.
8. **Normalisation (dp8).** Dividing by A instead of 2A.

The remaining numerical checks are minor. A vacuum of 11–20 Å changes φ by less than 0.001 eV, and 8 Å by +0.008 eV. MV smearing of 0.01 or 0.04 Ry changes φ by 0.02 eV or less.

## Deviations from the candidate design

- **Al(110) dropped.** The question and the grading need only (111) and (100). It would have added cost without a new decision point.
- **Dipole correction and one-sided relaxation dropped.** pw.x keeps the symmetry operations of the bulk-terminated slab, so relaxing one side only stops with "checkallsym: some of the original symmetry operations not satisfied". A version without symmetry would have multiplied the cost. The dipole question remains a rubric item (one-sided relaxations need a dipole correction).
- **Smearing decision point replaced by the band-count decision point.** Smearing between 0.01 and 0.04 Ry changes σ and φ by less than the quantum-size scatter, so it is not a genuine fork here. The default band count, by contrast, is a large effect that the agent has to discover.
- **Hidden thick slabs.** Instead of withheld thicknesses inside the purchasable range, the task uses 18-layer slabs beyond the lab limit of 12 layers. Interpolation inside the range is insensitive to the bulk reference; extrapolation to 18 layers is sensitive. The energies are graded as offsets from an anchor of −537.46 eV/atom, so that 4-significant-digit rounding still resolves 0.1 meV.
- **Tolerances.** σ is graded at ±0.008 J/m², tighter than the design's ±0.02, because the defensible methods agree within 0.004 and the Boettger and slope-fit paths sit at 0.01–0.02. φ is graded at ±0.03 eV and Δd12 at ±0.3 %. The 18-layer energies are graded at ±6 meV instead of ±3 meV per cell, because the quantum-size residual at N = 18 is 1.3–3.6 meV.

## Measured oracle length

The oracle makes 16 evaluations for 52 of the 58 units:

- 2 oriented bulk cells (kz = 18);
- the N = 7–12 ladders of both faces at k18 with nbnd = ⌈2.5N⌉ + 4;
- symmetric relaxations of the 9- and 10-layer Al(111) slabs.

**Gate (check_task, cpuset 8,9,24,25, 2026-10-20).** Passed. The oracle container ran 2870 s; all 12 ablations fail, the 6 robustness variants and 4-significant-digit rounding pass, and grading is deterministic. The oracle errors equal the calibration row `oracle_plan_N7_12`, because the fresh lab runs reproduced the construction energies exactly. The final `check.json` comes from a full rerun on the finished package, written at 12:07 UTC. That run passed with an oracle container time of 3187 s; the per-job times were 28–236 s for the SCFs, and 514 and 599 s for the relaxations.

**Harbor (campaign build-fresh-al-surface-workfunction, trial oracle-r1-fresh-al-surface-workfunction).** Reward 1.0, with every numerical check passing and a charge of 52. The trial took 62 minutes of wall time, 06:12–07:14 UTC. The CPU-pinning daemon pinned the lab container to cores 2,18,3,19. The job times were:

- SCFs: 43–269 s;
- the two relaxations: 628 and 666 s, against 451 and 528 s on dedicated cores in the gate, because the host was busy.

**Relaxation speed check.** Rerunning the reference relaxations with second-order potential and wavefunction extrapolation gave identical positions and energies and no speed-up, so the lab keeps the gated input. The reference runs in the cache are the original ones.

## Open risks

- **Relaxation times.** The 9- and 10-layer relaxations take 450–530 s on two dedicated cores. In the Harbor run, with a busy host, they took 628 and 666 s, slightly over the 10-minute guideline. SCFs stay at 4.5 minutes or less. The job timeout is 1500 s, and agents poll asynchronously. Shortening the relaxations would mean loosening the BFGS force threshold, which leaves the soft d12 mode unconverged; for N = 11–12, d12 creeps by 0.2 % after the forces fall below 1·10⁻³ Ry/bohr.
- **Recall.** PBE work functions of Al are in the literature (about 4.0–4.3 eV), and so is the Al(111) expansion (about 1 %). The graded σ, φ and 18-layer energies are protocol-specific, and the tolerances are tighter than the literature scatter.
- **Narrow margin.** The primitive-bulk ablation fails σ(100) by only 0.0016 J/m². It also fails σ(111) and both 18-layer energies by wider margins.
- **Mesh-dependent Boettger effect.** Its size depends on which primitive mesh an agent picks. A 24³ primitive reference happens to agree with the k18 slab reference, so an agent using it would pass that check honestly.
