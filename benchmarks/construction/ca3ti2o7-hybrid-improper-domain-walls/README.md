# Construction of ca3ti2o7-hybrid-improper-domain-walls

Released start with a released evaluator: the lab runs the authors' NEP potentials (Skogvoll et al. 2026, Zenodo 22801877, CC BY 4.0) on structures the agent builds.

## Rebuild

1. `./fetch.sh` downloads `Data.zip` (5.98 MB, md5 4007d33af8f3b2b26055b4fbc711abeb) into `~/.cache/ai-theorist/ca3ti2o7-hybrid-improper-domain-walls/raw`.
2. `./queue.sh jobs1.txt`, `jobs2.txt`, `jobs3.txt` run the construction relaxations (4 workers on CPUS="6 7 22 23", runtime image, lab code mounted read-only) and write `runs/<tag>.json` (energies, settings). About 15 minutes in total.
3. `python3 -B make_task.py` writes `environment/agent/catalog.json`, `environment/lab/private/hidden*.json` and `solution/ablation_values.json`, and prints the calibration.
4. `.venv/bin/python -m benchmarks.tools.check_task ca3ti2o7-hybrid-improper-domain-walls --cpuset 6,7,22,23`.

Files: `builder.py` (bulk-to-supercell mapping with parent-symmetry domain variants), `run_wall.py` (one relaxation), `analyze.py`/`xyzio.py` (octahedral order parameters, extxyz reader), `test_ref.py` (compares `nep_numpy.py` with the reference outputs shipped in NEP_CPU's `test_qnep/` folder, https://github.com/brucefan1983/NEP_CPU commit 015cdcd; only its reference text files were read, its code was not run), `test_train.py` (training-set errors), `test_dw.py` (energies of the authors' archived walls).

## Deviations from the catalog design

- **Evaluator.** No NEP library exists in the runtime image, and compiling NEP_CPU was not an option in this build. `environment/lab/nep_numpy.py` is a numpy implementation of NEP4 / NEP4-charge written for this task. Checks: NEP_CPU qNEP test forces and virials reproduced to 1e-13; descriptors identical; on the released training set qNEP-w-BECs gives 1.8 meV/atom, 0.07 eV/Å and BEC RMSE 0.19 e (the BEC-free qNEP gives 2.1 e, as expected); the authors' 11 qNEP-w-BECs wall structures have residual forces ≤ 0.02 eV/Å. Speed: about 0.5 s per evaluation for 576 atoms, so relaxations take 5–40 s. No extra dependency or image layer.
- **Targets.** The design's five absolute wall energies became four graded quantities plus two profile classes. The parent-(100) walls of the design (the authors' "100-180" files) lie at 45° to P and are partially charged: their energy grows with the cell length (ψ: 156, 158, 162, 171, 181 mJ/m² for 8, 12, 16, 24, 32 cells; the authors' 12-cell value 157.9 is reproduced), so no converged value exists. They are a rubric point. The 110 charged wall has no converged energy either, so its graded quantity is the slope of E_DW with domain width (12.2 mJ m⁻² nm⁻¹). It is consistent with the 100 slope (0.59 vs 1.21 mJ m⁻² Å⁻¹, ratio 0.49 ≈ cos²45°), confirming the electrostatic origin.
- **Wall width** is not graded (walls are one to two slabs wide; the width depends on definition). The profile type is graded instead: ψ walls Ising (min |ψ|/ψ₀ = 0.09), φ walls rotational (min |φ|/φ₀ = 0.71).
- **Mechanical boundary condition** (relaxing the normal cell length) changes energies by < 3% (strain ≤ 0.004), so it is a robustness variant, not a decision point. **Wall position** for in-plane walls is irrelevant (body-centring makes both plane types equivalent); it matters for the layer wall (rock-salt 13.4 vs bilayer 49.5 mJ/m²; one bilayer start relaxed into the rock-salt wall, 13.2).
- **90° walls and the switching path** were not used (budget and compactness).

## Calibration

| field | reference | robustness | nearest ablation | tolerance |
|---|---|---|---|---|
| neutral ψ wall (mJ/m²) | 133.6 (16 cells) | 131.1 (normal relaxed), 134.9 (8 cells), 134.1 (oracle, 12 cells) | 121.1 (NEP), 176.5 (axis swapped) | ±5% |
| neutral φ wall | 195.1 | 194.8 (normal), 196.7 (4 cells) | 183.0 (NEP), 206.5 (axis swapped) | ±5% |
| layer ψ wall | 13.66 | 13.07 (normal), 13.99 (6 c), 13.38 (oracle) | 49.5 (bilayer), −8.8 (NEP), 26.8 (÷A) | ±4 |
| charged ψ slope (mJ m⁻² nm⁻¹) | 12.16 (fit 8–24 cells) | 12.15 (normal 12/24), 11.99 (16–24), 12.22 (8–16) | 17.3 (4–8 cells), 6.08 (per L), −2.25 (NEP), 0 | ±15% |
| profiles | ψ ising, φ rotational | — | φ "ising" from P profile | exact |

The NEP ablation is the closest for the φ wall (2.4 mJ/m² outside the window); for the ψ wall it is 5.8 outside.

## Oracle

12 evaluations, charge 17 of 24 (71%), 70 s in the gate container: bulk single point; neutral ψ (8, 12 cells) and φ (8) walls; layer ψ walls (rock-salt and bilayer placement) and layer φ; charged ψ walls with 8, 12, 16, 24 cells; Born charges of the ψ wall. Its route uses no literature values. Gate (2026-10-26): passed, all 10 ablations fail, 3 robustness variants pass, rounding passes. Harbor oracle: reward 1.0 (campaign `build-ca3ti2o7-hybrid-improper-domain-walls`).

## Open risks

- **Construction burden.** Agents must build Ruddlesden–Popper wall supercells in numpy (no ASE). The oracle shows it is ~150 lines; a weaker agent may fail before reaching the physics forks.
- **Categorical fields** are guessable (50% each) from Landau intuition; the numerical fields are not.
- **Local minima.** Agent-built walls with other interpolations may relax to other minima. All our constructions of the same wall type agreed to < 1 mJ/m² except the layer wall placement (intended fork).
- **Evaluator.** A re-implementation: validated, but not the authors' binary. Targets and lab use the same code, so the task is self-consistent.
- **Paper memorization.** The paper (2026) may list wall energies for its 12-cell cells; the graded neutral 110 values coincide with the archived structures' energies (134.5 / 195.0 mJ/m²), so a model that recalled them would pass those two fields, but not the layer-placement and charged-slope fields.
