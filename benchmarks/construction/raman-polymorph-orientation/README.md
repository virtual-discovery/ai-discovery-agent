# Construction of raman-polymorph-orientation

This is a released-start task built on the Computational Raman Database v1.0.1 (Bagheri & Komsa, Sci. Data 10, 80 (2023); Materials Cloud record 0m9r0-g3486, CC BY 4.0). Three archived SnF2 records are used unchanged. The "measurements" are simulated from a documented forward model of a hidden sample built on those records (task `SOURCE.md`).

## Rebuild

1. `python3 -B fetch.py [DEST]` downloads the 152.8 MB zip (md5 e2aaea03…) into `~/.cache/ai-theorist/raman-polymorph-orientation`. It extracts records 3216, 2564 and 2390 from `crd.json`, asserts their sha256 and writes `pool.json` (204 kB).
2. In the runtime image, from this folder:
   ```sh
   docker run --rm --cpuset-cpus 6,7,22,23 --user $(id -u):$(id -g) -e HOME=/tmp -v <repo>:/repo \
       -w /repo/benchmarks/construction/raman-polymorph-orientation ai-theorist-benchmark-runtime:20261007 python3 -B make_task.py
   ```
   This writes `catalog.json`, `private/lab_model.json`, `hidden.json`, `hidden_values.json`, both `public.json` files, `solution/solve.py` (`analysis.py` plus a driver) and `solution/ablation_values.json`. It runs a leak scan, replays the reference purchases with the lab's own `task_lab.synthesize` and seeds, and runs `variants.py` on that replay (about 2 min). The variant run writes `calibration.json` and prints the calibration table.
3. `.venv/bin/python -m benchmarks.tools.check_task raman-polymorph-orientation --cpuset 6,7,22,23` writes `check.json`.
4. `python3 -B write_design.py` writes `design.json`, with the calibration table taken from `calibration.json`.

Files:

- `model.py`: the construction physics, used to build the hidden sample. It parses the records, finds the space-group operations by brute force, expands the Born charges with an assert on the acoustic sum rule and the site symmetry, computes LO mixing and the face-B targets.
- `analysis.py`: the analysis library shared by the oracle and the variants.

## The task

The sample is a β-SnF2 crystal (mp-7456, P2_12_12_1) with 3 vol % of randomly oriented α-SnF2 (mp-2843, C2/c), oriented by Q = Rz(48°)Ry(53°)Rz(62°). The frequency scale is s = 1.0237 relative to PBEsol. The spectrometer axis is offset by +2.120 cm⁻¹ (laser at 532.06 nm, nominal 532.00 nm). The third candidate is γ-SnF2 (mp-7457, P4_12_12), the tetragonal parent of β.

The agent buys 320 cm⁻¹ windows of face-A backscattering spectra with any polarizer/analyzer angles, at 1 unit per 60 s accumulation, with a budget of 20 units. The graded fields are:

- the main and minority phases;
- s (±0.2 %);
- the c axis (≤ 6°);
- for face B (normal at 30° in the lab xy plane), I(t,t)/I(z,z) of the strongest mode (within a factor 1.10) and the highest optical phonon for q along n_B (±0.8 %).

## What the data show (and the decision points rest on)

- **Calibration offset.** Stokes and anti-Stokes lines sit asymmetrically about zero on the axis, by 2.12 cm⁻¹, and so do the notch edges. A scale through the origin gives s +0.75 % (dp2).
- **Line density.** Nearest-line matching prefers the line-dense α record, both raw (mean deviation 4.9 cm⁻¹ against 8.8 for β) and scaled (3.5 against 6.1). The forward-model chi²/point is 141 for α, 57 for γ and 2.6 for β (dp3, dp4). γ's strong lines nearly coincide with β's. The discriminators are:
  - the 99/104 cm⁻¹ doublet, whose two components have different polarization patterns;
  - the lowest modes (45 against 36 cm⁻¹, from the soft mode of the ferroelastic transition);
  - the in-plane anisotropy of the A-mode tensors.
- **Mirror branch.** With parallel/crossed settings along the lab axes only, the mirror orientation diag(−1, 1, −1)·Q fits identically (chi² equal to 1e-9) and puts the c axis 65° off (dp5). Rotated settings remove it. The 180° turn about the face normal can never be resolved from face A. Grading accepts it, and both face-B targets are invariant under it.
- **Oblique LO modes.** For q along the face normal, polar modes mix, producing lines at about 404 and 472 cm⁻¹ on the axis and a near-coincident line next to the main line. These lines are absent from the archived TO list. Fitting with TO modes biases the orientation by 12°, and the face-B phonon taken from the TO list is −8.4 % off (dp6).
- **Born charges.** The archived Born charges have 3 rows for 12 atoms (phonopy BORN form, symmetry-independent atoms only). Assigning them unrotated violates the acoustic sum rule and shifts the face-B LO frequency by +5.7 %. Isotropic charges shift it by −7.9 %, and omitting ε∞ by +56 % (dp7).
- **Minority phase.** Adding an orientation-averaged α component lowers chi² by 39 %, against 0.6 % for γ. The test still finds α with only the low windows. A replay of the same sample without the minority phase returns "none" (α drop 0.01 %) (dp8).

## Deviations from the design (candidates.json / DETAILS.md section 7)

- **One composition (SnF2) instead of two.** SnF2 offers a real near-twin pair: β is the ferroelastic distortion of γ, with a transition at 66 °C. It also has a line-dense room-temperature phase (α) and polar Raman modes with large LO–TO splitting. SnF2 is low-profile, and its Raman spectra are not memorable.
- **No IR-reflectance endpoint.** The oblique LO lines in the Raman spectra already test the Born-charge/ε∞ model. The LO target was redefined as the highest optical phonon for q along an unmeasurable face's normal, rather than "the LO of the strongest IR mode", which is direction dependent in a biaxial crystal.
- **The "unmeasured geometry" target** is I(t,t)/I(z,z) of the strongest mode in face-B backscattering. It replaces I_par/I_cross of the two strongest modes:
  - the cross term of the A mode is small for most orientations, which makes I_par/I_cross ill-conditioned (23 % change for 3°);
  - the second-strongest mode is not unique (1.30 against 1.22 in Σ R²).
- **Frequency scale tolerance tightened** from ±0.4 % to ±0.2 %. The defensible spread is ≤ 0.011 %, and the calibration-offset path is +0.75 %.
- **Axis tolerance tightened** from 10° to 6°. The defensible spread is ≤ 0.3°, and the TO-only path is 12.2°.
- **Forks not graded:**
  - Degenerate partners: β has no degenerate modes. A single-partner treatment of γ leaves γ's chi² unchanged (55).
  - Intensity prefactors: a fit without (n+1)/ω is visibly poor (chi²/point 84) but recovers the graded values to 1.1°.
  - Quality flags: γ and α both carry `negative_freq_path = true`, and α is in fact stable.

  All three stay in the rubric.
- **Laser heating considered and rejected.** It was the coordinator's example confounder, but SnF2 is transparent at 532 nm, so strong heating would not be honest. A calibration offset of the laser line (a common, physically honest confounder) is used instead.
- **"Released start" caveat.** The graded truths follow from our hidden sample parameters through the archived tensors; they are not author-archived numbers.

## Calibration

The full table is in `design.json` under `calibration`; the source data are in `calibration.json`.

| Field | Oracle error | Defensible variants | Nearest wrong paths | Tolerance |
|---|---|---|---|---|
| s | +0.001 % | +0.011 % (centroid calibration), −0.0002 %, −0.0002 % | γ −0.24 %, α matching +0.51 %, no offset +0.75 %, raw −2.3 % | ±0.2 % |
| c axis | 0.14° | 0.28°, 0.12°, 0.16° | γ 7.1°, TO-only 12.2°, mirror 64.7° | ≤ 6° |
| face-B ratio (ln) | 0.001 | ≤ 0.0011 | γ 0.128, TO-only 0.216, mirror −0.432 | ≤ ln 1.10 = 0.095 |
| face-B phonon | +0.017 % | ≤ 0.011 % | no offset +0.94 %, TO-only fit +1.30 %, γ +1.91 %, raw −2.3 %, mirror +4.85 %, unrotated Born +5.7 %, isotropic Born −7.9 %, TO list −8.4 %, no ε∞ +56 % | ±0.8 % |
| main phase | β | — | α (line matching), γ (near twin) | exact |
| minority | α (chi² drop 0.389, γ 0.006) | — | none | exact |

The defensible variants are:

- orientation from areas of isolated non-polar lines, with s from the Stokes/anti-Stokes centroid;
- the low windows only;
- a constant linewidth.

## Oracle

The oracle makes 15 purchases for 17 of the 20 units (85 %):

- windows at 200, −200 and 470 cm⁻¹;
- parallel settings at 0, 45, 90, 135 and 30°;
- crossed settings at 0/90, 45/135 and 60/150;
- one 3-accumulation window at 360 cm⁻¹ without analyzer.

The analysis proceeds as follows:

1. Calibrate from the main line and its anti-Stokes partner.
2. For each candidate: run a coarse orientation grid on isolated non-polar line areas (6000 rotations), then a global spectral fit (rotation, ln C, s, offset, width law, and per-spectrum amplitude with a 1.5 % drift prior and a linear background), from 4 starts. LO mixing is included.
3. Pick the phase with the lowest chi²/point.
4. Run the minority test.
5. Compute the face-B predictions.

The analysis takes about 30 s, and about 50 s in the gate container including purchases. The estimated agent length is 45–75 tool calls.

Results:

- **Gate (2026-10-20, cpuset 6,7,22,23):** passed. All 11 ablations fail, the 3 robustness variants pass, rounding to 4 significant digits passes, and grading is deterministic. There are 8 decision points and 8 rubric items, and the package is 0.29 MB with 8 hidden values.
- **Harbor oracle (`runs/build-raman-polymorph-orientation/oracle-r1-raman-polymorph-orientation`):** reward 1.0, runtime 1 min 22 s. The graded errors were s +0.001 %, c axis 0.14°, ratio +0.1 % and LO +0.017 %.

## Open risks

- **Synthetic sample.** The physics is honest within stated simplifications:
  - thin-plate polarization (no birefringence);
  - LO Raman tensors without the electro-optic term;
  - one global frequency scale;
  - no soft-mode temperature dependence near β's 339 K transition.

  An expert agent may question the last point; it cannot change the graded values, which are defined through s.
- **Dev-easy risk.** An agent that writes a full forward model with LO mixing and checks Stokes/anti-Stokes symmetry will pass. Difficulty rests on:
  - noticing the axis offset;
  - expanding the Born rows;
  - modelling oblique LO modes;
  - buying rotated polarizations;
  - not trusting line matching.
- **Thin margins.**
  - γ passes the axis check narrowly (7.1° against 6°). It fails the main-phase, s (−0.24 % against 0.2 %), ratio (0.128 against 0.095) and LO (+1.9 %) checks.
  - The no-offset path fails the LO check by 0.14 percentage points but fails s by 0.55.
- **Chance.** An agent with only lab-axis settings has a 50 % chance of the right branch; the ablation is the wrong branch. The minority phase has three answers (α, γ, none).
- **Leak scan.** The catalog has about 10⁴ six-digit floats, and the first orientation choice (47°, 53°, 62°) collided with a c-axis component. The final orientation (48°, 53°, 62°) has no collisions, including sign-flipped components.

## Disk

- **Peak:** about 160 MB, of which 153 MB is the zip in `~/.cache/ai-theorist/raman-polymorph-orientation` and the rest a few MB of scratch. The zip and the scratch were deleted after the gate and the Harbor run.
- **Kept:** `fetch.py` with checksums, and `pool.json` (204 kB).
