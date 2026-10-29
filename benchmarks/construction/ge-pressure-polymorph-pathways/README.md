# Construction of ge-pressure-polymorph-pathways

This task is a released start with a released evaluator. Its source is the Fantasia et al. (2024) Ge database and DeePMD potential: Materials Cloud 238ad-1bg44, CC BY 4.0.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007`, with `--cpuset-cpus 2,3`, `--network none` and `--user UID:GID`. The cache is mounted at `/h/.cache/ai-theorist/ge-pressure-polymorph-pathways` with `HOME=/h`.

1. `./fetch.sh` downloads `allfiles.zip` (330 MB, md5 d5faa3ad7e44d8b797ff59d1d7ea66c8) into `~/.cache/ai-theorist/ge-pressure-polymorph-pathways/`.
2. `extract.py` parses all 2867 OUTCARs from the zip in place into `pool.json.gz` in the cache. This takes about 1 minute.
3. `build_records.py` writes `environment/lab/records.json.gz` (922 records: CD 219, bSn 152, ST12 245, BC8 150, R8 156) and `environment/agent/catalog.json`.
4. `dump_weights.py` writes the frozen graph's weights to `environment/lab/potential/`. It uses `pbparse.py`, a minimal protobuf reader, so no TensorFlow is needed.
5. `validate_dp.py` writes `dp_validation.json`, with checks against finite differences and against the DFT data.
6. `reference.py` writes `ref_dp.json`. It runs through the lab code path in about 2 minutes:
   - internal-relaxation corrections for every near-equilibrium frame;
   - surrogate-only relaxed enthalpy curves.
7. `python3 -B make_task.py` writes:
   - `private/hidden.json` and `private/hidden_values.json`;
   - `solution/ablation_values.json`;
   - `calibration.json`.
8. `.venv/bin/python -m benchmarks.tools.check_task ge-pressure-polymorph-pathways --cpuset 2,3,18,19` writes `check.json`.

## What the archive contains (findings behind the decision points)

- **Most records are training configurations, not equilibrium points.**
  - The majority of single points are rattled (forces 1–3 eV/Å, shear stress) or strained cells.
  - Each phase also has one relaxation series: ISIF = 3 at 0 GPa, then ISIF = 6 at PSTRESS = 10, 30, 60 and 100 kB (β-Sn also at 250 and 500 kB).
  - Each run starts from the previous run's final structure.
- **ISIF = 6 relaxes the cell only.**
  - CD and β-Sn have no free internal coordinates, so this does not matter for them.
  - For ST12, BC8 and R8 the internal coordinates stay at their 0 GPa values. The compressed frames, and the scaled single points, carry forces up to 0.5 eV/Å.
  - Relaxing them with the potential lowers the energy by up to ~15 meV/atom (ST12 at 10 GPa) and ~3–6 meV/atom near 6 GPa.
  - The correction moves the β-Sn/ST12 crossing from 6.60 to 7.53 GPa, β-Sn/R8 from 6.28 to 6.60 GPa, and β-Sn/BC8 from 5.72 to 5.88 GPa.
- **Printed stresses are unreliable in later relaxation steps.**
  - They deviate from −dE/dV of the E(V) fit by up to 3.8 GPa (CD), with RMS 0.4–1.1 GPa.
  - First steps and single points agree with the fit to 0.2–0.5 GPa.
- **Not every relaxation reached its target pressure.** Several runs exhausted NSW; for example, the CD run at 10 GPa stopped at 5.6 GPa.
- **Smearing differs between phases, but this is not a fork.** The settings are ISMEAR −5, 0 (0.05 eV) and 1 (0.2 eV), but energy(σ→0) and TOTEN differ by < 0.2 meV/atom. Cells of different size agree to 0.5 meV/atom.
- **The DeePMD potential is good for pressures, not for the metastable ordering.**
  - Energy RMSE is 2–4 meV/atom on the five phases.
  - Surrogate-only transition pressures are 0.15–0.36 GPa high.
  - It gives E(ST12) − E(BC8) = 5.2 meV/atom against 2.9 meV/atom from DFT.

## Deviations from the catalog design

- **Graded fields.**
  - Graded:
    - the CD→β-Sn, β-Sn/ST12 and β-Sn/BC8 crossings;
    - E(ST12) − E(BC8) and E(R8) − E(BC8) at 0 GPa;
    - the 0 GPa order.
  - Recorded only as hidden extras: β-Sn/R8 (6.60 GPa) and the ST12/BC8 crossing (0.91 GPa).
- **No graded barrier (design target (c)).**
  - The `Transitions` folder (536 single points) has no path labels.
  - It mixes tetragonal (a, c) grids, rattled images and several unrelated paths (CD↔β-Sn Bain-like, BC8↔hd, high-pressure β-Sn/Imma).
  - Endpoint images are not equilibrium phases and some use different smearing.
  - No barrier is uniquely recoverable, so kinetics is a rubric point (items 5–6).
- **Tolerances are tighter than the design's ±1–1.5 GPa.**
  - The defensible spread is ≤ 0.07 GPa.
  - Several wrong paths sit 0.25–1 GPa away.
- **Truth for the metastable phases combines two sources:** the archived DFT plus the released potential's internal-relaxation energy, run through the lab's code path. This is required by the convention "fully relaxed", which the archive itself does not satisfy. Without the correction, the archived data alone would give a well-defined but physically wrong answer, and that answer is ablation dp4.

## Calibration

| field | reference | robustness range | nearest ablation | tolerance |
|---|---|---|---|---|
| P CD→β-Sn (GPa) | 8.785 | 8.777–8.814 | 8.539 (PSTRESS as P), 8.926 (surrogate) | ±0.2 |
| P β-Sn/ST12 (GPa) | 7.533 | 7.516–7.586; harmonic-correction variant 7.35 | 7.455/7.597 (pass; those ablations fail elsewhere), 6.604 (uncorrected) | ±0.35 |
| P β-Sn/BC8 (GPa) | 5.880 | 5.870–5.903 | 5.723 (uncorrected, passes), 6.071 (surrogate for metastable phases, passes), 6.245 (surrogate only) | ±0.2 |
| E(ST12) − E(BC8) (meV/atom) | 2.86 | 2.74–2.90; harmonic variant 3.2 | 5.18 (surrogate), −2.10 (rattled), 3.80 (printed stress) | ±1.0 |
| E(R8) − E(BC8) (meV/atom) | 5.68 | 5.64–5.71 | 7.74 (surrogate), 1.68 (printed stress), 4.14 (rattled) | ±1.0 |
| order at 0 GPa | BC8, ST12, R8 | same | ST12, BC8, R8 (rattled); BC8, R8, ST12 (printed stress) | exact |

The two variants not in `ROBUSTNESS` were computed in construction scratch. The harmonic variant estimates the correction as −½ F_DFT·K⁻¹·F_DFT, with K the potential's Hessian. It falls within tolerance for every field.

Tightest margins:
- dp7 (PSTRESS as pressure) fails only on P CD→β-Sn, by 0.046 GPa beyond the tolerance.
- dp4 (uncorrected) fails on ST12 by 0.58 GPa.
- dp5b (surrogate for metastable phases) fails on E(ST12) − E(BC8) by 1.3 meV.

## Oracle

The oracle makes 28 evaluations for a charge of 28 of 34 units (82%), and runs in 46 s in the gate container. It uses no literature values:
- 24 relaxation records, the series of every phase;
- 4 `dp` relax_ions requests covering 75 metastable frames.

Results:
- **Gate (2026-10-26):** passed. All 8 ablations fail, the 6 robustness variants pass, and rounding to 4 significant digits passes.
- **Harbor oracle:** reward 1.0, campaign `build-ge-pressure-polymorph-pathways`.

The estimated tool-call count is 45–75.

## Open risks

- **The correction is a model.** The internal-relaxation correction comes from the potential, not from PBE. Its size is consistent between full relaxation and the harmonic estimate from DFT forces to ~20%, which the ST12 tolerance (±0.35) absorbs. A PBE re-relaxation could shift β-Sn/ST12 by ~0.2 GPa.
- **Dev-easiness.** A model that fits E(V) to the relaxation series and never inspects the forces fails only through dp4 (the ST12 crossing). Whether the residual forces get noticed is the main discriminator.
- **Memorization.** Known literature values (PBE Pt ≈ 8–10 GPa, ST12 or BC8 as the lower phase) are too coarse for ±0.2 GPa or ±1 meV.
- **Catalog design.** The catalog shows `ionic_steps`, so the relaxation series are easy to find, and the rattled-frame trap (dp1) bites only agents that buy single points.
