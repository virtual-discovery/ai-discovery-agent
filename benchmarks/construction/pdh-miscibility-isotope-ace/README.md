# Construction of pdh-miscibility-isotope-ace

Source: Zenodo record 20591788 (CC BY 4.0), the released Pd–H atomic cluster expansion PdH_PBE_v03a with its training and test pickles (Qamar, Ngoipala, Mrovec, Drautz et al., arXiv:2606.09341).

## Rebuild

Set `CACHE=~/.cache/ai-theorist/pdh-miscibility-isotope-ace`, then:

1. `./fetch.sh` downloads the 102 MB zip (md5 checked), unpacks `potential_files/ace/`, and builds a host venv with python-ace 0.4.0rc1, pandas and ase. These are construction tools only.
2. `$CACHE/venv/bin/python -B convert_potential.py $CACHE/potential_files/ace/PdH_PBE_v03a.yaml $CACHE/work/PdH_ace_ctilde.json` converts the B-basis to the C-tilde basis, using the library's exact transformation.
3. `$CACHE/venv/bin/python -B validate_ace.py $CACHE` checks `ace_numpy.py` against python-ace. Maximum differences:
   - on 60 test structures: 6e-12 eV/atom in energy, 4e-8 eV/Å in force, 1e-9 eV/Å³ in stress;
   - ZBL close-contact cases and a finite-difference force check: 1e-7 eV/Å.
4. `$CACHE/venv/bin/python -B extract_dft.py $CACHE $CACHE/work/dft_subset.json` extracts the DFT validation records.
5. Reference calculations run in the runtime image, with `docker run --cpuset-cpus 6,7,22,23 ... -v $PWD:/c -v $CACHE/work:/w`:
   - `python3 -B /c/reference.py base`
   - `python3 -B /c/reference.py dilute`
   - `python3 -B /c/reference.py ordered`
   - `python3 -B /c/reference.py grids`
   - `python3 -B /c/reference.py random 21 3 4,8,13,19,27,40 2`
   - `python3 -B /c/reference.py random 22 3 54,70,85,97,108 2`
   - `python3 -B /c/reference.py random 23 2 2,4,8,12,16,20,24,28,32 2`

   The outputs go to `ref/` (about 1.5 h of compute on four cores).
6. `python3 -B /c/make_task.py`, in the runtime image with the task folder mounted at `/t`, writes:
   - the lab files (evaluator, potential JSON, DFT subset);
   - the catalog and the public budget;
   - the private targets;
   - the ablation and robustness values.

   It also prints the calibration.
7. `.venv/bin/python -m benchmarks.tools.check_task pdh-miscibility-isotope-ace --cpuset 6,7,22,23`

## Deviations from the design (DETAILS.md §48)

- **Pure-numpy ACE evaluator.** The lab uses a numpy port instead of python-ace or LAMMPS, so it needs no extra dependency. `scan` uses an exact incremental update: only the moved atom's environment is recomputed, at about 10 ms per point.
- **No ACE–H₂ trap.** The ACE H₂ molecule is accurate: its minimum is 3 meV from DFT and its bond is 0.754 Å against 0.750 Å. The "ACE H₂ wrong by orders of magnitude" fork is therefore not real here. Instead, the prompt supplies the experimental gas thermochemistry (Δμ and ZPE of H₂ and D₂ at 300 K and 1 bar). The static ½E(H₂) comes from the potential.
- **Fixed-volume interactions are attractive, not repulsive.** In this ACE the fixed-volume H–H pair energies are weakly attractive (about −3 to −6 meV, 108-Pd cell). The design expected them to be repulsive. The fork still holds: fixed-cell configurational energetics are convex and give no gap (T_c ≈ 6 K, ablation dp1). The attraction that opens the gap appears only with zero-stress relaxation.
- **T_c graded as a mean-field quantity.** T_c is graded as the mean-field (ideal configurational entropy) spinodal top of a model defined in the prompt: random O-site occupation, zero-stress static Δe(x), and Cahn coherency for the coherent value. A Monte Carlo or cluster-expansion T_c is not uniquely recoverable (the Fe lesson), so the "mean-field vs MC" fork was dropped. In its place, dp7 (elastic-only Wagner–Horner mean field, 1836 K) and dp8 (ordered small cells, 634 K) test what sets T_c.
- **Sieverts ratio instead of the D/H plateau ratio.** The D/H plateau-pressure ratio is replaced by the dilute Sieverts ratio K_H/K_D. Under the stated x-independent vibrations, the plateau ratio is just its square, so grading both would grade a trivial transform. The absolute H₂ plateau pressure (log10) is graded instead.
- **Extra graded fields:** the T-site fraction at 300 K and the static T–O difference.
- **Provenance of the truth.** None of the graded numbers is archived by the authors (the paper reports MD of nanoparticles, not T_c or isotope effects). All truth comes from our own calculations with the released potential, through the lab code path (fresh-style truth on a released start).

## Calibration

The full table is in `design.json` under `calibration`; values come from `make_task.py`.

| field | reference | tolerance | robustness range | nearest ablation |
|---|---|---|---|---|
| solution_energy_0K_eV | −0.1391 | ±0.012 | −0.1403 … −0.1375 | fixed-volume 4-Pd cell −0.1248; harmonic −0.1746 |
| OT_difference_static_eV | 0.0262 | ±0.010 | 0.0262 | 4-Pd cell 0.0025 / 0.0801 |
| OT_difference_with_zpe_eV | 0.1252 | ±0.015 | 0.1235 … 0.1269 | 4-Pd cell 0.1014; harmonic 0.1568 |
| lattice_expansion_A3_per_H | 3.037 | ±8 % | 3.033 … 3.044 | T site 4.42 |
| tetrahedral_fraction_H_300K | 0.0131 | factor 2 | 0.0122 … 0.0137 | harmonic 0.0020; static 0.42 |
| sieverts_ratio_H_over_D_300K | 2.27 | ±20 % | 2.21 … 2.27 | ZPE only 3.16; harmonic 3.32 |
| log10_plateau_pressure_H_300K_bar | −3.11 | ±0.4 | −3.17 … −3.09 | ordered cells −3.54; harmonic −4.89 |
| Tc_incoherent_K | 419 | ±15 % | 410 … 426 | ordered 634; elastic only 1836 |
| Tc_coherent_K | 282 | ±20 % | 264 … 295 | <111> modulus 178; no coherency 419 |

## Notes on the physics in this potential

- **Pd host:** a0 = 3.9399 Å, C11/C12/C44 = 204/150/65 GPa.
- **Dilute H:** the O-site static solution energy is −0.092 eV. In the 4-Pd cell it comes out at −0.153 eV at zero pressure and −0.078 eV at fixed volume, because of the elastic image interaction. The expansion is Δv = 3.04 Å³/H.
- **Quantum H motion:**
  - The O site is very soft: the harmonic ħω is 36 meV, but the 3D-grid ZPE is 88 meV for H and 58 meV for D, and the first excitation is 70 meV for H.
  - The T site is stiff and nearly harmonic: its ZPE is 187 meV for H.
- **Coherency:** η = 0.0662, the minimum directional modulus is Y_min = 133 GPa (along <100>), and the coherent shift in the spinodal curvature is 0.111 eV.
- **Grid convergence:** the O-site grid (±0.84 Å, h = 0.06) is converged to 0.1 meV. The T grid is limited by the ZBL inner cutoff of the neighbouring Pd; cropping it by one point changes the ZPE by 1.7 meV, well inside the tolerances.

## Oracle

The oracle spends 23 of the 30 units in 21 evaluations (gate container time 912 s):

- 1 cell relaxation of Pd;
- 1 batch energy for the elastic constants;
- 1 H₂ relaxation;
- 3 dilute relaxations (O and T in the 108-Pd cell, O in the 256-Pd cell);
- 2 scans (O grid ±0.70 Å with h = 0.07; T grid ±0.45 Å with h = 0.045);
- 12 random PdHₓ cell relaxations in the 108-Pd cell;
- 1 PdH relaxation.

Its grids and its random samples differ from the reference set, so it also serves as an independent-route check. Container time for the gate is given in check.json.

## Risks

- **Coarse tolerance on the T-site fraction.** The T-site fraction is exponentially sensitive to the T–O free-energy difference, so it is graded within a factor of 2, which corresponds to about 18 meV.
- **The mean-field T_c definition is stated in the prompt.** An agent that runs Monte Carlo with short-range order instead will get a different T_c and fail; the prompt is explicit about this.
- **The agent must write a 3D Schrödinger solver** (sparse finite differences in scipy); the budget allows about two to four grid scans.
- **Dropped DFT record.** One DFT record (DFT-155) was dropped from the catalog because its fractional coordinates happen to contain a graded number's digits.
