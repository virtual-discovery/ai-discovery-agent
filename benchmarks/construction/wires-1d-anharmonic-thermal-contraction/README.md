# Construction of wires-1d-anharmonic-thermal-contraction

This is a released-start task. It is built from C. Cignarella, L. Bastonero, L. Monacelli and N. Marzari, *Extreme anharmonicity and thermal contraction of 1D wires* (Nano Lett. 2025; arXiv:2508.07971), and from Materials Cloud record k45dq-ebf35 (DOI 10.24435/materialscloud:fj-20, CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. `python3 -B fetch.py`
   - Downloads README.txt and the 145.7 MB zip into `~/.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction/` and checks their md5.
   - Unpacks only `SSCHA/` and `QHA/` (about 190 MB).
2. `python3 -B extract.py`
   - Standard library only (sqlite3 reads the AiiDA archive).
   - Writes `pool.json` (6.8 MB) to the cache.
3. `bash run_in_image.sh make_task.py /cache/pool.json`
   - Runs in the runtime image, on cpuset 2,3,18,19.
   - Writes the lab data (6.5 MB), the catalog, the private targets, `ablation_values.json` and `calibration.json`, and prints the references.
4. `.venv/bin/python -m benchmarks.tools.check_task wires-1d-anharmonic-thermal-contraction --cpuset 2,3,18,19`
   - Writes `check.json`.
5. `python3 -B write_design.py`
   - Writes `design.json` (decision points, calibration table, length).

Supporting files:

- `dyn.py`: QE/CellConstructor dynamical-matrix reader, 1D Fourier interpolation.
- `parselog.py`: python-sscha log parser.
- `model.py`: heat capacity, stiffness, stress slope, Hessian softest mode.
- `calib.py`: all reference, robustness and ablation numbers.

## What the archive contains, and what it does not

- **SSCHA:**
  - CuC2 (0–2000 K) and TaSe3 (0–1000 K): auxiliary and Hessian matrices on the 1×1×6 grid, 10k configurations, with 100k runs at high T.
  - AuSe2 (0–800 K): auxiliary matrices only, 1×1×10 double-cell grid, 10k configurations.
  - Logs carry the final free energy, Kong–Liu size, centroid forces, the SSCHA stress with errors and the ab initio average stress.
  - TaSe3 has a 400 K log but no 400 K matrices. The TaSe3 1000 K 100k run has matrices but no log.
- **QHA:** `qha_1D.aiida` holds **one** TaSe3 pw.x + ph.x calculation, at c = 3.358 Å (+4.68 %; AiiDA label "0.15"), not a set of strains. The QHA α_L of the design cannot be reconstructed.
- **sscha_1D.aiida:** 100 DFT training configurations of TaSe3 at 600 K. Not used.

## Archive findings the task is built on

- **Stress sign convention.**
  - The SSCHA stress is pressure-like (positive = compressive, as in pw.x). The SSCHA stress exceeds the ab initio average by a positive kinetic term that grows linearly with T. For TaSe3 this term is +0.0023 GPa at 0 K and +0.021 GPa at 1000 K.
  - σ_zz decreases with T for all wires, which means contraction.
  - TaSe3 (200–400 K): dσ_zz/dT = −3.45e−5 GPa/K for the SSCHA stress and −5.42e−5 GPa/K for the average stress (+57 %).
- **Stiffness.**
  - The 1D-ASR DFPT long-wavelength LA branch gives k = 106.4 eV per cell (104.7–107 for q = 0.005–0.05 × 2π/c). The T = 0 SSCHA Hessian gives 100.1.
  - The plain ph.x matrices have Γ residuals of 3–8 cm⁻¹. Their small-q LA slope diverges: 142 eV at q = 0.01, 997 eV at q = 0.002.
  - The strained calculation gives a finite-strain secant of 84.3 eV.
  - α_L(TaSe3, 300 K) = −5.07e−6 /K.
- **The paper's α_L appears to be twice as large.**
  - The paper's α_min values are −1.26e−5 /K (TaSe3) and −7.07e−5 /K (CuC2). Both are about twice what V(dσ/dT)/k gives with k = ∂²E/∂ε²: TaSe3 −5.1e−6 at 300 K, CuC2 about −3.3e−5 at 100–200 K. This is consistent with the "unnormalised Young modulus" Y[0] (the quadratic polyfit coefficient, k/2) of the archived QHA notebook being used as the modulus.
  - It cannot be verified, because the strain energies are not archived. The factor-2 path is ablation dp7.
- **Sign error in the paper's C_V formula.**
  - The SSCHA entropy is the harmonic entropy of the auxiliary frequencies. The Gaussian thermal family satisfies the momentum stationarity exactly, so dF/dT = −S_aux. Hence C_V = Σ c(ω)(1 − (T/ω)∂ω/∂T).
  - Paper Eq. (2) and `cv_example_script.ipynb` (`x = 1 + T/w*derivative`) use +. The classical quartic oscillator checks the sign: the minus sign gives the exact C = 3k_B/4, the plus sign gives 5k_B/4.
  - At 300 K the graded ratios are CuC2 0.680 (frequencies harden; below the harmonic 0.731) and TaSe3 0.979 (frequencies soften; above the harmonic 0.956). The + sign inverts both deviations: 0.782 and 0.932.
- **The Γ rotational mode in the Hessians.** TaSe3 Hessians have a fourth Γ mode at −1.9 cm⁻¹ (0 K) to −9.6 cm⁻¹ (1000 K). This is the axial rotation, which is locked in the auxiliary matrices but not in the Hessian. The genuine instability is the q = 1/3 flexural pair:
  - 11.1 cm⁻¹ at 600 K;
  - +9.2 (100k) or +0.5 (10k) cm⁻¹ at 650 K;
  - −7.0 (100k) or −28 (10k, a Kong–Liu-limited run, at q = 1/2) cm⁻¹ at 700 K;
  - about −11.4 cm⁻¹ at 800 K.
- **CuC2 Hessians** show split, imaginary flexural values (down to −13 cm⁻¹ at 300 K) at q = 1/6 and 1/3. The paper reports CuC2 stable to 2000 K. These values are a rubric point (do not claim a CuC2 instability); they are not graded.
- **Convergence.** Several runs stopped at the Kong–Liu limit (KL ratio 0.5) without meeting the gradient criteria:
  - CuC2 400, 1200, 1600, 1800 and 2000 K (10k), and 1400, 1600 and 1800 K (100k);
  - TaSe3 700 K (10k);
  - every AuSe2 run at T ≥ 150 K (10 populations).

  CuC2 400 K keeps centroid forces of ±0.008 eV/Å.

## Deviations from the design (candidates.json / DETAILS.md section 4)

- **QHA α_L targets dropped.** The QHA archive holds a single strain, so the QHA-versus-SSCHA comparison is a rubric point (item 4). The design's grading of SSCHA α_L for AuSe2 at 200 K was dropped too: all AuSe2 runs at T ≥ 150 K end at the Kong–Liu limit, and the stress slope is noisy (±35 %).
- **Stiffness defined in the instruction and graded.** The design assumed QHA strain energies. Instead, the instruction defines k as the harmonic, relaxed-ion, infinitesimal-strain stiffness. It is graded (±10 %), and α_L uses it. A different stiffness changes α_L, which is the cascade.
- **C_V at 300 K, not at 1000 K, and per vibrational mode of the simulation supercell.**
  - At 1000 K the derivative is not recoverable: defensible fit windows give 0.91–1.03 for CuC2, because the 1200 K run is unconverged and the C–C stretch is noisy.
  - A thermodynamic-limit (Fourier-interpolated) C_V is not recoverable either: interpolated flexural branches go imaginary near Γ.
  - The instruction therefore fixes the supercell and the per-mode normalization, (3N − n₀). TaSe3 was added as a second wire because its correction has the opposite sign.
- **Instability temperature kept.** It is graded against withheld 650/700 K Hessians, ±75 K around the midpoint of the 10k (650.0 K) and 100k (681.6 K) crossings.
- **Dropped decision points.**
  - 10k vs 100k sampling: not separable within the tolerances.
  - Supercell normalization factor 6/3 for C_V: fixed by the stated convention.
- **New decision points.** dp4 (average vs SSCHA stress), dp6 (finite-strain secant), dp8 (derivative sign) and dp10 (interpolation across the withheld temperatures).

## Calibration

The full table is in `design.json` → `calibration`.

| Field | Reference | Defensible variants | Nearest ablations | Tolerance |
|---|---|---|---|---|
| TaSe3 k (eV/cell) | 106.4 | 104.7 (q = 0.05), 100.1 (T = 0 SSCHA Hessian) | finite strain 84.3 (−21 %); plain ASR 142.2 (+34 %) | ±10 % |
| TaSe3 α_L 300 K (1e−6/K) | −5.068 | −5.03 (quadratic 100–500 K), −5.03 (global cubic, Hessian k), −5.15 | finite strain −6.40 (+26 %); plain ASR −3.79 (−25 %); average stress −7.95 (+57 %); k/2 −10.1; sign +5.07 | ±15 % |
| CuC2 C_V ratio 300 K | 0.6798 | 0.6958 (cubic), 0.6899 (linear 200–500), 0.6765 (all converged, cubic) | harmonic with SSCHA ω 0.7307 (+7.5 %); DFPT harmonic 0.7364; per 3N 0.6421 (−5.6 %); + sign 0.7816; Hessian ω 0.7576 | ±3.5 % |
| TaSe3 C_V ratio 300 K | 0.9794 | 0.9787, 0.9794, 0.9779 | harmonic with SSCHA ω 0.9559 (−2.4 %); DFPT harmonic 0.9542 (−2.6 %); + sign 0.9324; Hessian ω 1.064 | ±1.5 % |
| TaSe3 instability (K) | 665.8 | 698.9 (600/800 K, 10k), 695.7 (600 K/800 K 100k), 713.5 (quadratic fit 500–1000 K) | first released unstable 800 (+134 K); Γ rotation 0; auxiliary: none | ±75 K |

The tightest margins are:

- TaSe3 C_V, harmonic with SSCHA frequencies: −2.4 % against ±1.5 %;
- CuC2 C_V per 3N: −5.6 % against ±3.5 %;
- k, finite-strain secant: −21 % against ±10 %;
- α_L, plain ASR and finite strain: ±25 % against ±15 %;
- instability, quadratic fit (robustness): +48 K inside ±75 K.

## Oracle and length

The oracle makes 28 purchases (28 of 36 units) in 3–6 s of container time:

- TaSe3 scf.out;
- 1D-ASR DFPT;
- TaSe3 logs at 0, 200, 300 and 400 K;
- TaSe3 Hessians at 0, 500, 600 and 800 K;
- auxiliary matrices at 100, 200, 300, 500 and 600 K for TaSe3 and CuC2;
- logs for the convergence check.

The catalog is worth 137 units. An agent needs roughly 45–75 tool calls: parsing QE dynamical matrices, Fourier interpolation and the stress/heat-capacity fits.

- **Gate:** passed on 2026-10-26. All 13 ablations fail, all 8 robustness variants pass, and the rounding check passes.
- **Harbor oracle:** reward 1.0 (`runs/build-wires-1d-anharmonic-thermal-contraction`).

## Open risks

- **Graded values are protocol definitions, not archived outputs.** α_L and C_V are computed from archived ingredients with conventions stated in the instruction (stiffness definition, per-mode normalization, interpolation in ω²).
  - The graded α_L is about half of the paper's, and the C_V deviations have the opposite sign to the paper's. A model that recalls the paper fails, which is intended, but a reviewer could argue that the paper's numbers are "the" answer.
  - The README findings above document why: the sign error is derived and checked analytically; the factor 2 is inferred.
- **The quadratic fit including the unconverged CuC2 400 K run passes** (0.687). The convergence fork is therefore only a rubric item, not an ablation.
- **The instability window is wide (±75 K)** because released data (600/800 K) bias interpolation upwards by about 30 K relative to the withheld bracket.
- **Difficulty.** Conventions (sign, units) may be routine for strong models. The discriminating forks are:
  - dp1: the Γ rotation in the Hessian;
  - dp4: average vs SSCHA stress;
  - dp5/dp6: the stiffness route;
  - dp8: the C_V sign. Models may recall the paper's formula.
- **Response sizes.** Matrices are served as raw text, 25–220 kB per purchase (AuSe2 is largest).
