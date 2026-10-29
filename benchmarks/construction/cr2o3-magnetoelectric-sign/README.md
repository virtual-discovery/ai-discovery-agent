# Construction of cr2o3-magnetoelectric-sign

This task is a released start. It is built from the Bousquet et al. (J. Phys.: Condens. Matter 2024) DFT archive, Materials Cloud record xd06z-0pp29 (doi 10.24435/materialscloud:ek-fp, CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. `./fetch.sh` downloads `DFT_calculations.tar.gz` (213.6 MB, md5 dd3a0fa4aa1e7b990d9fad58f5f02ca8), `README.txt` and `files_description.md` into `~/.cache/ai-theorist/cr2o3-magnetoelectric-sign/` and unpacks the tarball there (1.4 GB). The archive's object store (rgw.cscs.ch) was returning 503 for about 30 minutes during the build, so the script retries.
2. `PYTHONPATH=<dir with h5py> python3 -B extract.py` writes `pool.json` (1.1 MB). h5py is needed only for the two phonopy hdf5 files (`pip install --target <dir> h5py`). The script asserts the two archive inconsistencies listed below.
3. `python3 -B make_task.py` writes `catalog.json`, the lab `data.json`, `private/hidden*.json` and `solution/ablation_values.json`, and prints the reference, robustness and ablation values. `model.py` holds the formulas.
4. `.venv/bin/python -m benchmarks.tools.check_task cr2o3-magnetoelectric-sign --cpuset 2,3,18,19` writes `check.json`.

## What the archive contains (findings that became decision points)

- **Force-constant eigenmodes.** `Cr_ph_eigenval_pol.hdf5` stores unitless values (Eu_1..4: 11.19, 18.46, 24.69, 35.40). They equal m_H(2πf)² of the all-H phonopy run in eV/Å², so they are force-constant eigenvalues k in eV/Å². The real-mass Eu phonons (9.49, 13.56, 16.85, 19.39 THz) are strong mixtures of the force-constant modes (the lowest phonon is 0.61 FC1 − 0.55 FC2 + 0.40 FC3).
- **Born charges.** VASP prints dipoles in "electrons Angst" (p_ion = −ZVAL·r), so P = −(p_ion + p_elc). They are linear to 0.05 Å.
  - The `Disp_0_0` reference has an electronic dipole of (−0.00027, 0.00014, −0.00738) e·Å, which the displaced runs lack. One-sided differences against it put offset/u into every atom's z row: the acoustic sum rule is violated by 3.7 e at 0.02 Å, and Z*_zz(Cr1) becomes 2.78 instead of 3.15. The in-plane response is unaffected. Central differences or the sum rule remove the offset.
  - **The `Disp_-0_005` folders are inconsistent.** Their POSCARs are displaced by −0.005 Å, but the OUTCARs (positions and p_ion) belong to +0.05 Å runs for Cr1–O4 and to +0.005 Å runs for O5 and O6. A ±0.005 Å central difference with the POSCAR displacements gives Z*_xx(Cr1) = −11.6.
  - Reference values: Z*(Cr1) has xx = yy = 3.005, zz = 3.146 and an antisymmetric in-plane part ∂P_x/∂u_y = +0.209. The O tensors are low-symmetry.
- **Domains.**
  - The nearest Cr–Cr pair is the face-sharing pair Cr1–Cr2 (2.62 Å along the trigonal axis). The next pair, at 2.85 Å, is lateral; the next pair along the same column is 4.13 Å.
  - **VASP** (MAGMOM +5 on Cr1, converged +2.60 μB on Cr1) is in-pointing.
  - **Elk** has `bfcmt` +1 on Cr1 in elk.in, but its converged moment on Cr1 is −2.60 μB, so it is out-pointing.
  - The two codes therefore give opposite raw slopes for every mode. The paper reports both for the out-pointing domain.
  - QE (inputs only) starts out-pointing.
- **Induced moments.**
  - Elk: 28 amplitudes per mode (±0.003 to ±0.05 Å). 111 runs end at the 500-loop limit with RMS potential changes of about 3e-10, because their energy target of 1e-15 is unreachable.
  - **Elk Eu_3 +0.03 Å is anomalous.** It stopped at 129 loops with "Convergence targets achieved", and its moment is 26% below the linear trend (M/a = −0.222 against −0.300).
  - Two Elk folders (Eu_1 −0.008 and −0.01) hold MPI-crash logs next to complete INFO.OUT files.
  - VASP: ±0.025 and ±0.05 Å per mode. 8 of the 16 OUTCARs stop without "General timing", because the runs used EDIFF = 1e-20 and NELM = 800. Their total magnetizations change by about 1e-5 μB or less over the last iterations.
  - Sphere-projected moments are 8–10% below the net moments: the Elk muffin-tin total and the VASP per-ion sums.
  - VASP orbital moments are printed to 3 decimals only.
- **Reference result (out-pointing domain).**
  - α⊥ = −0.9313 ps/m from Elk; VASP gives −0.9445 after the domain flip.
  - Fractions: Eu_1 −0.0002, Eu_2 −0.017, Eu_3 0.485, Eu_4 0.532.
  - dm/du: Eu_3 −0.2996 and Eu_4 −0.1922 μB/Å.
  - The cell volume is Ω = 93.35 Å³ (rhombohedral cell).
  - For comparison, the paper's numbers are Elk −0.921 and VASP −0.857. The paper's VASP number is close to our species-averaged-Born-charge ablation (−0.818) and diagonal-only variant (−0.844), which suggests that the authors averaged Z* over species, as their appendix says.

## Deviations from the catalog design (candidates.json / DETAILS.md section 40)

- **Cell-rotation decision point dropped.** Only the structural relaxation uses the rotated cell (+150° about z). All Born, phonon and moment runs share one frame. α⊥ and the Cr Born tensor are invariant under rotations about z, so no graded value depends on the rotation. The relaxation cell remains purchasable.
- **"Tensor component in a lab frame" replaced.** It is replaced by Cr1 Born-charge components, which test the index convention and the sign of the antisymmetric part, and by the two moment derivatives.
- **Graded dataset domains added** (VASP in-pointing, Elk out-pointing). They are not recallable: the paper presents both codes for the out-pointing domain.
- **QE cross-check not possible.** No QE outputs are archived; only the inputs are purchasable.
- **Orbital moments not usable** (3-decimal output). Spin versus orbital is a rubric point.
- **Annealed domain graded as a categorical.** The instruction gives the experimental temperature dependence (α∥ reverses sign on warming), so the domain follows from the computed sign and thermodynamics, without recall.
- **Two decision points added from discovered evidence:** dp7 (the Elk Eu_3 +0.03 run) and dp11 (the reference-run dipole offset). The design's "linear regime / polarization quantum" fork is replaced by dp5 and dp11: the Born data are linear and far from a quantum.

## Calibration

| field | reference | robustness (pass) | nearest wrong paths (fail) | tolerance |
|---|---|---|---|---|
| α⊥ out-pointing (ps/m) | −0.9313 | VASP flipped −0.9445; Elk ±0.01–0.015 with Born ±0.05 −0.9309; mean −0.9379 | Elk Eu_3 ±0.03 −0.873 (−6.3%); VASP per-ion −0.858; Z* transposed −0.842; muffin-tin −0.841; species-averaged Z* −0.818; phonon modes −1.246; input-field domain +0.931; codes averaged +0.007; hexagonal volume −0.310 | ±5% |
| α⊥ in-pointing | +0.9313 | as above | as above, with the sign reversed | ±5% |
| fraction Eu_3 / Eu_4 | 0.485 / 0.532 | VASP 0.473 / 0.539 | Eu_3 ±0.03: 0.451 / 0.567; species-averaged 0.442 / 0.596; transposed 0.428 / 0.605; phonon modes 0.764 / 0.253 | ±0.025 |
| dm/du Eu_3 (μB/Å) | −0.2996 | VASP −0.2961, Elk window −0.2998 | muffin-tin −0.2718; VASP per-ion −0.2717; ±0.03 pair −0.2609; domain flip +0.2996 | ±5% |
| dm/du Eu_4 | −0.1922 | VASP −0.1976 (+2.8%) | VASP per-ion −0.1778 (−7.5%); muffin-tin −0.1729 | ±5% |
| Z*(Cr1) xx / xy / zz | 3.005 / 0.209 / 3.146 | one-sided 0.02 with ASR 3.008 / 0.206 / 3.149 | mislabeled ±0.005: −11.6 / −0.67 / −12.3; transposed xy −0.209; species-averaged xy 0.000; one-sided without ASR zz 2.78 | ±0.08 / ±0.05 / ±0.08 |
| dataset domains | VASP in, Elk out | — | input field (Elk in); axial second pair (both flipped); same-domain assumption | exact |
| dominant mode, annealed domain | Eu_4_x, out-pointing | — | phonon modes → Eu_3; ignoring the α∥ sign change → in-pointing | exact |

## Oracle

The oracle makes 39 evaluations for 39 of the 44 units (89%), and uses no literature values:

- the force-constant eigenmodes plus the all-H phonons, to establish the units of k;
- the Born reference, +0.02 Å for all 10 atoms, and −0.02 Å for Cr1 and O1. It checks the POSCAR against the OUTCAR positions, checks in-plane linearity and imposes the acoustic sum rule;
- Elk ±0.01 and ±0.02 Å for each mode;
- VASP ±0.025 Å for each mode.

It identifies the domains from the converged moments and checks that VASP and Elk agree to within 5% after the flip.

- **Gate (2026-10-26, cpuset 2,3,18,19):** passed. All 14 ablations fail, the 4 robustness variants pass, rounding to 4 significant digits passes and grading is deterministic. The oracle took 5 s in the container. The package is 0.9 MB.
- **Harbor oracle:** reward 1.0, runtime 47 s (campaign `build-cr2o3-magnetoelectric-sign`).

## Open risks

- **Memorization.** The paper's headline (out-pointing → negative α⊥, Elk −0.92 ps/m) is public. A recalled −0.92 passes the α fields. The dataset domains (the VASP runs are in-pointing, contrary to the paper's presentation), the Born components, the moment derivatives and the fractions are not in the paper, so the whole check cannot be recalled.
- **Traps that need purchases.** The Disp_-0_005 and Elk +0.03 traps bite only if those items are bought. They are natural choices: the smallest displacement for "linear regime", and the larger amplitudes for signal-to-noise.
- **Tight margin on α.** The VASP per-ion-projected path lies at −7.9% against a ±5% window (defensible spread ≤1.5%). The Elk muffin-tin path fails more clearly, at −9.7%.
- **Budget.** A full central-difference Born set (21 units) plus Elk at 4 points per mode (16) plus the eigenmodes (1) and a phonon run (1) is 39 units, which leaves only 5 for VASP. Agents must use one-sided differences with the sum rule, symmetry, or fewer amplitudes.

## Disk

Peak use was about 1.65 GB: 214 MB of downloads plus 1.4 GB unpacked. After the build the raw downloads and the unpacked tree were deleted. `fetch.sh`, `pool.json` and the checksums in `SOURCE.md` are kept.
