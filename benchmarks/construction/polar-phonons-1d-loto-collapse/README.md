# Construction of polar-phonons-1d-loto-collapse

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** N. Rivano, N. Marzari, T. Sohier, npj Comput. Mater. 9, 194 (2023).
- **Archive:** Materials Cloud record 08gsg-gpd20 (DOI 10.24435/materialscloud:46-wj, CC BY 4.0).
- **Provenance, transformations and what is hidden:** the task's `SOURCE.md`.

## Rebuild

From this folder:

1. **Fetch.** `python3 -B fetch.py [DEST]` downloads the 30.8 MB zip and `files_description.md` through the archive's per-file links, checks their md5, unpacks the zip (134 MB) and verifies `members.json` (sha256 of 270 members). DEST defaults to `~/.cache/ai-theorist/polar-phonons-1d-loto-collapse`.
2. **Generate.** `bash run_in_image.sh make_task.py` runs in the runtime image on cpuset 2,3,18,19 (about 45 s). It writes:
   - `catalog.json`;
   - the lab data (`data/index.json` and gzip blobs, 4.6 MB);
   - the private targets and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json`.

   `model.py` holds the dyn-file parser, the axial Fourier interpolation (QE q2r/matdyn conventions, Wigner–Seitz weights) and the 1D and 3D long-range terms.
3. **Design.** `python3 -B write_design.py` writes `design.json` from `calibration.json` and `check.json`.
4. **Gate.** From the repository root, `.venv/bin/python -m benchmarks.tools.check_task polar-phonons-1d-loto-collapse --cpuset 2,3,18,19`.

## What the archive shows (and the task is built on)

- **Two boundary conditions on the same supercells.**
  - `1D_code` is the authors' QE with a 1D Coulomb cutoff (flagged `assume_isolated='2D'` in the inputs).
  - `3D_code` is standard QE.
  - At Γ both give the same analytic frequencies (no macroscopic field at q = 0).
  - At small q the 3D-periodic LO stays near a finite q → 0 value, about 100 cm⁻¹ above Γ for the tubes. The 1D LO starts at Γ and rises as q² ln(1/qt).
  - The two agree for q ≳ 2π/d. For BN(5,5): 1470.6 against 1474.0 cm⁻¹ at 1/6 |ΓZ|, and 1493.0 against 1493.2 at 1/3.
- **The authors' 1D interpolation.** `model.py` reproduces the archived 1D `matdyn.modes`. The LO agrees within:
  - BN(5,5): 0.03 cm⁻¹;
  - BN(6,6): 0.15 cm⁻¹;
  - BN chain: 0.50 cm⁻¹;
  - GaAs wire: 0.52 cm⁻¹.

  The scheme behind it:
  - **Kernel:** filled dielectric cylinder, Δ(x) = 2 I₁K₁ / [x (I₀K₁ + ε₁D I₁K₀)], with x = qt.
  - **Prefactor:** W = (1 − Δ)/(π t² ε₁D q²), times 4πe²/L.
  - **Ewald damping:** exp(−q²/4α) with α = 1 (2π/alat)², gmax = 14, G along the axis.
  - **Parameters:** t = d/4 and ε₁D = A(ε_zz − 1)/(π t²).

  The parameters were identified from the debug lines of the archived `matdyn.out` (`prefac` = 1/(π t² ε), `eps`). The paper's small-qt expansion of Δ is the expansion of this kernel.
- **The truth (explicit runs) differs from that interpolation by up to 2.2 cm⁻¹** (BN(5,5) at 0.05 |ΓZ|), which sets the scale of the defensible spread.
- **Archive defects:**
  - The BN(4,4) `1D_code/matdyn.modes` is a plain Fourier interpolation without the 1D term (within 1.2 cm⁻¹), although its `matdyn.in` sets `loto_2d`. It is 33 cm⁻¹ below the tube's explicit 1D runs at 0.05 |ΓZ| (task dp3).
  - BN(5,5) `longwave/dyn10` is a copy of `dyn100` (0.01 |ΓZ|), so there is no 0.10 target for that tube.
  - `GaAs_bulk/ZB_phase/dyn1` is a copy of the wurtzite file.
  - `GaAs_nanowire/3D_code` has no `dyn0`.
  - The BN(6,6) and GaAs-wire titles say "GaAs nanowire wz 84" (24 atoms in fact).
- **Radius.** The effective radius is poorly determined.
  - Within the authors' scheme, every t from 9 to 12.5 bohr passes all BN(5,5) checks, and 11.25 to 16 bohr all BN(6,6) checks. That is 1.4–2.0 times the atomic-shell radius (6.59 and 7.89 bohr); t = d/4 (11.25, 12.5) and the paper's Table I radii (11.49, 12.84) both lie inside.
  - The atomic-shell radius fails (dp4).
  - An unconstrained fit to the BN(4,4) runs has a spurious minimum at t ≈ 30 bohr, where ε₁D < 1. The physical minimum is at 8.3 bohr = 1.57 R.
- **GaAs wire.** The GaAs-like LO (244 cm⁻¹) disperses by only 5 cm⁻¹ up to 0.15 |ΓZ|, and plain interpolation is within 1 cm⁻¹. Its strongest IR mode (mass-scaled mode charge) is an H-wagging mode at 516 cm⁻¹.

## Deviations from the design (candidates.json / DETAILS.md section 11)

- **Targets.**
  - (1) The LO at small q for BN(5,5) under 1D boundary conditions is kept and extended. The task grades all five withheld explicit runs, at 0.01, 0.025, 0.05, 0.075 and 0.125 |ΓZ|, with tolerances of ±2–4.5 cm⁻¹ instead of ±1.5. It adds BN(6,6), with six runs from 0.01 to 0.125. The design's ±1.5 cm⁻¹ is below the authors' own interpolation error (2.2 cm⁻¹).
  - Hiding BN(6,6) as well blocks a family-interpolation shortcut. With (4,4) and (6,6) released, the average of their shift curves gives BN(5,5) within 1.5 cm⁻¹ without any 1D physics.
  - (2) The effective radius / dielectric parameter is not graded. It is not uniquely recoverable: the fit is shallow and model dependent (see above). It is a rubric point.
  - (3) The spurious 3D-boundary-condition q → 0 value is kept, as the LO frequency (truth: the archived 3D `matdyn.modes`) with ±3 cm⁻¹ instead of ±15%. The design's ±15% would pass most wrong screenings.
  - (4) The GaAs-nanowire target is dropped: the polar shift of the GaAs-like LO is under 5 cm⁻¹, and plain interpolation already passes. GaAs (wire and bulk) remains released data.
- **Lab.** It serves the archived dyn files verbatim, the released interpolated dispersions (frequencies and axial weights) and program outputs. There is no lab-side interpolation with a user-supplied long-range term: an agent that can build it can also run it, and serving it would reduce the task to choosing parameters.
- **Decision points.**
  - Kept from the design:
    - dp1: boundary conditions;
    - dp2: asymptotic model class;
    - dp3: plain interpolation, now carried by the archive's own defective BN(4,4) interpolation;
    - dp4: radius and dielectric definition;
    - dp6: budget and which q-points, now the Γ-only model against buying the grids;
    - dp7: units.
  - Added:
    - dp5: full kernel against its asymptotic form;
    - dp8: what sets the 3D artefact;
    - dp9: family transfer from BN(4,4).

## Decision points (9)

| dp | Fork | Kind |
|---|---|---|
| dp1 | 1D_code against 3D_code data; the ph.x Γ frequency is not the 3D q → 0 LO | competing explanations |
| dp2 | q² ln q (1D) against a linear (2D) collapse | model class |
| dp3 | 1D long-range term in the interpolation; archived BN(4,4) "1D" interpolation lacks it | defect found by cross-check |
| dp4 | Cylinder radius (electronic extent, not the atomic shell) and ε₁D from the polarizability per length | physical definition; cascade |
| dp5 | Full Bessel kernel against the small-qt expansion inside the interpolation | regime judgement; cascade |
| dp6 | Γ-only model (constant short-range part) fails beyond ~0.075 \|ΓZ\|; the grids must be bought | discriminating data under budget |
| dp7 | \|ΓZ\| = π/L; e² = 2 in Rydberg units | conventions |
| dp8 | 3D artefact set by the cell ε_zz and the supercell volume | which quantity sets the effect |
| dp9 | Shift curve is not transferable across tubes (with or without mode-charge scaling) | transfer |

## Calibration (errors in cm⁻¹ against the archived explicit runs; full table in `design.json` → `calibration`)

| Field (q/\|ΓZ\|) | Truth | Oracle | Defensible variants (range) | Nearest wrong paths | Tol. |
|---|---|---|---|---|---|
| BN55 0.01 | 1325.92 | −0.38 | −0.50 to −0.31 | family +2.4, e² −3.0, linear +3.6, plain −5.6 | ±2.0 |
| BN55 0.025 | 1345.47 | −0.95 | −1.51 to −0.65 | family −3.8, e² −10.6, plain −20.3 | ±3.0 |
| BN55 0.05 | 1381.59 | −0.96 | −2.23 to −0.29 | asymptotic kernel +4.2 (passes here), shell radius +6.6, family −6.9 | ±4.5 |
| BN55 0.075 | 1411.72 | −0.09 | −2.06 to +0.69 | Γ-only +4.2 (passes here), family −7.8, shell radius +8.8 | ±4.5 |
| BN55 0.125 | 1451.34 | +1.04 | −1.03 to +1.48 | cell ε +3.9 (passes here), shell radius +6.0, Γ-only +14.1 | ±4.5 |
| BN66 0.01 | 1332.83 | +0.04 | −0.02 to +0.11 | family +1.9 (passes here), e² −3.0, linear +3.4, plain −6.0 | ±2.0 |
| BN66 0.025 | 1354.59 | −0.75 | −1.06 to −0.41 | shell radius +3.1, family −6.5 | ±3.0 |
| BN66 0.05 | 1392.69 | −0.48 | −1.23 to +0.24 | shell radius +7.8, asymptotic kernel +9.9, family −11.5 | ±4.5 |
| BN66 0.075 | 1423.07 | +0.78 | −0.24 to +1.59 | cell ε +4.0 (passes here), Γ-only +5.9, shell radius +10.2 | ±4.5 |
| BN66 0.1 | 1445.30 | +1.76 | +0.76 to +2.44 | cell ε +5.9, shell radius +9.7, Γ-only +11.5 | ±4.5 |
| BN66 0.125 | 1461.52 | +1.75 | +1.02 to +2.18 | cell ε +5.3, shell radius +6.8, Γ-only +17.2 | ±4.5 |
| BN55 3D q→0 | 1420.41 | +0.03 | −0.38 (1D_code Γ run) | perpendicular ε +24.7, cylinder ε −37.3, analytic Γ −100.4 | ±3.0 |

- **Defensible variants:** the authors' t = d/4; no Ewald damping; Ewald α = 2; t = 1.5 R. Their maximum error is 2.44 cm⁻¹. An Ewald α = 0.5 split reaches +3.94 at BN66 0.1 and still passes (not listed).
- **Other wrong paths:**
  - the 3D-periodic dispersion: +95 at 0.01;
  - plain interpolation: −40 to −50 at 0.05–0.075;
  - the 2π/L zone convention: +54 at 0.05;
  - the asymptotic kernel: +19 at BN66 0.075.
- **Rounding** to 4 significant digits passes.
- **Tightest margins:**
  - the cell-ε ablation fails only through BN66 at 0.1 and 0.125, by 1.4 and 0.75 cm⁻¹;
  - the family transfer fails at BN55 0.01 by 0.4 cm⁻¹, but also at 0.05–0.125 by 2–8 cm⁻¹;
  - the e² ablation at 0.01 by 1.0, but also at 0.025–0.125 by more than 7.

## Oracle and length

- **Purchases:** the oracle makes 27 purchases (27 of 30 units):
  - BN(4,4) 1D: grid `dyn1`–`dyn7` and the explicit runs at 0.025, 0.05 and 0.10 |ΓZ|;
  - BN(5,5) and BN(6,6) 1D: grids and the explicit 0.15 runs;
  - BN(5,5) 3D: Γ.
- **Method:**
  - It calibrates t/R on BN(4,4) by a grid search over R ≤ t ≤ t(ε₁D = 1), giving 1.567.
  - It transfers t = 1.567 R: BN(5,5) 10.32 bohr, BN(6,6) 12.36 bohr.
  - It predicts with the long-range-corrected interpolation and checks it on the 0.15 runs.
- **Run time:** 8–14 s in the gate container.
- **Catalog:** worth 170 units (129 dfpt, 6 interpolations at 2 units, 29 outputs). The budget is 18% of it.
- **Agent length:** 45–75 tool calls expected.

## Results

- **Gate (2026-10-20, cpuset 2,3,18,19):** passed.
  - 9 decision points; all 15 ablations fail.
  - All 5 robustness variants and 4-significant-digit rounding pass.
  - 8 rubric items, 12 hidden values, no leaks.
  - Package 4.7 MB.
- **Harbor oracle (`runs/build-polar-phonons-1d-loto-collapse/oracle-r1-polar-phonons-1d-loto-collapse`):** reward 1.0, solver charge 27, total runtime 50 s (rerun on the final package).

## Open risks

- **Dev-easy risk.** An agent that knows the paper, buys the hidden grids and implements a long-range-corrected axial interpolation with a sensible radius passes. The discriminating forks are:
  - the Γ-only model, which fails at 0.075–0.125 and is the paper's own Eq. 2;
  - the asymptotic kernel;
  - the shell radius;
  - plain interpolation, which the archive's own BN(4,4) interpolation endorses.
- **Implementation burden.** Agents must parse 200–400 kB QE dyn files and build the interpolation themselves. Errors there fail the task without a physics mistake.
- **Truth is DFPT, not a model.** Every route carries an interpolation error of 1–3 cm⁻¹. The tolerances leave about 1–2 cm⁻¹ beyond the worst listed defensible variant at 0.05–0.125.
- **Thin margins:** the cell-ε path (dp4, second ablation) fails by about 1 cm⁻¹.
- **Memorization.** The q² ln q law and the cylinder model are published (2023). No graded number appears in the paper: Table I lists only Γ frequencies.

## Disk

- **Peak:** about 175 MB of my own files:
  - the cache held the 30.8 MB zip and the 134 MB extraction;
  - the scratch held 1.4 MB (the arXiv PDF and exploration scripts);
  - the gate and Harbor run folders held 7 MB each.
- **Cleanup after the Harbor run:** the raw downloads and the gate's lab image were removed. `fetch.py` with `members.json` refetches everything.
