# Construction of kagome135-correlation-origin

Released-start task with a computing lab. Source: Crispino et al., *Tunable electronic correlations in 135-kagome metals*, arXiv:2512.22576; data Zenodo record 17951035 (CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

1. `python3 -B fetch.py`: downloads the seven zips (25.6 MB) into `~/.cache/ai-theorist/kagome135-correlation-origin`, verifies the record md5s, extracts and checks the 48 members against `members.md5`.
2. `bash run_in_image.sh make_model.py`: converts the six `turn_*_hr.dat` into `/data/models/*.npz`.
3. `bash run_in_image.sh make_refs.py /data/refs_X.json MODEL:N:U:JU[:mesh] ...`: lab-path SSMF references. The jobs used are listed under `raw_refs` in `calibration.json` (27 jobs, 0.6–9 min each on one core; about 25 min on 4 cores).
4. `bash run_in_image.sh make_task.py`: stages files, writes catalog, lab data, `kag_ssmf.py` copy, private targets, ablation values and `calibration.json`.
5. `python3 -B write_design.py`, then the gate (`check_task --cpuset 2,3,18,19`).

`kag_ssmf.py` is the evaluator: it is copied unchanged into the lab and is the construction engine. Exploration scripts (`explore*.py`, `test*.py`, `validate.py`, `uc_scan.py`, `sigma_variants.py`, `prof*.py`) are kept for the record.

## What the archive shows

- **Six models, 33 Wannier functions each** (3 × 5 d in site-local frames, 5 × 3 pnictogen p, 3 × M s), 9 × 9 × 6 Wannier mesh, exactly Hermitian.
- **Electron counts.** Counting states below the archived `Fermi_level.dat` gives 28.0 (TiSb), 30.9 (VSb), 33.95 (CrSb), 28.0 (TiBi), 34.0 (CrBi), 33.95 (CrBi_9.2). The DFT example (`scf.out`, CsV3Sb5) has 123 valence electrons; the 46 excluded bands (Cs 5s5p, V 3s3p, Sb 4d, Sb 5s) hold 92, leaving 31. The paper's counts (27/30/33) omit the Cs 6s electron donated into the d–p manifold. **The task defines the pristine count by the archive (28/31/34).**
- **Model identity.** Folder `CrBi` has a = 5.49 Å, c = 9.24 Å (the CsCr3Sb5 lattice, the paper's main CsCr3Bi5); `CrBi_9.2` has a = 5.787 Å, c = 9.206 Å (the CsTi3Bi5 lattice, the paper's SM variant). The name encodes neither lattice uniquely.
- **Orbital order and frame.** The `.wout` projection table gives l = 2 with mr = 5, 3, 1, 2, 4 per site, i.e. dxy, dyz, dz2, dxz, dx2-y2, with local z = c and local x radial (hexagon centre → site). The onsite block mixes only functions 3 and 5 (Ag pair), confirming the order. The most renormalized orbital is **dxz** in this frame. The paper calls it B3g = dyz, so its frame has x and y exchanged relative to the archive's projections. The task grades the archive-frame label and states the frame.
- **Kinetic scale.** σ² = Σ_{(b,R)≠(a,0)} |H_ab(R)|² exactly. dxz has the smallest σ in every compound (TiSb 1.715, VSb 1.546, CrSb 1.431, TiBi 1.415, CrBi 1.331, CrBi_9.2 1.180 eV), matching the paper's hierarchy.

## Evaluator and validation

Z2 SSMF (de' Medici gauge c = 1/√(n(1−n)) − 1, λ0 term derived from the energy functional), density-density Kanamori interaction in the S^z form, FLL double counting in that form (d levels shifted by −[(2M−1)U − 5(M−1)J](n̄ − ½) + (U − 3J)/2). Local inter-orbital terms are renormalized as hoppings. Solver: Anderson warm-up, then Powell-hybrid Newton on (Z, d-level shift); the metallic branch is continued from U = 1 eV in 0.5 eV steps (halved on failure, minimum 0.02 eV). 9 × 9 × 6 mesh with time-reversal reduction, 10 meV Fermi smearing. Checks: Z = 1 exactly at U = 0 for any filling; 9 × 9 × 6 and 12 × 12 × 8 agree to ≤ 0.002 in Z.

Comparison with the paper (U = 6 eV, J = 0.18 U, at the paper's counts):

| Case | Paper | This evaluator |
|---|---|---|
| CsTi3Sb5, N = 27 | 0.83–0.90 | 0.866–0.917 |
| CsV3Sb5, N = 30 | (0.74–0.84 quoted at N = 27) | 0.767–0.855 |
| CsCr3Sb5, N = 27 | Ag, B2g ≈ 0.67 | 0.72–0.84 |
| CsCr3Sb5, N = 33 | 0.07–0.35, B3g smallest | 0.135–0.462, dxz (= paper B3g) smallest |
| CsTi3Bi5, N = 27 | 0.77–0.86 | 0.833–0.893 |
| CsTi3Bi5, N = 33 | B3g 0.22, others 0.39–0.51 | dxz 0.66, others 0.73–0.81 (**not reproduced**) |
| CsCr3Bi5 (Sb lattice), N = 33 | 0.01–0.16 | 0.011–0.196 |

The Cr and CsCr3Bi5 regimes and the orbital hierarchy are reproduced; overdoped CsTi3Bi5 is not. Graded values are those of this evaluator, under the conventions stated in the instruction, and never the paper's numbers.

At the archive counts, CsCr3Bi5 has no metallic solution at U = 6 eV: the branch ends near U ≈ 5.8 eV (Sb lattice) and ≈ 5.0 eV (TiBi lattice). This end point is path dependent by up to about 0.3 eV (N = 33 reached U = 6 by one path but not by the lab path), and probes near it cost up to 9 minutes. U_c is therefore not graded; CsCr3Bi5 is graded at U = 4.5 eV, well inside the metallic branch, and its U = 6 eV state is a rubric item.

## Deviations from the design (candidates.json / DETAILS.md section 31)

- **U and J are stated, not hidden** (U = 6 eV, J = 0.18 U, the paper's protocol). A hidden U is not recoverable without an external calibration.
- **CsCr3Bi5 is not withheld from the service.** Its graded Z̄ (U = 4.5 eV) requires the right folder and the right count. Z at U = 6 is not graded (Mott/beyond the metallic branch), and neither is U_c (see above).
- **Added archive-forensics forks:** electron counts (Cs 6s), folder identity, projection order and frame. The "orbital selectivity" fork became two decision points (projection order; paper frame vs archive frame).
- **Decomposition:** the filling share uses ln Z̄, averaged over both exchange orders (Shapley). Its value is insensitive to the count convention (0.463 nominal vs 0.464), so it does not carry a dp1 ablation.
- **Not used:** ARPES/optical comparison (no data in the archive); the hybridization function (rubric only).

## Decision points (7) and calibration

| dp | Fork | Ablation (value) | Reference | Tolerance |
|---|---|---|---|---|
| dp1 | pristine count: archive 28/31/34 vs nominal 27/30/33 | CrSb Z̄ 0.337, Z_min 0.135; CrBi Z̄ 0.728 / 0.538 | 0.268, 0.070; 0.678 / 0.395 | ±0.025 |
| dp2 | CrBi vs CrBi_9.2 identity | swapped Z̄ 0.395 / 0.678; σ 1.180 | 0.678 / 0.395; 1.331 | ±0.025; ±3% |
| dp3 | projection order vs Wannier90 default | dx2-y2 | dxz | exact |
| dp4 | archive frame vs recalled paper label | dyz | dxz | exact |
| dp5 | counterfactual average vs one path / regression | 0.054, 0.874 / 0.846 | 0.464 | ±0.08 |
| dp6 | σ: full moment vs window / d-block / orbital mean | 1.116 / 0.551 / 1.698 (CrBi) | 1.331 (CrBi), 1.715 (TiSb) | ±3% |
| dp7 | Hund ratio at the right filling | 2.50 (nominal) | 3.14 | ±10% |

Robustness: the 12 × 12 × 8 mesh (all SSMF fields, |Δ| ≤ 0.0006 in Z̄, 0.0006 in share, 0.002 in the ratio) and σ from a coarse-mesh PDOS (exact to 1e-4) both pass. The full table is in `design.json` and `calibration.json`.

## Oracle

The oracle uses 13 evaluations of the 20-unit budget (65%). It derives the counts from the archived Fermi levels, the folder identity from the `.wout` lattices, and the orbital names from the `.wout` projection table. It uses no literature values.

- **Gate (2026-10-26, cpuset 2,3,18,19):** passed. All 14 ablations fail, both robustness variants pass, rounding to 4 significant digits passes, grading is deterministic. The oracle took 1132 s in the container (13 evaluations, longest about 9 min: CsCr3Bi5 at U = 6 eV, which has no metallic solution). The package is 30.4 MB.
- **Harbor oracle (2026-10-26):** reward 1.0 (campaign build-kagome135-correlation-origin). solve.py's report wording for the folder name was edited after the gate run; the numbers and lab calls are unchanged, and the Harbor run used the final file.

## Open risks

- **Evaluator, not the authors' code.** Absolute Z depends on our SSMF conventions (FLL form, local off-diagonal terms renormalized). These are stated in the catalog only in outline. Agents cannot reproduce them, but they do not need to: every graded Z comes from lab evaluations.
- **Dev-easy risk.** dp6 (exact σ) and dp5 (explicit averaging) are fairly mechanical. The discriminators are likely the archive forensics: dp1 (Cs electron, counter to the paper), dp2, dp3/dp4 (frame), and budget use.
- **Lab latency.** Evaluations take 0.5–1.5 min, but a request beyond the metallic branch (for example CsCr3Bi5 at U = 6) takes up to about 9 minutes before reporting no metallic solution.
- **The paper's counts conflict with the archive.** The instruction defines the count by the charge-neutral DFT; an agent that follows the paper is graded wrong on dp1.

## Disk

Peak use was about 0.27 GB: zips 26 MB, extracted members 220 MB, models 11 MB and reference JSON. All were deleted after the build; `fetch.py`, `members.md5` and `calibration.json` (with every raw reference evaluation) are kept. The task package is 30.4 MB.
