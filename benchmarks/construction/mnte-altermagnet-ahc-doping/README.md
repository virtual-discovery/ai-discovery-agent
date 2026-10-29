# Construction of mnte-altermagnet-ahc-doping

This is a released-start task. It is built from Kluczyk et al., *Coexistence of anomalous Hall effect and weak net magnetization in collinear antiferromagnet MnTe*, Phys. Rev. B 110, 155201 (2024) (arXiv:2310.09134). It uses Zenodo record 18863304 (DOI 10.5281/zenodo.18863304, CC BY 4.0), files `dft_data.zip` and `fig6_script.zip`. Provenance and checksums are in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

1. `python3 -B fetch.py [DEST]` downloads the two zips (4.1 MB) through the record's per-file links and verifies them and their 34 members. The default DEST is `~/.cache/ai-theorist/mnte-altermagnet-ahc-doping`.
2. `bash run_in_image.sh make_model.py` converts `wannier90_hr.dat` to `DEST/model/archived.npz` (exact).
3. `bash run_in_image.sh make_cache.py 48 72 100 150 200 300 400` writes the mesh caches `DEST/cache/mesh_N.npz`. It takes about 2.6 h on cpuset 6,7,22,23: the 300 mesh takes 53 min and the 400 mesh 88 min.
4. `bash run_in_image.sh make_task.py` (about 1 min) stages the archive files, writes the catalog, the lab data (`archived_hr.npz` as integer micro-eV, the seven caches, `mnte_tb.py`), the private targets, `solution/ablation_values.json`, the budget in both `public.json` files and `calibration.json` here.
5. `.venv/bin/python -m benchmarks.tools.check_task mnte-altermagnet-ahc-doping --cpuset 6,7,22,23` writes `check.json`.
6. `python3 -B record_oracle.py` writes `oracle_run.json`, then `python3 -B write_design.py` writes `design.json`.

`mnte_tb.py` is the evaluator: it is copied unchanged into the lab and is also the construction engine. All references use it on the lab's own data.

## What the task is

The lab evaluates the archived 32-band spinor Wannier Hamiltonian with numpy, using the archived WannierTools 2.6.1 protocol. I read that protocol from the v2.6.1 source (`sigma_AHC.f90`):

- the Kubo formula in the atomic gauge, with the Wannier centres at the atomic positions given in `wt.in`;
- Fermi–Dirac occupations with kT = `Eta_Arc` (1 meV);
- σ_ab = +(e²/ħ)∫Σ f Ω_c, a sign opposite to the Wannier90 convention;
- a Γ-centred mesh.

The model can be evaluated for the archived Néel vector (+y) or its exact time reverse (−y).

| Request | Cost | Returns |
|---|---|---|
| `bands` | 1 | band energies at ≤ 1000 k, or energies and ⟨σ⟩ at ≤ 500 k |
| `berry_curvature` | 1 | energies and Ω_x, Ω_y, Ω_z of every band at ≤ 500 k |
| `ahc` | 1 (N = 48, 72, 100), 2 (150), 3 (200), 5 (300), 8 (400) | σ_yz, σ_zx, σ_xy and the electron count on an N × N × N/2 mesh, at ≤ 500 μ, kT 0–25 meV |
| `electron_count` | 1 | electron count on any menu mesh |

The mesh requests are evaluated exactly from precomputed histograms: Ω-weighted and state-count histograms in 0.1 meV bins from 4.80 to 8.80 eV, plus the summed curvature of all states below the window. A numpy Kubo sum at 400 × 400 × 200 (32 M k-points) takes 88 min on 4 cores, so it cannot run within a lab request.

The budget is 32 units. The meshes cost 21 units per Néel direction, and k-list requests can be bought without limit.

The graded fields are:

- σ_xy at hole densities of 2.0 × 10²¹ and 6.0 × 10²¹ cm⁻³;
- the number of sign changes of σ_xy on a 5 meV grid from VBM − 1.00 to VBM − 0.25 eV, and the hole density of the sign change;
- the hole density 0.30 eV below the VBM;
- the forbidden components (must be zero);
- the allowed components for Néel ∥ y and Néel ∥ x;
- the electron-doped σ_xy (must be `null`).

## What the archive shows (and the task is built on)

- **The archived Hamiltonian is not the one the authors used.** `wt.in` reads `wannier90_symmed_hr.dat`, the WannSymm output, which is not archived. The archived `wannier90_hr.dat` breaks the magnetic symmetry:
  - E(k) and E(gk) differ by up to 28 meV for the antiunitary operations and by 9 meV for inversion;
  - `bands_symmed.dat` differs from `bands_ori.dat` by up to 26 meV.

  I symmetrized the model myself over the 8 magnetic operations. The result is exactly symmetric, reproduces the DFT grid eigenvalues to 2 × 10⁻⁵ eV, and leaves Γ–A unchanged, because the archived model is already symmetric there. It does **not** reproduce `bands_symmed.dat`, which changes Γ–A by up to 18 meV. I could not identify what WannSymm 1.2.2-alpha did (I tested spin, conjugation, orbital and degeneracy variants).

  Whatever the cause, the archived and the self-symmetrized models give σ_xy within 1–5 S/cm of each other at 200³ in the graded range. The lab therefore serves only the archived model.
- **The archived Hamiltonian gives nonzero σ_yz and σ_zx.** On the 400 mesh the spikes reach |σ_yz| up to about 3600 S/cm and |σ_zx| up to about 17000 S/cm. Away from the spikes they are of order 1 S/cm. The archived output gives |σ_yz|, |σ_zx| < 10⁻⁷ S/cm. C2z forbids both.
- **`3_wannsymm/symmetries.dat` is garbled.** Its 8 matrices are E, I, a sixfold screw, S6-type, two threefold operations ×T, C2z without c/2 ×T and m_z ×T. Ops 5–8 do not map the atoms onto themselves, and ops 3–4 do not map the moments. They contradict the operation labels in `wannsymm.out`, which are correct (E, I, C2z·t, m_z·t, C2y·t·T, m_y·t·T, C2x·T, m_x·T). Read literally, the C3·T operation forbids σ_xy.
- **The VBM is off the high-symmetry lines.** Band 22 peaks at 6.13177 eV at k_frac ≈ (−0.019, −0.001, ±0.420) and its symmetry images. For comparison:
  - the maximum on the WannSymm/Wannier90 band path is 6.11716 eV;
  - the maximum on the DFT 6 × 6 × 8 grid is 6.10700 eV;
  - `wt.in` has E_FERMI = 6.13013 eV, the authors' zero, 1.6 meV below the VBM;
  - the OUTCAR E-fermi is 6.2098 eV, in the gap.
- **The Wannier conduction bands are unfaithful.** The frozen window ends at 6.7 eV (`dis_froz_max`), below the DFT CBM. On the fitting grid the valence bands agree with DFT to 3 × 10⁻⁵ eV, while bands 23–24 are off by up to 1.56 eV. The DFT CBM is 7.441 eV at K; the Wannier CBM is 7.854 eV, and 9.002 eV at K.
- **Convergence of σ_xy.** Coarse meshes give spikes of hundreds to tens of thousands of S/cm, and many spurious sign changes on the 5 meV grid:

  | Mesh | Sign changes |
  |---|---|
  | 48 | 10 |
  | 72 | 11 |
  | 100 | 8 |
  | 150 | 3 |
  | 200, 300, 400 | 1 |
  | archive (801) | 1 |

  At p = 2 × 10²¹ cm⁻³, σ_xy by mesh is −87.7 (48), 14.1 (72), −171.5 (100), −50.5 (150), −68.4 (200), −51.6 (300) and −50.0 (400) S/cm. The archive gives −48.3 S/cm.
- **Smearing.** On coarse meshes a large smearing is not a cure: on the 100 mesh, σ_xy at E_F − 0.05 eV goes from 1 S/cm at 1 meV to about 200 S/cm at 5–25 meV. On the 400 mesh the near-edge values hardly move (14, 13 and 16 S/cm at 1, 5 and 25 meV). In the graded range, kT = 25 meV moves σ_xy(2e21) by +32 % even on the 400 mesh.
- **Other archive facts.**
  - The DFT run has a weak net moment along z (OSZICAR mag_z = 6 × 10⁻⁴ μB/cell), allowed by the same symmetry as σ_xy.
  - The `fig6_script` maps are Ω_z of the top valence band for Néel 0° and 90°. The 0° map has zero mean, and my evaluator correlates with the 90° map (r = 0.83 at kz = 0.45 b3, the file's 0.450).
  - `WT.out` prints Tmin = 100 K, which this routine does not use.

## Deviations from the design (candidates.json / DETAILS.md section 27)

- **σ_xy at E − E_F = −0.2, −0.5 and −0.8 eV became σ_xy at two hole densities (2.0 and 6.0 × 10²¹ cm⁻³), about 0.39 and 0.76 eV below the VBM.**
  - −0.2 eV sits on a 3 meV spike of the archived curve (−79 S/cm, with sign changes at −0.196 and −0.1935 eV) that no affordable mesh resolves. The near-edge region is not converged anywhere, so it is ungradeable.
  - Densities make the p → μ → σ cascade the task and remove the E_F-convention ambiguity.
- **"Number of sign changes in [−1, 0] eV" (six in the archive, five of them near the edge) became the count on a 5 meV grid in [VBM − 1.00, VBM − 0.25]** (one). The 5 meV grid makes the count independent of smearing ≤ 1 meV and of VBM offsets of a few meV for converged meshes. On a 2 meV grid, the 300 mesh with step occupation gives three.
- **The hole density "at −0.2 eV" became the hole density 0.30 eV below the VBM.** This tests the VBM location (off-path), the E-fermi trap and the spin counting.
- **"Allowed components" is graded for Néel ∥ y and ∥ x.** The forbidden components are graded numerically (≤ 0.5 S/cm), and an electron-doped field (must be `null`) was added.
- **"Behaviour under Néel reversal (odd)" is a rubric item, not a graded field.** The service's reversed model is the exact time reverse, so the check is trivial; it is still available as a request.
- **No adaptive-refinement request.** The menu of cached uniform meshes replaces it, because the numpy cost makes on-the-fly dense meshes impossible within 10 min.

## Decision points (8)

| dp | Fork | Kind |
|---|---|---|
| dp1 | Mesh convergence: 100/150/200 fail (sign changes, σ at 2e21); 300/400 pass | convergence, budget |
| dp2 | VBM off the high-symmetry lines (path +8 %, DFT grid +14 % in p) | discovered evidence |
| dp3 | OUTCAR E-fermi in the gap (p −44 %) | metadata trap |
| dp4 | Spinor counting (factor 2) | convention, cascade into σ |
| dp5 | σ_yz, σ_zx of the unsymmetrized model are artefacts | data-quality judgement |
| dp6 | Garbled `symmetries.dat`; Néel ∥ x forbids all components | symmetry, archive defect |
| dp7 | Conduction bands unfaithful: electron doping not supported | regime judgement |
| dp8 | Room-temperature smearing to hide the noise (+32 % at 2e21) | method choice |

dp2, dp3 and dp4 cascade into the hole densities. dp4 also cascades into σ_xy, because the wrong μ(p) is used.

## Calibration

The full table is in `design.json` under `calibration` (from `calibration.json`).

| Field | Reference | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| σ_xy(2.0e21) (S/cm) | −48.25 (archive) | 400: −50.0; 300: −51.6; 400 step: −48.4; 300 at 0–2 meV: −43.7 to −51.6 | 200: −68.4 (+42 %); kT 25 meV: −63.7 (+32 %); spin ×2: −70.1; 100: −171.5 | ±15 % |
| σ_xy(6.0e21) (S/cm) | 108.46 (archive) | 300/400: 114.1–114.4 (+5.2 to +5.4 %) | spin ×2: −77.0 | ±15 % |
| sign changes | 1 | 200/300/400 at kT 0–25 meV: 1 | 100: 8; 150: 3; 72: 11; 48: 10 | exact |
| sign-change density (cm⁻³) | 4.171e21 | 400: 4.138e21; 300: 4.132e21; 200: 4.114e21 | spin ×2: 8.34e21 | ±8 % |
| p(VBM − 0.30 eV) (cm⁻³) | 1.313e21 | meshes 100–400: 1.312–1.313e21; wt.in E_FERMI as VBM: 1.325e21 (+0.9 %) | band-path VBM: 1.418e21 (+8.0 %); DFT-grid VBM: 1.492e21 (+13.6 %); OUTCAR: 7.36e20; spin ×2 | ±5 % |
| forbidden max (S/cm) | 0, 0 | 0 | raw 400: 3589 / 17347 | ≤ 0.5 |
| allowed | y: [xy]; x: [] | — | `symmetries.dat` literal: [] / []; isotropic: [xy] / [xy] | exact |
| electron-doped | null | — | from model: −12.6 S/cm | null |

The model-versus-archive offset (3–7 % at the graded densities) reflects the different, unarchived symmetrized Hamiltonian; it sits well inside ±15 %.

kT = 8.6 meV (the 100 K printed in `WT.out`) passes (+11 %). It is recorded but not an ablation.

## Oracle

The oracle makes 15 evaluations for 28 of 32 units (88 %):

1. VBM search: one 10 × 10 × 8 scan and three 9³ refinements.
2. The DFT-grid comparison against `wannier90.eig`.
3. The band path and a spin-expectation check.
4. The electron count on the 400 mesh.
5. `ahc` on the 400, 300 and 200 meshes.
6. Both Néel signs on the 48 mesh.
7. Two Berry-curvature probes: around the VBM, and the kz plane of the VBM.

It uses no literature value. Every request returns within about 1 s; the gate container took 18 s.

Results:

- **Gate (2026-10-25, cpuset 6,7,22,23):** passed. All 12 ablations fail, the 4 robustness variants pass, rounding to 4 significant digits passes, and grading is deterministic. The task has 8 decision points and 7 rubric items. The package is 9.6 MB, with 4 hidden values and no leaks.
- **Harbor oracle (`runs/build-mnte-altermagnet-ahc-doping/oracle-r1-mnte-altermagnet-ahc-doping`):** reward 1.0, runtime 56 s, 28 of 32 units.
- **Oracle errors:** σ_xy +3.6 % (2e21) and +5.2 % (6e21); sign-change density −0.8 %; p(−0.30) +0.001 %.

The estimated agent length is 45–75 tool calls.

## Open risks

- **Dev-easy risk.**
  - The symmetry fields are memorable from the MnTe altermagnet literature (σ_xy allowed for Néel ⊥ a, forbidden for Néel ∥ a). They count as two exact checks out of eleven.
  - dp3 and dp7 can be defused by reading the files carefully.
  - The strongest discriminators should be:
    - dp1: non-monotonic convergence (the 150 mesh looks converged at 2e21 while 200 does not);
    - dp2: the VBM off the band path;
    - dp5: the artefact components;
    - the cascade.
- **Lucky passes.**
  - The 150 mesh gives one sign change for some VBM offsets (5 meV grid alignment).
  - The 200 mesh with step occupation passes σ_xy(2e21) (+4 %).
  - An agent can thus pass dp1 with less than the converged mesh, though the 100 mesh and the 200 mesh at 0.5–2 meV fail.
- **Archive-versus-model offset.**
  - Graded σ comes from the authors' unarchived symmetrized model; the lab's converged values differ by 3–7 %. The ±15 % tolerance covers this with 8–12 points of margin, and the nearest wrong path (200 mesh) is at +42 %.
  - Tightest margins:
    - σ_xy(6e21): converged +5.4 % against 15 %;
    - p(−0.30): the band-path VBM fails at +8.0 % against 5 %, and the robust spread is under 1 %;
    - σ_xy(2e21): 300 with step occupation passes at −9.4 %.
- **Ambiguity risk on the electron-doped field.** An agent may argue that a qualitative estimate is possible. The evidence (1.56 eV conduction errors on the fitting grid, 0.41 eV CBM offset) is unambiguous in my view.
- **Spikes.** Raw σ_yz and σ_zx on dense meshes reach thousands of S/cm at isolated μ, from near-degenerate pairs at single k-points; σ_xy has the same spikes. An agent scanning finely may be confused, but the graded σ_xy points and the 5 meV-grid count are robust to them.

## Disk

- **Peak:** about 0.15 GB. That is 11 MB of zips (including the unused `fig.zip`), 35 MB extracted, 6 MB of model, 4 MB of caches and 85 MB of scratch (exploratory meshes, the paper PDF, WannierTools sources), plus the gate's lab-image layer (about 10 MB on the shared base).
- **After the build:** raw downloads, extracted members, caches and scratch were deleted. The refetch scripts, `members.md5` and `calibration.json` are kept.
