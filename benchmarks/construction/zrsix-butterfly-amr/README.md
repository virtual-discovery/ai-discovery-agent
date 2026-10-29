# Construction of zrsix-butterfly-amr

Released-start task built from Zhang and Yazyev, *Magnetoresistance in ZrSiX (X = S, Se, Te) nodal-line semimetals*, Phys. Rev. B 112, 165147 (2025) (arXiv:2504.12547), Materials Cloud record tyan3-q3j67 (CC BY 4.0). Provenance is in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` through `run_in_image.sh` (cpuset 12,13,28,29; cache mounted at `/data` = `~/.cache/ai-theorist/zrsix-butterfly-amr`).

1. `python3 -B fetch.py` downloads the zip (8.4 MB) and README, verifies them and the 32 members, skips the POTCARs.
2. `bash run_in_image.sh make_model.py` converts the three `wannier90_hr` files to npz.
3. `NTHREADS=4 bash run_in_image.sh make_spline.py` (about 8 min): exact eigenvalues on 128 × 128 × 64, spline coefficients of the band group.
4. `bash run_in_image.sh make_counts.py` (about 20 min): band-energy histograms on the 40/80/120/160 meshes.
5. `bash run_in_image.sh make_cache.py ZrSiSe`, then `ZrSiTe`, `ZrSiS` (about 2.5 h, 3.5 h, 0.5 h on 4 cores): Chambers tables. `TOP=240 bash run_in_image.sh make_cache.py ZrSiSe bc_+000.0 ...` reproduces the 240-mesh convergence check (not shipped).
6. `bash run_in_image.sh make_data.py` stages the archive and assembles the lab data (11 MB).
7. `bash run_in_image.sh make_task.py` writes the catalog, `private/hidden.json`, `solution/ablation_values.json` and `calibration.json` through the lab's own `task_lab.run`.
8. `.venv/bin/python -m benchmarks.tools.check_task zrsix-butterfly-amr --cpuset 12,13,28,29`.

`zrsix_lab.py` is copied unchanged into the lab. `tb.py`, `grid.py`, `spline.py`, `chambers.py` are the construction engine; `analyze.py` reads raw tables; `test_parabolic.py` and `validate_exact.py` are the solver checks.

## Evaluator and validation

- **Chambers solver (numpy).** σ/τ = (2e²/V N_k) Σ (−∂f/∂ε) ⟨v_i v̄_j⟩. Every start point of the 160 × 160 × 80 mesh whose energy lies in a window around the neutral point seeds an orbit (RK4 in k-space arc length, 0.004 Å⁻¹, projection back onto the energy surface). Its contribution is replaced by the orbit average of v_i v̄_j (Shockley tube integral): closed orbits exactly by Fourier series over one period, others by backward recursion to s_max = 10 × 20 T ps. One trajectory gives all 15 Bτ values. Tables are binned by energy (1 meV, cloud-in-cell) and band, so μ and T (10–40 K) are free at lookup; the 40 and 80 meshes are subsets of the 160 mesh.
- **Band interpolation.** Sorted bands have kinks on the nodal lines, which a spline of the bands cannot follow. The group of bands crossing E_F (12–13 for S/Se, 12–14 for Te) is isolated from the others by ≥ 0.28 eV near E_F, so its characteristic-polynomial coefficients are analytic; they are splined on a 128 × 128 × 64 grid of exact eigenvalues and the roots and implicit derivatives give E and v. Spurious speeds at in-group crossings are capped at 7.5 eV Å (the exact model never exceeds 6.9 eV Å near E_F). Without the cap the 80 mesh gave MR ~ 50× too large (a single spurious orbit).
- **Checks.** Parabolic band: periods and σ_xx, σ_xy(ω_cτ) match Drude to 0.1%. Against exact diagonalization of the archived H(k) (`validate_exact.py`, 24 random ZrSiSe Fermi-surface orbits, θ = 0, 45, 90°): same closure flags, period error median 1e-4 (max 3%), ⟨v v̄⟩ error median < 1e-4.
- **No archived curves.** The record contains only inputs (no WannierTools outputs), so the solver could not be validated against archived MR curves. Qualitatively it reproduces the paper: butterfly with maximum near 45–55° for ZrSiSe, peanut pattern with maximum at c for ZrSiTe, open orbits for in-plane fields from the c-axis tubes, non-saturating MR from compensation.

## What the archive shows

- **ZrSiS `input.dat` is the ZrSiSe file**: Se atoms, the ZrSiSe lattice (3.6273/8.3828 Å), E_FERMI 8.9117. The ZrSiS neutral point is 9.6928 eV, 0.78 eV higher.
- **E_FERMI values are not the neutral points.** Counting 12 of 24 spinless bands (24 valence electrons: Zr d 2×4, Si sp 2×4, X p 2×4): ZrSiSe 8.9222 eV (input 8.9117, n_e/n_h = 0.68 at the input value), ZrSiTe 8.8610 eV (input 8.8965, n_e/n_h = 1.90). `NumOccupied = 13` in the ZrSiTe input.
- **ZrSiTe** has a third band (14) at E_F and an unsymmetrized Hamiltonian (degeneracies 1/2/4); its ρ(θ) vs ρ(−θ) asymmetry is ≤ 1% and band 14 changes MR by ≤ 2%, so both are rubric items only.
- **Hall resistivity.** For the compensated crystals ρ_yx(B ∥ c) is a small difference that changes sign at different Bτ on the 40/80/160/240 meshes; at the archived E_FERMI it changes sign at Bτ ≈ 5 T ps (cf. the paper's 4.5 T ps for ZrSiS). Not graded.

## Deviations from the design (candidates.json / DETAILS.md §30)

- **No archived MR curves** to validate against; truth is the lab's validated solver on the authors' Hamiltonians (released-start: the authors' models, our calculation). Stated in SOURCE.md.
- **All three compounds are served** (ZrSiSe not withheld): with no Hamiltonian access an agent could not compute ZrSiSe at all. Graded quantities are derived (peak angle, ratios, exponent, neutral μ), never served.
- **Compensation ratio n_e/n_h** is 1 by definition at neutrality; replaced by the neutral μ (three compounds) and n (Se, Te).
- **Hall sign-change field dropped** (not converged; see above) and replaced by the field exponent of MR(B ∥ b), which tests compensation directly.
- **"Open orbits at 45° in the ac plane"** is answerable (none at the neutral point: open fraction < 0.1%) but a yes/no field is guessable; it is a rubric item, and dp5 tests the closed-orbit vs Chambers fork through the peak angle.
- **Added**: ZrSiTe anisotropy ratio and a two-sample Kohler comparison (dp6).
- **ZrSiTe asymmetry and band selection** are not decision points (effects ≤ 2%).

## Decision points (6)

| dp | Fork | Kind |
|---|---|---|
| dp1 | E_FERMI of the inputs vs electron counting (cascades into every transport field) | discovered evidence, cascade |
| dp2 | ZrSiS input.dat is a copy of the ZrSiSe file | archive defect |
| dp3 | Electron count of the spinless model (NumOccupied = 13; spin factor 2 in n) | counting convention |
| dp4 | Mesh convergence: 40 fakes saturation and a 36° peak, 80 is 18% off in MR(b)/MR(c) | convergence under budget |
| dp5 | Chambers on the real Fermi surface vs closed-orbit effective-mass model | model class |
| dp6 | Kohler: compare samples at their own Bτ | units / scaling |

## Calibration

| Field | Reference (160) | Defensible (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| μ (eV) | S 9.6928, Se 8.9222, Te 8.8610 | 120 mesh: ≤ 0.7 meV | input E_F: −10.5 / +35.6 meV; ZrSiS copy −781 meV; 13 bands +968 meV | ±3 meV |
| n (cm⁻³) | Se 1.509e20, Te 4.007e20 | 120: +0.2 / −0.3%; 80: +5.0 / +0.6% | one spin −50% | ±6% |
| peak angle (°) | 54.8 | cubic 54.7; 240 mesh 56.1 | 80 mesh 51.6; input E_F 46.8; 40 mesh 35.6; mass model 90 | ±3° |
| MR(b)/MR(c) | 2.607 | 240 mesh 2.845 (+9%) | 80 mesh 2.147 (−18%); input E_F 1.937; 40 mesh 0.931 | ±15% |
| exponent (B ∥ b) | 1.610 | 240 mesh 1.529; 80 mesh 1.415 | input E_F 0.264; 40 mesh 0.267 | ±0.3 |
| ZrSiTe MR(c)/MR(b) | 12.83 | 240 mesh 11.67 (−9%); 120 μ 12.40 | none fails this field alone (40 mesh 9.72, −24%; 80 mesh 10.13); those paths fail other fields | ±25% |
| Se/Te sample ratio | 2.202 | 240 mesh 2.134; 120 μ 2.248 | equal Bτ 0.469; swapped 0.142; input E_F 1.650; 80 mesh 2.719 | ±20% |

## Oracle

22 evaluations, 44 of 50 units (88%): two electron counts per compound (80 mesh coarse, 160 mesh fine), ZrSiSe B ∥ c and B ∥ b on 160 and 80, a coarse 80-mesh scan (30, 45, 60°), the 160-mesh scan 40–70° in 5° steps, ZrSiTe B ∥ c and B ∥ b on 160. No literature value is used (24 valence electrons follow from the projections). Oracle errors < 1e-5 relative.

- **Gate (2026-10-27, cpuset 12,13,28,29):** passed; all 9 ablations fail, 3 robustness variants pass, rounding passes, deterministic; 11.0 MB, 10 hidden values, no leaks.
- **Harbor oracle** (`runs/build-zrsix-butterfly-amr/oracle-r1-zrsix-butterfly-amr`): reward 1.0, 49 s.

## Open risks

- **Mesh-80 lucky passes.** The 80 mesh passes the exponent and the ZrSiTe ratio; it is caught by MR(b)/MR(c) (−18% vs 15%) and narrowly by the peak angle (−3.2°). An agent that buys only 80-mesh data fails one field.
- **The 160 mesh itself is not fully converged**: 240-mesh values differ by +9% (MR(b)/MR(c)), +1.3° (peak), −9% (ZrSiTe). All lie inside the windows, but the windows are wider than the 160→240 change by only ~6 percentage points for MR(b)/MR(c).
- **Broad maximum.** MR(θ) varies by only 5% between 45° and 65° at Bτ = 5, so the peak angle depends on a parabola through flat data; a parabolic fit on the 10° grid of the 80-mesh coarse scan would miss.
- **Dev-easy risk.** dp1–dp3 are discoverable by reading the inputs carefully; dp6 is a standard Kohler argument. The strongest discriminators should be dp4 (the 40 mesh fakes a physically plausible saturation) and the dp1 cascade.
- **Orbit request** cost ~17 s for long orbits on 2 CPUs (returns a job id; within limits). Points with v ∥ B give tiny "closed" loops with meaningless masses.

## Disk

Peak about 0.3 GB in the cache (zip 8 MB, extracted 87 MB, models, splines 90 MB, tables 25 MB incl. the 240 check). After the build the raw downloads, extracted files, splines, histograms and tables were deleted; `fetch.py`, `members.md5`, `calibration.json` and the scripts are kept.
