# Construction of cu2wse4-orbital-edelstein

This is a released-start task. It is built from Nakazawa, Yamaguchi and Yamakage, *Orbital-related gyrotropic responses in Cu2WSe4 and chirality indicator*, Phys. Rev. Materials 10, 064603 (2026) (arXiv:2512.16387). It uses Zenodo record 20131063 (CC BY 4.0), file `CWS_PRM.zip`. Provenance and checksums are in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

1. `python3 -B fetch.py` downloads `CWS_PRM.zip` (24.6 MB) and verifies its 34 members into `~/.cache/ai-theorist/cu2wse4-orbital-edelstein/x`.
2. `bash run_in_image.sh make_model.py` converts the two `.HWR` files to npz.
3. `bash run_in_image.sh make_cache.py wosoc 16 24 32 48 64 96 128`, then the same with `soc`, writes the histogram caches with `CWS_DEGEN = 5e-3`, the default. On cpuset 6,7,22,23 under shared load:
   - without SOC: about 11 min;
   - with SOC: about 90 min, of which the 128 mesh is 55 min.
4. `bash run_in_image.sh make_task.py` (about 2 min) writes the staged files, the compact lab data, the catalog, the private targets, the ablation values, both budgets and `calibration.json`.
5. `.venv/bin/python -m benchmarks.tools.check_task cu2wse4-orbital-edelstein --cpuset 6,7,22,23`.
6. `python3 -B record_oracle.py`, then `python3 -B write_design.py`.

`cws_tb.py` is the evaluator. It is copied unchanged into the lab and is also the construction engine.

## What the archive shows

- **Models.**
  - Both Hamiltonians are in Hartree, on a 7 × 7 × 3 Wigner–Seitz box. The weights sum to 108, so the Wannier grid was 6 × 6 × 3.
  - Basis: Cu-d ×4, W-d ×2 and Se-p ×8, in CIF order with OpenMX real harmonics. The SOC model is spin-block ordered (uu = dd* exactly).
  - The two header Fermi levels (−0.180300 and −0.180817 Ha) lie mid-gap. Relative to them:

    | Model | VBM (eV) | CBM (eV) |
    |---|---|---|
    | SOC | −0.5988 (0.5 meV above M, on the M–Γ line) | +0.6437 (off Γ along Γ–X, 1.5 meV below Γ) |
    | no SOC | −0.6445 (at M) | +0.6776 (at Γ) |

  - The Wannier bands reproduce the archived DFT bands (Fig1b; path Γ–X–M–Γ–Z–R–A–Z|X–R|M–A, k in Bohr⁻¹) to 2.5 meV near the gap.
- **Symmetry breaking off the grid.** S4 breaks by up to 2 meV; σd(xy) is exact. Folding degeneracies are split off-grid by up to 3.3 meV (no SOC) or 0.35 meV (SOC), at Z, R and A and along Γ–Z and X–M. With a 1e-6 eV degeneracy cutoff these splittings produce huge Berry-curvature-dipole spikes and a 7–14% bias in the near-edge orbital Edelstein response. A 5 meV cutoff (degenerate-subspace treatment) removes both and reproduces the archive best.
- **Edelstein validation.** Fig2 is in units of 10⁻⁵ μ_B/f.u., and 1 A/m equals 1.8689 × 10⁻⁵ μ_B/f.u. At kT = 30 meV and γ = ħ/2τ = 30 meV:
  - the spin curve matches the evaluator in absolute units (scale 1.002, rms residual 2–6%);
  - the orbital curves match the shape to 2–3% rms, but with a factor −2.0 (×1.01–1.02, no energy shift); the no-SOC curve includes the factor 2 for spin.
- **The orbital sign question.** The evaluator's orbital moment is the paper's Eq. (3): m = (|e|/2ħ) Im⟨∇u|×(H−ε)|∇u⟩. It equals the direct wave-packet moment (q/2)⟨r × v⟩, q = −|e|, exactly on a model triangle molecule, in sign and magnitude. The archived orbital curves therefore carry a factor −2 relative to the paper's own formula. The archive's orbital and spin curves have the same sign in the valence band; physically they have opposite signs. Graded values use the evaluator under the conventions stated in the instruction.
- **Berry-curvature dipole.** With the 5 meV cutoff, Fig4c,d are reproduced in shape (residual about 0.3 at kT = 10 meV on meshes ≤ 64), with a scale of about 1.92–1.95. The likely cause is a factor 2 in the definition of σ_ijl for j ≠ l (counted once versus both orders) or in τ. BCD magnitudes are not graded.
- **Intra-atomic moment.** −μ_B⟨L⟩ is validated by fitting the SOC model's own on-site blocks to λL·S (correct signs on Cu, W and Se). It gives an orbital-to-spin ratio of +0.27 at the valence edge and −2.41 at the conduction edge. The Berry-phase moment gives −11.2 and −0.34, so the intra-atomic operator reverses both conclusions.

## Deviations from the design (candidates.json / DETAILS.md section 34)

- **Graded values come from the validated evaluator, not from the archive,** because of the archive's orbital factor of −2. The archive validates the shapes and the spin units.
- **The SOC/no-SOC comparison is graded at equal hole density** (5 × 10²⁰ cm⁻³), not as "agreement within 20%". This adds the spinless-model counting cascade (dp3, dp4) and the energy-axis trap (dp5).
  - The value is 0.786: SOC reduces the valence-band orbital response by about 20%.
  - On the paper's common header axis the ratio is 1.07 ("agree well").
- **Nonlinear Hall components are graded as exact symmetry content.** The allowed set is graded, and σ_zxy/σ_xyz = −2 for the BCD. "σ_xyz, σ_zxy magnitudes ±30%" is dropped (see above). "BCD dominance over other terms" cannot be computed from the released models and is a rubric item.
- **Added fields:** the absolute orbital magnetization in A/m (units, τ and sign), α_yy/α_xx and the chirality indicator (both exact by symmetry).
- **NCTE Hall (Fig3) is not used.**

## Decision points (10)

| dp | Fork | Kind |
|---|---|---|
| dp1 | P1 CIF: find D2d with the twofold axes along a and b; not the rotated −4m2 setting or a chiral class | symmetry, discovered from the coordinates |
| dp2 | Berry-phase (wave-packet) versus intra-atomic orbital moment | model choice; reverses both conclusions |
| dp3 | The spinless no-SOC model counts twice in α | metadata (54 functions, collinear) |
| dp4 | The spinless model counts twice in the hole density (cascade into μ(p)) | cascade; the two errors partly cancel |
| dp5 | Header Fermi levels mid-gap at different distances from the band edges | discovered evidence; the paper's own comparison axis |
| dp6 | The VB-edge spin term converges slowly, while M_orb converges by N = 24 | convergence under budget |
| dp7 | τ = ħ/(2γ), q = −|e|, per-volume normalisation | units, sign |
| dp8 | CB spin dominance (W-d SOC) | transfer |
| dp9 | σ_zxy = −2σ_xyz for the BCD, not permutation symmetry | symmetry plus mechanism |
| dp10 | tr G = 0 (achiral, gyrotropic) | classification |

## Calibration

| Field | Reference (mesh 128) | Defensible (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| orb/spin, VBM − 0.10 | −11.19 | 96: −10.60; 1/N extrapolation: −13.40 | 32: −7.60; 48: −9.00 (passes); intra-atomic: +0.27 | ±25% |
| orb/spin, CBM + 0.10 | −0.339 | 96: −0.337; all meshes −0.32 to −0.37 | from VB: −11.2; intra-atomic: −2.41 | ±25% |
| M_x^orb (A/m) | −0.09412 | 96: −0.09402; 16–128: −0.0905 to −0.0941 | τ = ħ/γ: −0.188; q > 0: +0.094; intra-atomic: +0.0022 | ±15% |
| SOC/no-SOC at p = 5e20 | 0.786 | 96: 0.785; 32: 0.787 | both spin factors omitted: 0.928; header axis: 1.066; intra-atomic: 0.562; holes only: 0.464; α only: 1.572 | ±10% |
| symmetry fields | [xx, yy], −1, 0; 6 NLH components, −2 | raw mesh-128 BCD ratio: −2.000005 | rotated, chiral, permutation-symmetric variants | exact / ±0.05 / ±0.15 |

## Oracle

The oracle makes 6 evaluations for 26 of 30 units (87%):

- one `bands` request per model (path plus local grids around M and Γ);
- `fermi_surface` at meshes 96 and 128 for both models.

It derives the point group from the CIF by brute force, and uses no literature values. Its errors are below 0.03%.

- **Gate (2026-10-26, cpuset 6,7,22,23):** passed. All 13 ablations fail, 3 robustness variants pass, rounding passes and grading is deterministic. The package is 13.0 MB.
- **Harbor oracle:** reward 1.0, runtime 49 s.

## Open risks

- **Low oracle evaluation count.** 6 evaluations is below the 15–60 target. Each mesh request is expensive; agents will spend more on convergence studies and k-list checks.
- **Dev-easy risk.** dp1, dp7, dp9 and dp10 are textbook. The likely discriminators are:
  - dp2: the operator, with both offered;
  - dp3–dp5: the spinless-model cascade and the energy axis;
  - dp6: the spin convergence, where the 48 mesh lucky-passes at −19.6%.
- **Mesh uncertainty in the VB ratio.** It is not fully converged at 128 (≈ −11 to −13). The ±25% window covers the extrapolation; the 32 mesh fails at −32%.
- **Evaluator conventions.** The 5 meV degeneracy cutoff and the atomic gauge are fixed by the lab and stated in the catalog. A different cutoff (1e-6 eV) would move the VB orbital values by about 10%. That cutoff is not exposed to the agent.
- **Archive conflict.** An agent that remembers the paper's Fig. 2 (orbital and spin of the same sign in the VB) could be misled. The graded sign follows the paper's own Eq. (3) and the instruction's stated convention.

## Disk

Peak use was about 0.25 GB:

- zip: 25 MB;
- extracted: 113 MB;
- model npz: 30 MB;
- caches: 35 MB, including the 1e-6 and test-cutoff caches;
- scratch.

After the build, raw downloads, extracted members and caches were deleted. The refetch script, `members.md5` and `calibration.json` are kept.

## Fix after the dev pilots (2026-10-26)

Sonnet 5 and GPT-5.6 Terra both reported orbital_to_spin_ratio["CBM+0.10"] ≈ −2.77 against the reference −0.339. Both used the **sum** of the Berry-phase and intra-atomic moments. Reproduced through the lab's caches, this path gives:

| Field | Summed moments | Reference |
|---|---|---|
| CB ratio | −2.7507 (Terra: −2.75067) | −0.339 |
| VB ratio | −10.93 | −11.19 |
| M_x (A/m) | −0.0919 | −0.0941 |
| SOC/no-SOC ratio | 0.794 | 0.786 |

The summed path therefore passes every graded field except the CB ratio.

The old wording ("Wannier functions sit on their atoms") did not exclude the summed moment. The rubric itself called the sum the "complete" moment. So the instruction was ambiguous.

The instruction now reads: "r is the position operator of the tight-binding model, diagonal in the Wannier basis, with each Wannier function at the position of its atom. The orbital magnetization is the response of this moment alone; do not add any other orbital contribution."

- Rubric item 3 now treats the intra-atomic part as a stated limitation.
- The summed path is ablation `berry_plus_intra_atomic_moment` (dp2): computed by `add_summed_ablation.py` from the lab caches, and included in `make_task.py`.
- References and tolerances are unchanged.
- Gate re-run on cpuset 14,15,30,31: passed, with all 14 ablations failing.
- Harbor oracle in campaign `build-cu2wse4-orbital-edelstein-fix`: reward 1.0.
