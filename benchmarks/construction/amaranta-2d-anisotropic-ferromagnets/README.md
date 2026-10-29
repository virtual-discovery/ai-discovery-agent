# Construction of amaranta-2d-anisotropic-ferromagnets

Released-start the benchmark task (see `benchmarks/tools/BUILDING.md`). Source: Putatunda, Cuono, Orlando and Picozzi, AMaRaNTA results, Materials Cloud record wygeh-jrc57 (DOI 10.24435/materialscloud:he-3k, CC BY 4.0); method paper arXiv:2509.20226. Provenance and checksums are in the task's `SOURCE.md`. Sibling task (same archive, triangular dihalides, no Cr compounds): `amaranta-2d-triangular-helimagnets`. Its model code was copied here and generalised; nothing is imported from it.

## Rebuild

Run everything from this folder. Compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) -e HOME=/tmp -v <repo>/benchmarks:/b -w /b/construction/amaranta-2d-anisotropic-ferromagnets ...`).

1. `python3 fetch.py [CACHE]` downloads about 57 MB into `~/.cache/ai-theorist/amaranta-2d-anisotropic-ferromagnets/`:
   - `data_readable.zip`, checked against its md5;
   - the `db.sqlite3` member of `group1`, `group4`, `group5`, `group7` and `group10.aiida`, read with HTTP range requests (central directory, local header, member). The archives are 0.9–3.5 GB and are never downloaded in full. Each member is checked against its CRC32 and sha256.
2. `python3 -B extract.py [CACHE]` writes `pool.json` (1.0 MB): structures, J_dict, couples and all 120 four-state sets (480 VASP energies) of the ten Cr compounds, with consistency asserts.
3. `python3 -B mc_reference.py full mc_CrI3_full.json` and `python3 -B mc_reference.py nn_only mc_CrI3_nn_only.json 12,18,24 20,21,22,23,24,25,26,27 20000` run the reference Monte Carlo (about 1.5 hours on 4 CPUs in total; `mc_CrI3_full.json` was accumulated in two runs, first 30×30 at 34–37 K, then 12–24 at 34–39 K).
4. `python3 -B make_task.py` (about 2 minutes) writes:
   - `environment/agent/catalog.json`;
   - `environment/lab/private/{model.json, hidden.json, hidden_values.json}`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json` in this folder.
5. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task amaranta-2d-anisotropic-ferromagnets --cpuset 2,3,18,19` writes `check.json`.

The code is split as follows:

- `model.py`: the teacher model and `validate()` against the archived DFT energies. It builds the layer space group from the full structure, symmetrises the archived tensor over the stabiliser of its bond, and maps it to every bond with K = 2A.
- `reference.py`: ground states (Luttinger–Tisza and enumeration of collinear patterns), anisotropy matrices and their single-ion/exchange split, the ideal octahedral frame and the J–K–Γ–Γ′ decomposition, DM vectors, the zone-centre gap (uniform precession, cross-checked against linearised Landau–Lifshitz dynamics of the two-site cell) and the dipolar sum.
- `mc.py`, `mc_reference.py`, `tc_analysis.py`: classical Monte Carlo (Metropolis plus over-relaxation, colour-class vectorised) and Binder-crossing analysis.
- `solution/solve.py`: an independent, agent-like re-implementation that does not import this folder.

## Archive files used (shared source group)

The task uses CrCl3, CrBr3, CrI3, CrGeTe3, CrSiSe3 and CrSBr (PKs 29719, 24863, 44376, 36128, 36985, 40899). CrSiTe3, CrPS4, CrClO and CrBrO are extracted for the validation table only. No compound of the sibling task is used.

## What the data showed (design inputs)

- **K = 2A holds for the Cr compounds too.** AMaRaNTA's single-ion formula (E₊z − E₊x − E₋x + E₋z)/4 returns half of the coefficient of (e^z)². With K = 2A the model fits the archived spin–orbit four-state energies better than with K = A for every compound; for example, CrBr3 gives 0.0031 against 0.0113 meV rms.
- **The teacher model reproduces the archived DFT.** Within-set rms residuals of the spin–orbit sets, in meV (K = A in brackets):

  | compound | rms | max | collinear J2/J3 single-flip residual |
  |---|---|---|---|
  | CrCl3 | 0.0022 (0.0042) | 0.004 | −0.17 |
  | CrBr3 | 0.0031 (0.0113) | 0.008 | 0.17 |
  | CrI3 | 0.053 (0.060) | 0.16 | 0.94 |
  | CrGeTe3 | 0.130 (0.132) | 0.24 | −3.2 to −3.6 |
  | CrSiSe3 | 0.0022 (0.0055) | 0.005 | 7.5 |
  | CrSBr | 0.040 (0.040) | 0.12 | 2.1 |
  | CrSiTe3, CrPS4, CrClO, CrBrO | 0.021–0.058 | ≤ 0.19 | −3.2 to 5.3 |

  The non-spin–orbit collinear J2/J3 calculations differ from the model's exchange sums by up to 7.5 meV (CrSiSe3), as the sibling found for the dihalides. They do not enter the task: the teacher model defines the truth.
- **The bond dependence is in the DFT.** Carrying the archived tensor to the other bonds by symmetry fits the archived energies to 0.053 meV rms for CrI3. Copying it unchanged to all bonds gives 2.4 meV (CrGeTe3: 0.13 against 0.98; CrBr3: 0.003 against 0.28). The space group is built from the full structure (D3d for the trihalides and CrGeTe3, C3v-like for the polar CrSiTe3, Pmmn-like for CrSBr). The archived tensor violates its bond's stabiliser by 0.0003–0.04 meV and is averaged over it.
- **CrCl3** has a 0.1 % monoclinic distortion: its three nearest-neighbour bonds are 3.4873 and 2 × 3.4907 Å, and the second and third shells split by 0.4–0.6 %, within AMaRaNTA's 1 % grouping. No exact operation relates the archived bond to the other two, so the tensor is carried by the pseudo-D3d operations (0.15 Å tolerance). The tensor is nearly isotropic, so this choice moves nothing graded.
- **CrSiSe3 is a Néel antiferromagnet in these data.** J1 = +5.13 and J3 = +1.85 meV are antiferromagnetic, while J2 = +0.08 meV. The archived structure straddles the cell boundary in z (Cr at z = 18.40 and 0 with c = 18.40 Å), and AMaRaNTA's 1NN and 3NN "couples" carry a −c offset (|d| = 18.7 Å). Positions are unwrapped; the four-state supercells themselves are correct. The anisotropy matrix of the Néel pattern has an easy plane (E(x) − E(z) = −0.109 meV/Cr: single-ion −0.039 plus exchange −0.070). A ferromagnet would instead have an easy axis (+0.032).
- **The anisotropy split (E(x) − E(z) per Cr of the ground state, single-ion + exchange, meV):**

  | compound | total | single-ion | exchange |
  |---|---|---|---|
  | CrCl3 | 0.0335 | 0.0279 | 0.0056 |
  | CrBr3 | 0.1595 | 0.0842 | 0.0753 |
  | CrI3 | 0.4197 | 0.2097 | 0.2100 |
  | CrGeTe3 | 0.5169 | −0.1883 | 0.7052 |
  | CrSiSe3 (Néel) | −0.1091 | −0.0387 | −0.0704 |
  | CrSBr (x = a) | 0.0628 | 0.0324 | 0.0304 |

  The exchange part grows with ligand spin–orbit coupling (Cl < Br < I; Te). In CrGeTe3 (and CrSiTe3) it overrides an easy-plane single-ion term.
- **CrI3 bond in the ideal octahedral frame:** J = −4.149, K = +1.228, Γ = 0.088, Γ′ = −0.114 meV (residual of the four-parameter form 2·10⁻⁷ meV after symmetrisation). The Cr→I directions to the upper I atoms are at −4.2°, 115.8° and 235.8°. The actual Cr–I–Cr plane normal is tilted 58.0° from c, against the ideal 54.7°. Using it instead of the ideal axis gives K = 1.250 and Γ = 0.081 (passes). The mirror-image frame gives K = −0.14 and Γ = 0.54.
- **CrSBr.**
  - The DM vector on the chain bond along a is (0, −0.264, 0) meV for the chain of site 0. The bond's mirror forbids the symmetric xz part, and the archived J_xz = +0.260, J_zx = −0.269 respect that. Inversion maps one chain onto the other, so D is staggered.
  - The archived four-state sets cannot distinguish staggered from uniform D: both fit to 0.040 meV. The staggering is a symmetry consequence, and it is graded through co- and counter-rotating cycloids.
  - The anisotropy is triaxial: easy c, intermediate a (+0.0628), hard b (+0.1176 meV/Cr). The zone-centre gap is (2/S)√(Δ1Δ2) = 0.1146 meV for S = 3/2. The uniform-precession formula equals the linearised Landau–Lifshitz frequency of the two-site cell (checked for CrI3, CrBr3, CrGeTe3 and CrSBr).
  - The model's easy axis (c) differs from the experimental b axis. The dipolar shape anisotropy of a 3 μB layer (0.13 meV/Cr favouring the plane) would move it into the plane, toward a, the model's intermediate axis. This is a rubric point, not graded.
- **Dipolar shape anisotropy (point dipoles of 3 μB, E(x) − E(z) per Cr):** CrCl3 −0.055, CrBr3 −0.046, CrI3 −0.035, CrGeTe3 −0.037, CrSiSe3 (Néel) +0.018, CrSBr −0.128 meV. Only CrCl3 (magnetocrystalline +0.0335) flips to easy-plane. The margin depends on the assumed moment, so it is a rubric item.
- **Excluded compounds.** CrPS4 has dimerised chains (3.60/3.73 Å) and an easy axis tilted 9° from c, which rests on off-diagonal elements of 0.016 meV (noise level). CrClO and CrBrO are single-ion dominated and nearly isotropic in exchange. CrSiTe3 repeats CrGeTe3.
- **`other_data/exchanges.json`** was not used; the AiiDA `J_dict` nodes and the text files agree (asserted).

## Deviations from the design (candidates.json / DETAILS.md §17)

- **Lab.** The lab is the designed teacher model: the authors' full Hamiltonian, with the tensor carried by the layer space group. It returns supercell energies of arbitrary configurations in diagonal cells of at most 36 sites, with at most 6 configurations per request. Noise σ = 0.02 meV per supercell energy, as in the sibling. Cost is ceil(N/12) per configuration; the budget is 245.
- **Compounds.** Six instead of ten (see "Excluded compounds" above). The physics question and the four competing explanations are kept. Dipolar shape anisotropy is a rubric item because the energies contain no dipolar term and the CrCl3 verdict depends on the assumed moment.
- **Graded targets.**
  - (1) The easy axis is graded as a polar angle for all six compounds, with the total anisotropy energy and its single-ion part for all six. This grades the design's central question directly.
  - (2) The Monte Carlo Curie temperature is graded for CrI3 only. Its anisotropy (0.42 meV/Cr) is four times the source precision, so Tc is recoverable within 15 %. For CrCl3 and CrBr3 it is a rubric item.
  - (3) K and Γ in the bond frame are graded for CrI3 only, with J (tolerances below). Γ′ is reported but not graded, because the defensible actual-plane-normal frame moves it by 0.047 meV.
  - (4) The Γ-point magnon gap is graded for CrSBr, the triaxial case, with the spin-length convention fixed in the prompt. For the uniaxial honeycombs it would be 4/3 of the anisotropy energy, which is already graded.
  - (5) Ten energy predictions in cells larger than the lab serves grade the full model: bond dependence, the single-ion/exchange split and the staggered DMI.
- **Kept forks:** bond dependence (dp4), symmetric vs antisymmetric off-diagonal parts (dp7), frame transformation (dp5), triaxial anisotropy (dp9), four-state configurations with SOC (dp6) and the Tc method (dp10).
- **Added forks:** CrSiSe3 being antiferromagnetic (dp1), separating single-ion from exchange anisotropy (dp2, dp3), DM aliasing and staggering (dp8), and the noise floor for small anisotropies (dp11).

## Decision points (11) and calibration

| dp | Fork | Nearest wrong path | Truth / tolerance |
|---|---|---|---|
| dp1 | CrSiSe3 assumed ferromagnetic | type FM, polar 0°, +0.032 meV | Néel, 90°, −0.109 |
| dp2 | whole anisotropy attributed to the single-ion term | single-ion 0.517 vs −0.188 (CrGeTe3), 0.420 vs 0.210 (CrI3); predictions 0.036–0.099 meV/site | ≤ max(0.005, 10 %); ≤ 0.03 |
| dp3 | exchange anisotropy neglected | CrGeTe3 in plane; totals = single-ion; predictions 0.09–0.46 | ≤ 10°; ≤ max(0.005, 10 %) |
| dp4 | CrI3 tensor copied to all bonds | polar 90°, total 1.236 meV, predictions 0.218 | ≤ 10°; ≤ 0.042 |
| dp5 | mirror-image octahedral frame | K −0.140, Γ 0.544 | K ± 20 %, Γ ± 0.06 |
| dp6 | collinear-only (no off-diagonal) mapping | K 0.544, Γ 0.316 | as above |
| dp7 | CrSBr J_xz read as symmetric | D = 0, polar 41.6°, gap 0.577 | D ± 0.05; ≤ 10°; ± 15 % |
| dp8 | CrSBr D copied to the second chain / mapped in a 2a cell | predictions 0.180 / 0.104 meV/site; D = 0 | ≤ 0.03; D ± 0.05 |
| dp9 | uniaxial gap from one in-plane axis | 0.0837 (−27 %) or 0.157 (+37 %) | 0.1146 ± 15 % |
| dp10 | Tc from mean field / nearest-neighbour model | 61.9 K / 22.0 K | 34.7 K ± 15 % |
| dp11 | single-site four-state rotations at the noise level | single-ion errors 0.017 (CrBr3), 0.013 (CrSiSe3), 0.046 (CrSBr) | ≤ max(0.005, 10 %) |

Defensible spread:

- anisotropy: oracle ≤ 0.003 meV; collective states with lab noise ≤ 0.0016 meV;
- CrI3 bond: oracle ≤ 0.002 meV; actual plane-normal frame K +1.8 %, Γ −0.007, J −0.007;
- Tc: Binder crossing of 12×12/18×18 35.8 K, susceptibility peak at 30×30 35.8 K, oracle 35.2 (+1.3 %) K;
- DM: oracle 0.004–0.010 meV;
- gap: oracle −1.2 %;
- predictions: oracle ≤ 0.0033 meV/site; unsymmetrised archived tensor ≤ 0.0035.

Tightest margins:

- predictions: tolerance 0.03 against the nearest wrong path at 0.036 (CrGeTe3, anisotropy attributed to the single-ion term). That path also fails the single-ion check by 0.70 meV;
- anisotropy single-ion: CrCl3 tolerance 0.005, against 0.0056 for the anisotropy-as-single-ion path (which fails elsewhere) and 0.0025 for the single-site path (which fails for CrBr3, CrSiSe3 and CrSBr);
- Tc: tolerance 15 % against −37 % (nearest-neighbour model, 22.0 K) and +78 % (mean field, 61.9 K) (nearest-neighbour model);
- gap: tolerance 15 % against −27 % (uniaxial, a axis). A uniaxial formula with the averaged in-plane anisotropy gives +5 % and passes.

Full numbers: `calibration.json`, and `design.json` → `calibration`.

## Oracle and length

- **Plan.** Shells by length with 1 % grouping, up to the fourth shell. The cheapest diagonal cell with at least 13 sites has no coupling reaching a site's own image and never lands +d and −d of a tensor bond on the same site: 3×3 (18 sites) for every compound.
- **Fit.** Least squares over random noncollinear configurations (10 for CrCl3, CrBr3 and CrSiSe3; 14 for CrI3 and CrGeTe3; 16 for CrSBr) plus ferromagnetic x, y, z and partner-pattern x, z states, with these parameters:
  - for the honeycombs, one tensor carried by C3 about Cr, after checking C3 on the ligands;
  - for CrSBr, free tensors for both chains, then a refit with the tie the free fit reveals (T₁ ≈ T₀ᵀ). CrSBr also adds six off-axis ferromagnetic states, from which its anisotropy matrix and gap are taken directly.
- **Tc.** Coarse 12×12 scan from 10 to 80 K, then the Binder crossing of 12×12 and 18×18 on a 7-point grid. Monte Carlo uses Metropolis plus over-relaxation on up to 4 processes.
- **Usage.** 22 evaluations, 220 of 245 units (89.8 %).
- **Monte Carlo crossings (reference).** The Binder crossings drift with size: 12/18 35.8, 12/24 35.4, 12/30 35.2, 18/24 34.8, 18/30 34.8 and 24/30 34.7 K. The susceptibility peaks are 36.7 (12), 36.5 (18), 36.2 (24) and 35.8 K (30). The graded reference is the 24/30 crossing, 34.72 K.
- **Time.** about 15.5 minutes (927 s) of container time in check_task on 4 contended CPUs, most of it the oracle's own Monte Carlo; the Harbor oracle took 14 min 20 s with reward 1.0.
- **Expected agent length.** 55–85 tool calls.

## Open risks

- **Dev-easy risk.** A careful agent that inspects the structures, fits tensors by regression with symmetry, uses collective states for the anisotropy, cross-checks the second CrSBr chain and writes a correct Monte Carlo solves the task. The discovered-evidence forks are:
  - CrSiSe3 being antiferromagnetic;
  - the noise floor of single-site rotations;
  - CrSBr's 2a aliasing and DM staggering;
  - the sign change of CrGeTe3's single-ion term.
- **Teacher-model truth.** The CrSBr DM staggering follows from symmetry, not from the archived energies, which cannot distinguish it. The CrCl3 tensor relies on pseudo-symmetry. The anisotropies of CrCl3 and CrSBr (0.03–0.12 meV) are near the source precision. CrGeTe3 is reproduced worst (0.13 meV rms within sets).
- **Monte Carlo.** The graded Tc comes from Metropolis plus over-relaxation chains of 30,000 sweeps at 24×24 and 30×30. Binder-crossing statistical error is about 0.3–0.4 K, well inside the 15 % tolerance. Agents with short or poorly equilibrated runs may land outside it.
- **Category and convention edges.**
  - The polar angle is folded to 0–90°.
  - The ground-state names must match the list in the prompt.
  - The gap convention (S = 3/2, classical single-ion treatment) is stated explicitly. An agent using the quantum (2S − 1)/2S factor for the single-ion term would get 0.099 meV (−13 %), inside the tolerance.

## Disk

- Peak use was about 57 MB in the cache: the zip plus five `db.sqlite3` members of 4.6–21 MB. The raw downloads were deleted after the gate; `fetch.py` restores them with pinned checksums.
- The task package is 0.25 MB (0.15 MB counted by the gate); this folder is 1.2 MB, including `pool.json` (1.0 MB) and the Monte Carlo results.

## Gate

- **check_task (2026-10-20, 17:14 UTC):** passed.
  - 11 decision points; all 13 ablations fail.
  - 5 robustness variants pass; rounding to 4 significant digits passes.
  - 29 hidden values, no leaks.
  - Oracle: 22 evaluations, 220 of 245 units.
- **Harbor oracle:** reward 1.0 (`runs/build-amaranta-2d-anisotropic-ferromagnets/oracle-r1-amaranta-2d-anisotropic-ferromagnets`). The CrI3 Tc came out at +1.3 % from the reference.
