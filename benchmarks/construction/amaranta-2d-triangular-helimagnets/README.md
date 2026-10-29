# Construction of amaranta-2d-triangular-helimagnets

Released-start the benchmark task (see `benchmarks/tools/BUILDING.md`). Source: Putatunda, Cuono, Orlando and Picozzi, AMaRaNTA results, Materials Cloud record wygeh-jrc57 (DOI 10.24435/materialscloud:he-3k, CC BY 4.0); method paper arXiv:2509.20226. Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

Run everything from this folder. Compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) -e HOME=/tmp -v <repo>/benchmarks:/b -w /b/construction/amaranta-2d-triangular-helimagnets ...`).

1. `python3 fetch.py [CACHE]` downloads about 7 MB into `~/.cache/ai-theorist/amaranta-2d-triangular-helimagnets/`:
   - `data_readable.zip`, checked against its md5;
   - the `db.sqlite3` member of `group2.aiida`, `group3.aiida` and `group10.aiida`, read with HTTP range requests (central directory, local header, member). The archives are 0.9–4.5 GB and are never downloaded in full. Each member is checked against its CRC32 and sha256.
2. `python3 -B extract.py [CACHE]` writes `pool.json` (0.7 MB): structures, J_dict, couples and all 96 four-state sets (384 VASP energies) of the eight compounds, with consistency asserts.
3. `python3 -B make_task.py` (about 4 minutes on 4 logical CPUs) writes:
   - `environment/agent/catalog.json`;
   - `environment/lab/private/{model.json, hidden.json, hidden_values.json}`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json` in this folder.
4. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task amaranta-2d-triangular-helimagnets --cpuset 2,3,18,19` writes `check.json`.

The code is split as follows:

- `model.py`: the teacher model (site point group from the full structure, bond tensors, K = 2A) and `validate()` against the archived DFT energies.
- `reference.py`: bond tables, the exact supercell energy used by the lab, Luttinger–Tisza, equivalence classes of q, classification and the easy-axis search.
- `solution/solve.py`: an independent, agent-like re-implementation. It does not import this folder.

## Archive files used (shared source group)

This task uses only the eight triangular dihalides NiCl2, NiBr2, VCl2, VBr2, VI2, MnCl2, MnBr2 and MnI2, at PKs 11815, 12821, 16233, 17997, 14160, 4157, 3579 and 2705. It draws on:

- their `*.J_dict.txt` and `*.POSCAR.vasp`, plus `other_data/PK_list.txt`, from `data_readable.zip`;
- the `db.sqlite3` provenance databases of `group2` (V and Ni dihalides), `group3` (Mn dihalides) and `group10` (NiCl2; the CrSiTe3 workchain in it is not read).

No Cr compound, no NiI2 and none of the `repo/` VASP output blobs are used. The sibling catalog task `amaranta-2d-anisotropic-ferromagnets` (Cr compounds) can use group1, group4, group5, group7 and group10 without overlapping in materials or targets.

## What the data showed (design inputs)

- **The readable exchanges.json is corrupted.** Its parameters are carried over between entries; for example, the MnBr2 entry holds the MnI2 values, and the NiBr2, NiCl2, VCl2 and VBr2 entries hold others. The `*.J_dict.txt` files and the AiiDA `J_dict` nodes agree (asserted), so they are the source. The zip has no NiCl2 POSCAR; the structure comes from the AiiDA workchain input.
- **The pair geometry is in the AiiDA "couples" output.**
  - The 1NN tensor belongs to a specific bond, for example a2 in the hexagonal cells.
  - AMaRaNTA groups neighbour distances within 1 %. In the distorted rectangular cells of VBr2, VI2 and MnCl2, its "2NN" is therefore the second, longer nearest-neighbour bond (3.769, 4.071 and 3.735 Å; four plus two bonds) and its "3NN" is the bond along a. J2 is 7.78, 1.06 and 1.80 meV, against J1 = 12.27, 4.00 and 1.19 meV.
  - MnBr2 and MnI2 have rectangular cells but are undistorted.
- **The teacher model reproduces the archived spin-orbit four-state energies.** Within each set, the rms residual is 0.003–0.03 meV (VI2: 0.26 meV), provided that two conditions hold:
  - the tensor is carried to the other bonds with the true site point group (D3d or C2h, axial vectors). Using C6 rotations, which are not symmetries of the 1T layer, gives errors of meV;
  - the single-ion coefficient is K = 2A. AMaRaNTA's formula (E₊z − E₊x − E₋x + E₋z)/4 returns half the coefficient of (eᶻ)², and K = 2A fits better than K = A.
- **The collinear J2/J3 calculations disagree with the model's exchange sum.** Their single-flip energies differ from the model's ΣJ by up to 8 meV (VI2) or 22 %, which points to SOC-free J1 differences or couplings beyond the three shells. This does not affect the task, because the teacher model defines the truth, but it bounds how literally the parameters can be read.
- **Ground states (Luttinger–Tisza) of the model:**
  - NiCl2 and NiBr2: spirals along the bond direction, q = 0.0811 and 0.1192 of (b1 + b2);
  - VCl2, MnBr2, MnI2: 120°. MnBr2 and MnI2 would be stripes in a J1–J2 model (J2/J1 = 0.15 and 0.37), but the antiferromagnetic J3 stabilizes K;
  - VBr2 and MnCl2: incommensurate spirals along b (k = 0.789, 0.608), with cos(πk) = −J'/(2J);
  - VI2: collinear stripe (J' > 2J).
  - The anisotropic exchange shifts q by at most 0.006 Å⁻¹.
- **Anisotropy.**
  - VI2's stripe moments are tilted 30.7° from c toward a. Along the Cartesian axes E_z < E_y < E_x, so a three-axis comparison or the single-ion term (easy axis c) misses the tilt; the tilt comes from J_xz, which survives the mirror average.
  - The 120° plane is vertical in VCl2 (by 0.031 meV/site) and MnBr2 (by 0.012). In MnI2 it is horizontal by only 0.003 meV/site, so it is not graded.

## Deviations from the design (candidates.json / DETAILS.md §16)

- **Lab.** The lab is the designed teacher model: the authors' full Hamiltonian, with the tensor carried by the site group (C2h for the distorted cells, not the C3 the design assumed). It returns supercell energies of arbitrary configurations; cells are diagonal, with at most 36 sites and at most 6 configurations per request. Noise σ = 0.02 meV per supercell energy, which matches the scatter of the archived energies about the model. Cost is ceil(N/12) per configuration; the budget is 240.
- **Graded compounds.** All eight compounds are graded on type and q, instead of three held-out ones.
- **No Monte Carlo ordering temperature and no Γ magnon gap.** With isotropic exchange there is no finite-T order in 2D. Tc and the gap would be set, logarithmically for Tc, by anisotropies of 0.01–0.1 meV, below the source's stated precision. For the spirals the gap is ill-defined: with anisotropic exchange an incommensurate spiral has no simple LSWT. Neither is uniquely recoverable, so finite-temperature physics is a rubric item.
- **Anisotropy fork.** It is graded through the VI2 easy axis.
- **Energy predictions added.** Each compound has three configurations in 7×7 or 4×7 cells, larger than the lab serves, so that exchange magnitudes are graded where the ground state is insensitive (aliased J3 in Mn and VCl2, pair counting).
- **Dropped forks.** "Sign convention" was dropped because the agent maps energies and no sign is graded. "S vs unit vector" is fixed by the prompt (classical unit vectors). "Noise-limited anisotropy" is part of dp5: the tilt is resolved cheaply by global stripe energies, not by single four-state components.

## Decision points (7) and calibration

| dp | Fork | Nearest wrong path | Truth / tolerance |
|---|---|---|---|
| dp1 | six NN equivalent in distorted cells | 120° for VBr2, VI2, MnCl2; q error 0.10–0.51 Å⁻¹; predictions 0.43–5.6 meV/site | type exact; q ≤ 0.02 Å⁻¹ |
| dp2 | 2a pair in a 4-periodic cell (J3 doubled) | Ni q error 0.074, 0.146 Å⁻¹; predictions 0.89 (VCl2), 1.5, 2.3 (Mn), 6–11 (Ni) | q ≤ 0.02; RMSE ≤ 0.3 |
| dp3 | exchange truncated after J2 | Ni → FM; MnBr2, MnI2 → stripe | type exact |
| dp4 | compare FM, 120° and stripe only | Ni → FM, VBr2 → stripe, MnCl2 → 120° | type exact |
| dp5 | VI2 axis from x, y, z or SIA | c axis, 30.7° off | ≤ 7° |
| dp6 | VI2 tensor copied to all bonds | 39.9° off | ≤ 7° |
| dp7 | ordered pairs double counted in predictions | RMSE 3.6–32 meV/site | ≤ 0.3 |

Defensible spread:

- q: oracle ≤ 0.0026 Å⁻¹; anisotropic LT ≤ 0.0062; collinear J_zz mapping ≤ 0.0004.
- Predictions: oracle 0.013–0.144 meV/site. The largest is VBr2, because shared-reference four-state sets give a coherent error of about 0.08 meV/site (1σ) for six classes. The exact isotropic model is within 0.061.
- VI2 axis: oracle 2.8°; six noisy probes in a 3×3 cell 1.2°.

Tightest margins:

- predictions: VBr2 oracle 0.144 against tolerance 0.3, with the nearest wrong path at 0.43 (MnCl2, ideal lattice). The small J3 of the distorted compounds (0.07–0.17 meV) is deliberately below what this check resolves;
- q: tolerance 0.02 against 0.074 (NiBr2 with J3 doubled), with defensible spread ≤ 0.006;
- axis: tolerance 7° against 30.7°, with defensible spread ≤ 2.8°.

Full numbers: `calibration.json`, and `design.json` → `calibration`.

## Oracle and length

- **Plan.** Bond classes by length up to 2.05 × the NN distance. For each compound, the cheapest diagonal supercell in which every class has a representative without periodic images within that range:
  - 3×8 for the hexagonal cells;
  - 3×6 for the distorted rectangular cells (6 classes);
  - 2×6 for MnBr2 and MnI2.
- **Mapping.** Collinear four-state configurations sharing one reference spin, then Luttinger–Tisza with classification.
- **VI2 axis.** Nine stripe directions, off-axis probes included, in a 2×3 cell, then a fit of the quadratic form.
- **Usage.** 21 evaluations, 215 of 240 units (90 %). Mapping every class of all compounds in 36-site cells would cost about 250–290 units.
- **Time.** About 170 s of container time. The Harbor oracle took 2 min 44 s with reward 1.0.
- **Expected agent length.** 50–80 tool calls.

## Open risks

- **Dev-easy risk.** A careful agent that groups bonds by exact length, uses alias-free cells and probes off-axis directions solves the task with standard tools. The discovered-evidence forks are the 1.5 % distortion, which only a cross-check shows to matter (40–70 % different exchange), and the VI2 tilt, which a three-axis comparison misses.
- **Teacher-model truth.** The VI2 tilt and the anisotropies are model consequences below the source's 0.1 meV precision, and VI2 is the compound the model reproduces worst (0.26 meV rms within sets). The ground states depend only on the isotropic exchange and are robust.
- **Single-q assumption.** The ground states of NiBr2 and NiCl2 are taken as single-q spirals. Their anisotropic exchange (≤ 0.1 meV) is ten times smaller than in NiI2, for which a multi-q state has been proposed; this was not checked by Monte Carlo.
- **Category and tolerance edges.** The type check is exact and the q check uses symmetry equivalents. An agent with noisy J that lands MnI2 or MnBr2 slightly off K must still report "120-degree".

## Disk

- Peak use was about 110 MB in the cache: duplicate exploration copies of the three `db.sqlite3` members (55 MB in total) and the zip. The raw downloads were deleted after the gate; `fetch.py` restores them with pinned checksums.
- The task package is 0.46 MB; this folder is 0.8 MB, including `pool.json`.

## Gate

- **check_task (2026-10-20):** passed. 7 decision points; all 7 ablations fail; 4 robustness variants pass; rounding to 4 significant digits passes; 36 hidden values, no leaks.
- **Harbor oracle:** reward 1.0 (`runs/build-amaranta-2d-triangular-helimagnets/oracle-r1-amaranta-2d-triangular-helimagnets`).
