# Construction of 2d-qsh-shc-plateau

This is a released-start task built from J. Zhou, S. Poncé and J.-C. Charlier, *High-throughput calculations of spin Hall conductivity in non-magnetic 2D materials*, npj 2D Mater. Appl. 9, 39 (2025) (arXiv:2501.02110). It uses Materials Cloud record rfcdy-wcy32 (DOI 10.24435/materialscloud:yv-h7, CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

1. `python3 -B fetch.py [DEST]` downloads `Data.7z` (408.8 MB, md5 1d261e78…), `README.txt` and `files_description.md` through the archive's per-file links. `Data.7z` is a single solid LZMA2 block, so no range read is possible; the whole file is downloaded and the 73 members of `members.md5` (110 MB) are extracted with the host's `7z` and verified. The default DEST is `~/.cache/ai-theorist/2d-qsh-shc-plateau`.
2. `bash run_in_image.sh build_models.py /data/x/Data/Wannier_Hamiltonians /repo/benchmarks/tasks/2d-qsh-shc-plateau/environment/lab/data insulators/Bi2_P-3m1 insulators/Hf2Br2_P-3m1 insulators/CuLi2As_P-6m2 insulators/Ti2I2N2_Pmm2 insulators/Li2Tl2_P3m1 'insulators/Bi2In2_P4|nmm' insulators/Hg4O2_P1 insulators/KSnAs_P3m1` writes the eight lab models (4.6 MB) and checks them against the archived `aiida_band.dat` (`build_models.log`).
3. `bash run_in_image.sh reference.py /repo/benchmarks/tasks/2d-qsh-shc-plateau/environment/lab/data /data/x/Data/Wannier_Hamiltonians /scratch/ref <folder>` for each material (set `QSH_SCRATCH` to a scratch folder; one CPU per material, 5–60 min each), then copy `ref/*.json` to `reference/`.
4. `bash run_in_image.sh make_task.py` stages the archive files, writes the catalog, lab data, private targets, ablation values and budgets, scans the staged files for the graded numbers and writes `calibration.json`.
5. `python3 -B write_design.py` writes `design.json`; `.venv/bin/python -m benchmarks.tools.check_task 2d-qsh-shc-plateau --cpuset 2,3,18,19` writes `check.json`; `python3 -B record_oracle.py` writes `oracle_run.json`; rerun `write_design.py` to record the oracle length.

`qsh_tb.py` (evaluator) and `shc_worker.py` (two-process spin Hall worker) are copied unchanged into the lab. `convert_hr.py` holds the hr.dat reader and Wannier90's `use_ws_distance` image selection. `z2_ti2i2n2_fine.py/.txt` is the fine Wilson-loop check of Ti2I2N2.

## What the task is

The lab evaluates the archived Wannier Hamiltonians of eight insulators of `Wannier_Hamiltonians/insulators/` with numpy and returns archived spin Hall files of `Top_materials/`.

| Request | Cost | Returns |
|---|---|---|
| `bands` | 1 | Wannier energies (optionally ⟨σ⟩) at ≤ 1000 k |
| `band_edges` | 1 | max of band n, min of band n+1 and the smallest direct gap on an N² mesh (N ≤ 300) |
| `wcc` | 1 / 2 / 3 (len(k)·n_loop ≤ 6000 / 24000 / 60300) | hybrid Wannier centres of bands 1..n |
| `shc` | ceil(k-points / 40000) + extra spin components | σ^γ_xy(μ) in Wannier90 units, with mesh, η, position convention, μ list and optional adaptive refinement |
| `ribbon` | 2 | ribbon eigenvalues and edge weights |
| `archived_shc` | 1 | one archived Top_materials Fermi scan (`aiida-shc-fermiscan.dat`) and its `aiida.win` |

The budget is 40 units. Graded: Z2 of all eight models (exact), the plateau σ^z_xy of all eight in (ħ/e)(e²/h) (± max(3 %, 0.02)), and the plateau spin direction (θ ± 2.5°, φ ± 5°) of Bi2_P-3m1 and Hg4O2_P1.

## The evaluator and its validation

- **Wigner–Seitz image selection.** The archived runs used `use_ws_distance = .true.`. Applying Wannier90's image selection (mp_grid supercell, final centres of `cut_aiida.wout`, tolerance 1e-5 Å) reproduces every archived `aiida_band.dat` to ≤ 0.07 meV; the plain hr.dat Fourier sum is off by 1–22 meV.
- **Spin ordering.** Wannier90 interleaves spinors (1↑, 1↓, 2↑, …): ‖[H(k), S_z]‖/‖H‖ is 0.03 (CuLi2As) with interleaved ordering and 1.04 with blocked ordering. Kramers partners share their centres.
- **Spin Hall formula and units.** Qiao et al. (Wannier90 `berry_task = shc`), T = 0 step occupation. CuLi2As (S_z conserved by the horizontal mirror) gives +1.0000 e/2π; the archived Top_materials sign conventions are reproduced (CuLi2As +, Hf2Br2 −, Ti2N2I2 −).
- **Position convention.** With every Wannier function at its lattice point (the hr.dat phase convention, no Berry connection) Bi2 shows σ^x_xy = 0.335, which its mirror x → −x forbids. With the diagonal position operator at the Wannier centres it shows σ^x_xy = 0.0015 and σ^y_xy = 0.352; the archived Wannier90 run (full position matrix elements, true spin, relaxed structure 2.4 % larger) gives 0.009 and 0.356.
- **Known limitations.** The spin operator is the projected spinor operator (the archive has no `.spn`); the position operator is diagonal. The MLWFs of Ti2I2N2 and Bi2In2 break the crystal symmetry, so these models give small symmetry-forbidden in-plane components (0.04–0.09), which are not graded. σ^z_xy differs between the two position conventions by < 1 % except for Hg4O2 (14 %) and Bi2In2 (23 %).

## What the archive shows (and the task is built on)

- **Fermi energies.** The scf runs use cold smearing (0.01 Ry). For six of the eight models the printed scf Fermi energy lies inside the valence band of the Wannier model (Hf2Br2 13 meV, CuLi2As 18, Bi2In2 37, Hg4O2 55, Li2Tl2 62 and Ti2I2N2 63 meV below the valence-band maximum). The `.win` `fermi_energy` equals the valence-band maximum of the scf grid in seven folders (Hg4O2: the DFT mid-gap, one of the manually corrected cases).
- **Band edges of the `.win` comments come from the coarse scf grid.** For Bi2In2 (P4/nmm) the Wannier valence-band maximum lies at a generic k (0.24, 0.94), 87 meV above the `.win` VBM, so the `.win` mid-gap (−3.896 eV) is inside the Wannier valence band (VBM −3.872, CBM −3.835 eV). The Wannier model reproduces the scf eigenvalues on the scf grid to 0.08 meV.
- **Gaps of the Wannier models:** Bi2 576, KSnAs 615, Hg4O2 59, Hf2Br2 50, CuLi2As 46, Bi2In2 37, Li2Tl2 31 and Ti2I2N2 22 meV.
- **Top_materials are relaxed structures.** Their SHC files belong to structures relaxed with scalar-relativistic pseudopotentials (Bi2 lattice 2.4 % larger; Hf2Br2 cell rotated by 60°; Ti2N2I2 relaxed to Pmmn), not to the database models; there is no hr.dat for them. Their mid-gap plateaus (0.525, −0.915, 1.026, −0.947) match the paper's Table 2 (0.5, −0.9, 1.0, −1.0).
- **Z2.** Bi2, Hf2Br2, CuLi2As, Ti2I2N2, Li2Tl2 and Hg4O2 are Z2-odd; Bi2In2 and KSnAs are trivial. For Ti2I2N2 the centre pair near 1/2 winds between k = 0.028 and 0.036 (`z2_ti2i2n2_fine.txt`); a uniform sampling of [0, 1/2] needs ≥ 201 points (151 points give 0). Li2Tl2 needs ≥ 101 points along one axis.
- **Plateaus (reference).** Bi2 0.569, Hf2Br2 −0.977, CuLi2As 1.000, Ti2I2N2 −0.852, Li2Tl2 2.786, Bi2In2 0.823 (trivial), Hg4O2 −0.758, KSnAs −0.001. Li2Tl2's plateau far above 1 (Z2-odd, large S_z mixing) and Bi2In2's large plateau in a trivial insulator show that the plateau size does not determine Z2.
- **Convergence.** Ti2I2N2: −3.84 (60²), −1.33 (120²), −0.871 (240²) unrefined, −0.852 with refinement; Li2Tl2: 3.84 (60²), 2.81 (120²). η = 5 meV changes the plateaus by ≤ 2 % (Ti2I2N2 −0.835).

## Deviations from the design (candidates.json / DETAILS.md section 28)

- **Graded values are validated-evaluator values.** The archive has hr.dat files only for the database (MC2D, unrelaxed) structures and SHC files only for the relaxed Top_materials structures, with no `.spn`. The authors' SHC files cannot be recovered from any released model, so the targets are the lab evaluator's values on the database models; the archived relaxed-structure files are offered as a purchasable trap.
- **Candidate set.** The eight insulators are the four paper QSHIs (Bi2, Hf2Br2, CuLi2As, Ti2I2N2 — in the database Ti2I2N2_Pmm2), two of the paper's other plateau materials (Li2Tl2, Bi2In2, plus Hg4O2), and the trivial KSnAs. Bi2In2 turned out Z2-trivial with a large plateau (the design's competing explanation 3).
- **Y2C2I2 replaced.** Y2C2I2 is a metal, P1 in the database, in a cell rotated with respect to the paper's C2/m structure, and its archived SHC belongs to the relaxed structure. The canted-spin question is asked instead for the in-gap plateaus of Bi2 (symmetry-constrained: σ^x_xy = 0) and Hg4O2 (P1).
- **Units.** Plateaus are graded in (ħ/e)(e²/h) as designed; tolerance max(3 %, 0.02) instead of ±5 % (the relaxed-structure values are 6–11 % away).
- **No parity fork.** The lab does not expose eigenvectors or symmetry operators, so the design's dp1 (parity vs Wilson loop) became the Wilson-loop resolution fork (dp7).
- **Edge-state count not graded.** The ribbon request is offered as a cross-check; the count depends on termination and the ribbon widths allowed for the 56–64-orbital models are narrow compared with the decay length in 22–60 meV gaps. It is a rubric item.
- **Staged excerpts.** `aiida.out` without the final band listing and `cut_aiida.wout` without the iteration lines: three unrelated log numbers coincided with graded values at 4 significant digits.

## Decision points (9) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | scf Fermi level inside the valence band; `.win` fermi_energy = VBM | metadata trap, found by comparing with the Wannier band edges |
| dp2 | `.win` band edges come from the coarse scf grid (Bi2In2) | evidence found only with a dense band-edge search |
| dp3 | k-mesh aliasing in 22–31 meV gaps | convergence, visible only after buying a second mesh or refinement |
| dp4 | lattice vs Wannier-centre position operator | symmetry cross-check (Bi2 σ^x_xy) |
| dp5 | archived Top_materials files and paper table are relaxed structures | provenance in the purchased `aiida.win` (cell) |
| dp6 | Z2 vs plateau magnitude (Bi2In2 trivial, 0.82) | competing explanation |
| dp7 | Wilson-loop resolution near small gaps (Ti2I2N2, Li2Tl2) | budgeted resolution choice |
| dp8 | vacuum-padded Wannier90 units vs sheet value | units |
| dp9 | in-plane spin components of canted plateaus | regime/symmetry |

dp1/dp2 cascade into every plateau; dp4 into both spin directions and two plateaus; dp7 and dp3 compete for budget.

## Calibration

The full table is in `design.json` under `calibration`; `calibration.json` lists every ablation's failing checks.

| Field | Reference | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| Bi2 plateau | 0.5685 | η = 5 meV 0.5685; 240² 0.5685 | relaxed archive 0.525; paper 0.5 | ±0.0171 |
| Hf2Br2 | −0.9768 | −0.9767; −0.977 | relaxed −0.915; paper −0.9 | ±0.0293 |
| CuLi2As | 1.0000 | 0.9974; 1.0000 | scf E_F 0.380; e/4π 2.0 | ±0.030 |
| Ti2I2N2 | −0.8520 | −0.835 (η 5 meV); −0.871 (240² unrefined) | relaxed −0.947; 120² −1.327; paper −1.0 | ±0.0256 |
| Li2Tl2 | 2.7857 | 2.763; 2.785 | `.win` E_F 2.400; scf E_F 0.520; 60² 3.844 | ±0.0836 |
| Bi2In2 | 0.8233 | 0.8223; 0.8234 | lattice 0.630; DFT mid-gap 0.168 | ±0.0247 |
| Hg4O2 | −0.7579 | −0.7568; −0.7579 | lattice −0.863; scf E_F −0.523 | ±0.0227 |
| KSnAs | −0.0013 | — | e/4π −0.0026 passes (floor); thickness passes | ±0.02 |
| Bi2 θ / φ (deg) | 31.78 / 89.76 | 31.77 / 89.76 | lattice 38.81 / 42.43; relaxed archive 34.11 / 88.51 (passes; fails on the plateau) | ±2.5 / ±5 |
| Hg4O2 θ / φ (deg) | 150.24 / −5.79 | 150.24 / −5.79 | lattice 154.0 / 11.1; scf E_F 176.5 / −141.7 | ±2.5 / ±5 |
| Z2 | 1,1,1,1,1,0,1,0 | 201 × 200 either axis | coarse 51 × 80 (Ti2I2N2, Li2Tl2 → 0); plateau-magnitude rule (Bi2In2 → 1) | exact |

## Oracle

See `oracle_run.json` and `design.json` → `length_estimate`. The oracle buys, per material, a 300² `band_edges`, a Wilson loop (201 × 119 for gaps below 0.1 eV, 51 × 110 otherwise) and one spin Hall calculation at mid-gap with η = 0.5 meV and the Wannier-centre convention (60² for large gaps; 120² + 8 × 8 refinement of ≤ 400 points otherwise; three components for Bi2 and Hg4O2), plus a 180² check for Ti2I2N2. No literature value is used.

Results:

- **Oracle:** 25 evaluations, 36 of 40 units (90 %), 617 s of lab jobs (slowest: the Hg4O2 three-component spin Hall run, 89 s, and the Ti2I2N2 180² check, 80 s). Largest oracle errors: Bi2In2 plateau −0.0088 (tolerance 0.0247) and Li2Tl2 −0.0040 (0.0836); angles within 0.005°.
- **Gate (2026-10-25, cpuset 2,3,18,19):** passed. All 13 ablations fail, the 3 robustness variants pass, rounding to 4 significant digits passes, grading is deterministic; 9 decision points, 8 rubric items, 5.93 MB, 12 hidden values; oracle container 642 s. After the gate the example `position` value in the instruction and catalog was replaced by a neutral placeholder and the gate was re-graded with `--reuse-state` (passed).
- **Harbor oracle (`runs/build-2d-qsh-shc-plateau/oracle-r1-2d-qsh-shc-plateau`):** reward 1.0, 7.6 min, 36 of 40 units. It ran on the package before the placeholder edit (the runner does not rerun a terminal job in the same campaign); the edit changes no file the oracle or the verifier reads.

The estimated agent length is 50–80 tool calls.

## Open risks

- **Validated-evaluator targets.** The plateaus are those of the projected-spin, diagonal-position evaluator on the database models, not of a full Wannier90 calculation with `.spn`; for non-quantized plateaus the physical value may differ by several per cent. The instruction says graded values are the service's evaluation of the archived models.
- **Position-convention fork.** Its consequence is graded through the Bi2 and Hg4O2 directions and the Hg4O2/Bi2In2 plateaus. The only direct evidence is a symmetry test on Bi2; for P1 Hg4O2 the agent must transfer that conclusion.
- **Thin margins.** Ti2I2N2 240² unrefined at η = 1 meV passes by 0.7 % of margin; the relaxed-archive Bi2 direction passes (it fails on the plateau); KSnAs's plateau is graded with an absolute floor and does not discriminate.
- **Memorisation.** The four paper QSHIs are known Z2-odd; their paper plateaus fail for Bi2, Hf2Br2 and Ti2I2N2.
- **Runtime.** A 160 000-point Hg4O2 request takes about 6 min on two lab CPUs.

## Disk

- **Peak:** about 0.8 GB: 409 MB `Data.7z`, 168 MB of non-hr extraction plus about 110 MB of hr.dat files, and about 105 MB of scratch (models, paper PDF, logs); plus the 1.15 GB gate lab image layer (removed).
- **After the build:** raw downloads, scratch files and the gate image deleted (no Harbor `2d-qsh-shc-plateau__*` images remained); `fetch.py`, `members.md5` and `reference/*.json` (the per-material reference calculations) kept.
