# Construction of inmnti2-weyl-census

This is a released-start task. It is built from Grassano, Binci and Marzari, *Type-I antiferromagnetic Weyl semimetal InMnTi2*, Phys. Rev. Research 6, 013140 (2024) (arXiv:2208.11412). It uses Materials Cloud record hkprp-xw234 (DOI 10.24435/materialscloud:ph-3c, CC BY 4.0), file `upload.zip`. Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. Run `python3 -B fetch.py [DEST]`. It downloads `upload.zip` (160.5 MB, md5 27d1aa76…) and `files_description.md` through the archive's per-file links. It then extracts and verifies the 96 members listed in `members.md5` and streams excerpts (header and final state) of the two Wannier90 `.wout` logs. The default DEST is `~/.cache/ai-theorist/inmnti2-weyl-census`. The extraction uses about 205 MB besides the zip.
2. Run `bash run_in_image.sh make_task.py` in the runtime image on cpuset 6,7,22,23. It takes about 4 min; `--reuse-census` reuses `census.json` and takes 30 s. It does the following:
   - converts both `InMnTi2_hr.dat` files to `lab/data/*.npz` (complex64; band energies change by less than 0.001 meV);
   - writes `lab/data/archive.json` (DFT bands, Kubo rows with ω ≤ 1.5 eV, WannierTools outputs);
   - copies `tbmodel.py` into the lab;
   - stages 21 archive files plus excerpts under `agent/archive/` and writes `catalog.json`;
   - writes `private/hidden.json`, `hidden_values.json`, `solution/ablation_values.json`, the budget in both `public.json` files, and `census.json` and `calibration.json` here.

   It asserts the archive facts below. `model.py` holds the readers and estimators. All references use the lab's own evaluator on the lab's npz files.
3. Run `.venv/bin/python -m benchmarks.tools.check_task inmnti2-weyl-census --cpuset 6,7,22,23`. It writes `check.json`.
4. Run `python3 -B record_oracle.py` (writes `oracle_run.json`), then `python3 -B write_design.py` (writes `design.json`).

## What the task is

The lab holds the two archived Wannier Hamiltonians and evaluates them with numpy (`tbmodel.py`). It also serves archived DFT bands, Kubo outputs and WannierTools outputs.

| Request | Cost | Returns |
|---|---|---|
| `wannier_bands` | 1 | Wannier band energies at ≤ 500 k |
| `gap_scan` | 1 (mesh ≤ 24) or 2 (≤ 36) | local minima of E_{n+1} − E_n on a mesh³ grid |
| `find_nodes` | 1 | BFGS minimisation of the gap² from ≤ 8 seeds (Hellmann–Feynman gradient) |
| `berry_flux` | 1 | link-variable Berry flux of bands 1..n through ≤ 8 spheres |
| `dft_bands` | 3 | the band energies of one of `hubbard/nscf`, `hubbard/nscfXYZ`, `standard/nscf`, `standard/nscfXYZ`, with input and header lines |
| `kubo` | 2 | archived Kubo files of one folder `10`…`300` |
| `wannier_tools` | 2 | archived node or chirality outputs |

The budget is 28 units. The archived items alone are worth 32 units, plus any number of evaluations. My evaluator reproduces the archived WannierTools node positions and energies, and all 24 chiralities.

The graded fields are:

- the node counts (DFT+U and PBE);
- the sizes of the magnetic symmetry orbits;
- the energy, pair separation and k_z offset of the near-k_z = 0 family;
- the PBE node energy;
- σ_xy, σ_yz and σ_zx.

## What the archive shows (and the task is built on)

- **`standard` is non-magnetic.** Its scf input has no `starting_magnetization`, and its output has no magnetization, 24 symmetry operations, QE 6.5 and E_F = 15.6841 eV. It is PBE+SOC with time reversal. Its 24 nodes (bands 52/53 = Wannier 18/19) lie at +38.4 to +39.2 meV. Neighbouring opposite-chirality nodes are only 0.155 Å⁻¹ apart, so grids up to 24³ show 12 minima. Six further minima keep a 6.5–6.9 meV gap; they are avoided crossings. The PBE Wannier model matches its DFT lines to 1.7 meV and is symmetric to 3.5 meV.
- **`hubbard` is the magnetic state.** Ti moments are +2.0 and +1.9 μB, Mn −3.6 μB, along z (net 0.11 μB/cell); E_F = 16.3261 eV, QE 7.0, U(Mn) = 2.45 eV and U(Ti) = 2.27 eV. Its scf keeps 8 operations, 4 of them combined with time reversal (the magnetic point group of T_d with m ∥ z).
- **Band index.** `exclude_bands = 1-34` in both .win files. The DFT+U disentanglement window (`dis_win_min` 6.33 eV) removes two more bands, so DFT band 52 is Wannier band 16 for hubbard but band 18 for standard. WannierTools used `NumOccupied` = 16 and 18. Using offset 34 for hubbard gives 72 "nodes" between Wannier bands 18/19, at +0.66 to +1.44 eV.
- **The DFT+U Wannier model fails off its grid.** It was built on a 9×9×9 grid (the paper says 12×12×12, which is true only for `standard`), with random projections and spreads of 4–7 Å².
  - It reproduces DFT to 0.4 meV at the scf/nscf k-points, which belong to its fitting grid.
  - Between grid points it breaks the magnetic symmetry: up to about 0.1 eV, RMS 13–15 meV near E_F, even for C2z.
  - Its 24 nodes (bands 16/17, total chirality 0) scatter from −32 to −107 meV. Within the 8-node orbit near k_z = 0 they span −47.8 to −96.4 meV.
  - On the archived DFT+U lines (`hubbard/nscfXYZ`), DFT crosses (minimum gap 0.3 meV on x) at −97.1 to −97.8 meV for both nodes of the pair, symmetrically. The Wannier model leaves gaps of 9–26 meV on the same lines.
  - The paper's statement of "96 meV, all nodes equivalent within 1 meV" holds for DFT, not for the archived model. The paper's kW = (0.2361, 0.3712, 0) is in units of 2π/3.134 Å and coincides with one Wannier node.
- **Magnetic orbits.** Antiunitary operations act as k → −Rk. Nodes on the k_x = 0 and k_y = 0 planes have a C2x·T (C2y·T) stabiliser, so the 24 nodes form orbits of [4, 4, 4, 4, 8]. The 8-node orbit sits at k_z = ±δ, and the sign is tied to chirality: χ = +1 at k_z < 0.
- **Kubo folders.** The folder names 10, 30, 58, 120, 210 and 300 are k-mesh sizes. The single `kubo/InMnTi2.win` shows `berry_kmesh = 300 300 300` and `kubo_smr_fixed_en_width = 0.002`; the candidate design itself misread the names as broadenings.
  - σ_xy(0) by folder is 74.4, 108.8, 103.0, 110.9, 114.4 and 114.1 S/cm (non-monotonic).
  - σ_yz and σ_zx are archived at −0.20 and −4.46 S/cm (300³), but symmetry forbids them for m ∥ z.
- **Fermi energy of the standard line run.** The `standard/nscfXYZ` output prints E_F = 15.7968 eV, computed from its line k-points. The `standard/nscfXYZ` lines do not pass through the PBE nodes.
- **Lattice.** All electronic structures use the PBE/screening cell (3.1336 Å primitive), where the DFT+U pressure is 237.6 kbar. The archived DFT+U vc-relax converges at 3.364 Å, with V = 76.13 Å³ (+24 %), 0.9 eV lower. The paper says the relaxation reproduces 3.134 Å. The relax output ends before its final scf.

## Deviations from the design (candidates.json / DETAILS.md section 26)

- **Targets changed.**
  - The node energy, separation and k_z offset come from the archived DFT lines, not from the Wannier model, because the model breaks the symmetry. Only the family containing the near-k_z = 0 nodes is graded, because the other four orbits have no DFT check.
  - "12 pairs within ±0.2 eV of E_F" became a count of nodes between DFT bands 52/53. Within ±0.2 eV, Wannier bands 15/16 also cross near Γ, so a window definition is ambiguous.
  - Added the orbit sizes, the PBE node energy and the forbidden Hall components.
  - "Representative node in reduced coordinates" was dropped. The paper's kW is a Wannier node, about 0.008 Å⁻¹ from the DFT node, and the separation captures the same error.
  - Type I/II was dropped as a graded field. The tilt-to-velocity ratios are 0.2–0.26 on every line, so the classification is not discriminating; it is a rubric item.
- **Optical onset not graded.** The archived Re σ_xx is already 83 S/cm at 0.05 eV: the Γ pocket and the Wannier node scatter (2|E| from 64 to 214 meV) smear the onset, so no unique onset is recoverable. The expected 2|E_W| and the smearing are rubric items.
- **Hall conductivity is a lookup with a decision.** σ_xy(0) is read from the right archived folder after interpreting the folder names. A numpy Kubo at 300³ is beyond the lab's 10-minute limit.
- **PBE versus DFT+U.** The archive has no PBE magnetic or DFT+U non-magnetic run, so the ground-state energetics fork cannot be decided from data. It became dp1 (which run is magnetic), and the energy comparison is a rubric item.
- **No surface or Fermi-arc requests.** The archived spectral functions (23–62 MB each) were not used.

## Decision points (9) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | `standard` is non-magnetic (metadata) | model identification |
| dp2 | Wannier band index: 16/17 (hubbard) against 18/19 (standard); `exclude_bands` alone misleads | metadata / cross-check |
| dp3 | The DFT+U Wannier model is exact on its fitting grid but wrong between grid points; DFT lines needed | aliasing exposed only by an off-grid cross-check (consequence visible only after buying the DFT lines) |
| dp4 | Magnetic orbits with antiunitary operations and stabilisers | symmetry |
| dp5 | Chirality sign and k_z displacement | convention + DFT z-line |
| dp6 | Kubo folder names are meshes; non-monotonic convergence | mislabelled runs, found in the input |
| dp7 | Archived σ_yz and σ_zx are numerical artefacts | data-quality judgement |
| dp8 | E_F printed by a line nscf run | metadata trap |
| dp9 | PBE nodes merge on coarse grids; avoided crossings at 6.6 meV | true node versus gap minimum |

dp2 cascades into both censuses and both energies. dp3 cascades into the energy and the separation.

## Calibration

The full table is in `design.json` under `calibration`.

| Field | Reference | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| node energy (meV) | −97.57 | −96.55 to −98.27 (x, y, z lines; fit windows 0.02–0.05 Å⁻¹) | −91.89 (Wannier pair at the DFT line), −73.22 (Wannier family mean), +38.88 (PBE) | ±3 |
| pair separation (Å⁻¹) | 0.3689 | 0.3690 (x and z crossings) | 0.3852 (Wannier pair at the line, +4.4 %), 0.3961 (Wannier family), 0.1546 (PBE) | ±2 % |
| k_z offset, χ = +1 (Å⁻¹) | −0.00979 | −0.00974 (x-line passage) | +0.00979 (upper-band chirality), 0 (symmetric node) | ±0.002 |
| PBE node energy (meV) | 38.88 | 38.38–39.24 (single nodes) | −73.82 (line-run E_F), −361.3 (hubbard band index) | ±4 |
| σ_xy (S/cm) | 114.05 | 114.39 (210³), 110.91 (120³, −2.8 %) | 108.76 (30³, −4.6 %), 102.95 (58³), 74.41 (folder 10) | ±3.5 % |
| σ_yz, σ_zx (S/cm) | 0, 0 | 0 | −0.20 / −4.46 (raw 300³) | ±1.5 |
| node count | 24, 24 | — | 72 (offset 34), 119 (PBE at Wannier 16/17), 30 (avoided crossings), 12 (coarse PBE grid) | exact |
| orbit sizes | [4, 4, 4, 4, 8] | — | [24], [4×6], [8, 8, 8], [8, 16] | exact |

Two values are not ablations:

- **Memorised paper values:** −96 meV and 0.37 Å⁻¹ pass. The paper's k_z shift "0.0047 Å⁻¹" fails by 0.0051, and "all nodes equivalent" ([24]) fails.
- **A single Wannier node near the DFT pair:** −95.7 meV passes the energy, but its separation fails.

## Oracle

The oracle makes 18 evaluations for 23 of 28 units (82 %). In order:

1. Two `wannier_bands` calls against the staged scf eigenvalues give the band offsets.
2. Hubbard census: a 24³ scan, 3 `find_nodes` calls and 3 `berry_flux` calls.
3. Standard census: a 36³ scan and 5 `find_nodes` calls; the 6.6 meV minima are rejected.
4. `hubbard/nscfXYZ` and the Wannier bands on its x-line.
5. Kubo 300 and 210.
6. Orbits from the scf symmetry list.

No literature value is used, so every graded window is reached from purchased data. The slowest request is the 36³ PBE scan, at 39 s.

Results:

- **Gate (2026-10-20, cpuset 6,7,22,23):** passed. All 17 ablations fail, the 7 robustness variants pass, rounding to 4 significant digits passes, and grading is deterministic. The task has 9 decision points and 8 rubric items. The package is 36.3 MB with 10 hidden values, and the oracle took 64 s in the gate container.
- **Harbor oracle (`runs/build-inmnti2-weyl-census/oracle-r1-inmnti2-weyl-census`):** reward 1.0, runtime 1 min 24 s, 23 of 28 units.

The estimated agent length is 50–80 tool calls.

## Open risks

- **Dev-easy risk.**
  - dp1, dp2 and dp8 can be defused by reading headers.
  - The strongest discriminators should be dp3, dp4 and dp6:
    - dp3: grid validation succeeds while the model is wrong;
    - dp4: antiunitary operations and stabilisers;
    - dp6: folder names.
  - An agent that buys `hubbard/nscfXYZ` and compares it with the Wannier model will see the discrepancy at once.
- **Thin margins.**
  - σ_xy: 120³ passes by 0.7 points of margin, and 30³ fails by 1.1 points.
  - The separation of the Wannier pair at the DFT line fails by 2.4 points.
  - The k_z offset has the widest margins but discriminates only the sign and 0. The Wannier family mean (−0.0095) passes.
- **Ambiguity risk on σ_yz and σ_zx.** The instruction asks for the crystal's values. An agent reporting the model's −4.46 S/cm fails; I consider symmetry-zero the unique physical answer.
- **Memorisation.** Paper values pass for the energy and the separation; the counts (12 pairs) also pass.
- **Size.** `standard.npz` is 25 MB (limit 40), and the package is 36.3 MB (limit 60).

## Disk

- **Peak:** about 0.55 GB. That is 160 MB of zip, 280 MB extracted during exploration (205 MB in the final extraction), 60 MB of npz and about 20 MB of scratch (paper PDF, scans), plus the gate's lab image layer.
- **After the build:** raw downloads, scratch files and the gate image were deleted. The refetch script, `members.md5` and `census.json` (80 kB) are kept.
