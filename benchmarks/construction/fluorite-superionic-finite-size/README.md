# Construction of fluorite-superionic-finite-size

This is a benchmark task with a released start (see `benchmarks/tools/BUILDING.md`). The source is F. Grasselli, J. Chem. Phys. 156, 134705 (2022), with data in Materials Cloud record yzhmw-vje16 (CC BY 4.0).

## Rebuild

1. `python3 -B fetch.py` reads 1207 members out of the remote zips by HTTP range requests, without downloading whole archives. It verifies each member's CRC32 and its SHA256 against `members.json`. Output goes to `~/.cache/ai-theorist/fluorite-superionic-finite-size/raw/`, about 640 MB, with each `nvt_log` stored as a `.npy` array of its thermo columns.
2. `python3 -B extract.py` writes `pool.npz` and `pool_meta.json` (286 runs) to the same cache.
3. `python3 -B make_task.py` writes:
   - the catalog;
   - `data.json` (the lab payloads);
   - the private targets;
   - `hidden_values.json`;
   - `solution/ablation_values.json`.

   It also prints the calibration numbers.
4. `.venv/bin/python -m benchmarks.tools.check_task fluorite-superionic-finite-size --cpuset 4,5,20,21` writes `check.json`.

Supporting files:

- `analysis.py`: the reference protocol (D, c_V, T_c, Arrhenius fits).
- `remotezip.py`: the range-request zip reader.

## Archive files used (shared source group)

The sibling task `superionic-heat-transport-size-onsager` uses the same record. This task reads only the NVT diffusion folders:

- `CaF2.zip`: `CaF2/NVT_L_T/`
- `PbF2.zip`: `PbF2/NVT_L_T_6p056/` and `PbF2/NVT_L_T_6p056_DW/`
- `UO2.zip`: `UO2/NVT_L_T/`
- `alpha-AgI.zip`: `alpha-AgI/NVT_L_T_5.37ang/`

From each run it uses `msd_CM_B4.dat`, `msd_self_B4.dat`, `nvt_log`, the input deck, and the folder scripts.

It uses none of the `*_heat` folders and not `UO2_heat.zip`. Its targets are diffusion coefficients, c_V peaks and activation energies. It does not use heat flux or thermal conductivity, so its targets do not overlap with the heat-transport design.

## Design as built

- **Lab request:** `md_run` (material, supercell n, temperature).
  - Cost is max(1, ceil(atoms/300)): 1 for ≤ 300 atoms, 2 for 324/500, 3 for 768, 5 for 1500.
  - The budget is 140. The catalog total is 503.
- **Each purchase returns:**
  - the verbatim input deck;
  - 140 lag rows of `msd_CM_B4.dat` (220 for UO2);
  - a thermodynamic summary of the full log (means, 20-block standard errors, total-energy variance).
- **Withheld runs:** PbF2 6×6×6 and 7×7×7 (all T), AgI 6×6×6 (all T), CaF2 4×4×4 1250 K, and UO2 4×4×4 2450 and 3150 K.
- **Graded (15 checks):**
  - T_c of the 768-atom cells (±120 K);
  - CaF2 and PbF2 ordered and superionic activation energies (±0.25, ±0.25, ±0.15 and ±0.10 eV);
  - D of seven withheld runs (±0.15 dex each).

## Deviations from the design (and why)

- **Thermodynamic-limit T_c is not graded.** Only three sizes exist for CaF2 and UO2, and 1/L versus 1/N extrapolations of the c_V peak differ by more than 100 K, so the limit is not recoverable. The graded T_c is that of the 768-atom cell. The report must say that the limit is not supported.
- **The "which material shows weak finite-size effects" bit is not graded.** It is one bit and stated in the abstract, so it is a rubric item. AgI enters the grade instead through two 864-atom predictions, with a "fluorite correction transferred to AgI" ablation.
- **UO2 activation energies are not graded.** ln D(1/T) of UO2 is curved on both sides of T_c, and defensible windows give 3.4–5.1 eV (ordered) and 1.0–2.5 eV (superionic). The rubric asks for this to be stated.
- **The PbF2 superionic E_a is graded but does not separate cell sizes.** The 96-atom cell gives 0.29–0.48 eV against 0.26–0.34 eV for 768 atoms. It is kept as a robust measured quantity with a ±0.10 eV tolerance.
- **The block-averaging decision point was replaced.** The lab serves block standard errors, and the data contain no replicas: the DW set duplicates the main PbF2 runs byte for byte. Separating noise from signal is therefore a rubric item. Difficulty comes instead from evidence in the data and metadata:
  - **dp1:** the UO2 decks use a 2 fs timestep (0.02 ps frames, 8000-row MSD tables); the other decks use 4 fs.
  - **dp3:** the 768-atom PbF2 cell agrees with 1500/2592/4116 atoms to within 7% at 800–1200 K but is low by a factor of 1.9 at 600 K. The discriminating purchase is the 1500-atom run below T_c.
  - **dp6:** the 324-atom cell gives T_c within 30 K of 768 atoms for CaF2 and PbF2, but 200 K higher for UO2.
  - **dp2 and dp4:** the size effect changes sign with temperature. Small cells overestimate D in the ordered regime and underestimate it near T_c, which defeats a Yeh–Hummer 1/L correction.
- **NVE versus NVT is not offered.** Diffusion NVE runs exist only for PbF2 n = 2, 3, and the author shows they agree with NVT.

## Calibration

Values are from `make_task.py`. The full table is in `design.json` under `calibration`.

| Graded field | Reference | Defensible alternatives | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| T_c CaF2 (K) | 1495.6 | fluctuation c_V 1403.9 | 96-atom 1641.5 | ±120 |
| T_c PbF2 (K) | 801.8 | fluctuation c_V 725.7 | 96-atom 943.0 | ±120 |
| T_c UO2 (K) | 2742.9 | fluctuation c_V 2795.9 | 324-atom 2939.5; 96-atom 3314.1 | ±120 |
| E_a CaF2 ordered (eV) | 1.880 | oracle subset 1.815; 10–80 ps 1.828; windows 1.83–2.11 | single Arrhenius 1.422; 96-atom 1.232 | ±0.25 |
| E_a CaF2 superionic (eV) | 0.682 | 0.679 / 0.696 / 0.742; windows 0.64–0.83 | 96-atom 1.099; single 1.422 | ±0.25 |
| E_a PbF2 ordered (eV) | 0.619 | 0.620 / 0.613 / 0.528; windows 0.61–0.68 | 96-atom 0.220 | ±0.15 |
| E_a PbF2 superionic (eV) | 0.304 | 0.290 / 0.288 / 0.330; windows 0.26–0.34 | 96-atom 0.420; single 0.485 | ±0.10 |
| D predictions (dex) | released runs | ≤ 0.086 (oracle, 10–80 ps, Arrhenius fit, 1/N, momentum frame) | −0.19 to −0.32 (96-atom, UO2 frame time, 768 assumed converged, 324-atom UO2), +0.23 (AgI corrected), negative (Yeh–Hummer) | ±0.15 each |

The tightest margins:

- **T_c, wrong side:** 96-atom CaF2 and PbF2 fail by 26 K and 21 K.
- **T_c, defensible side:** fluctuation c_V for CaF2 passes with 28 K to spare.
- **PbF2 superionic E_a:** the 96-atom value fails by 0.016 eV, but other checks also catch that path.
- **PbF2 600 K prediction:** a linear 1/L extrapolation from the 768- and 1500-atom cells gives +0.15 dex, at the edge of the tolerance. This path is not counted as an ablation. Using the 1500-atom value (−0.03) or a 1/N extrapolation (+0.09) passes.

## Oracle

- **Purchases:** 43 runs for 121 of 140 units:
  - 768-atom scans for CaF2 (11), PbF2 (12) and UO2 (10);
  - PbF2 1500 atoms at 600 and 800 K;
  - AgI at 32 and 500 atoms at 600 and 900 K;
  - four small-cell comparisons.
- **Gate (container):** 8 s.
- **Gate (2026-10-19):** passed. All 10 ablations fail and all 5 robustness variants pass.
- **Harbor oracle:** reward 1.0 (`runs/build-fluorite-superionic-finite-size`), rerun on the final files after a cosmetic note change. The first attempt was blocked by another builder's `__pycache__` (in `fresh-inp-polar-phonons`) until the runner was changed to verify only the selected task.

## Open risks

- **Memorization.** The paper quotes T_c ≈ 1450/800/2800 K for 768 atoms. A model that recalls them falls inside the T_c tolerance, but the activation energies and the seven predictions still require the data.
- **Broad peaks.** The CaF2 and PbF2 c_V peaks of the 768-atom cells are broad (within noise over about 100 K), hence the ±120 K tolerance. The UO2 peak is sharp.
- **Fixed lattice parameter.** These are empirical force fields at a fixed lattice parameter (NVT). The graded values are protocol-specific and not experimental.
- **Peak disk use** was about 1.0 GB under the cache: raw members plus `pool.npz`. These were deleted after the build; the refetch scripts and checksums remain.
