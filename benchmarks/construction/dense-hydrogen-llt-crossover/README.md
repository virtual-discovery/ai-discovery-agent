# Construction of dense-hydrogen-llt-crossover

This is a released-start the benchmark task (see `benchmarks/tools/BUILDING.md`). It is built from Tenti, Jäckl, Nakano, Rupp and Casula, *Hydrogen liquid-liquid transition from first principles and machine learning*, Phys. Rev. B (2025), doi:10.1103/pbrk-3zgd. The data are in Zenodo record 14790037 (CC BY 4.0).

## Rebuild

All compute steps run in `ai-theorist-benchmark-runtime:20261007` with `--cpuset-cpus <yours> --user $(id -u):$(id -g) -e HCACHE=/cache -v ~/.cache/ai-theorist/dense-hydrogen-llt-crossover:/cache`. The fetch steps need network and run on the host.

1. `python3 -B fetch_thermo.py`: lists the 23 classical tars by walking their tar headers with HTTP range requests (`remotetar.py`). It then reads every `input.lammps`, `log.lammps` and `slurm-*.out` member (1491 members, about 1 GB) and checks each against `members_thermo.json`. Zenodo rate-limits (HTTP 429); the client honours Retry-After.
2. `python3 -B fetch_pimd.py`: reads all 175 PIMD tables from the nested tars inside `PIMD.tar` (about 0.8 GB), checked against `members_pimd.json`.
3. `python3 -B fetch_frames.py`: reads six byte ranges of every nested `traj/trajectory.tar` at 900–1100 K (221 runs, about 0.6 GB), checked against `members_frames.json`.
4. `extract.py` (container) writes `pool_thermo.npz`, `pool_pimd.npz` and `pool_meta.json`.
5. `structure.py 0 4` … `structure.py 3 4` (four containers, one CPU each, about 25 min), then `structure.py merge`. This writes `pool_structure.json`.
6. `make_task.py` (container, `-e TASK=/task` with the task folder mounted) writes:
   - `environment/lab/data/*`;
   - the catalog;
   - both `public.json` files;
   - the private targets;
   - `hidden_values.json`;
   - `solution/ablation_values.json`.

   It also asserts the evidence behind each exact-graded answer (loops by N and T, S(k)/N at N = 2048, and the PIMD label duplication).
7. `python3 -B write_design.py` writes `design.json`, including the calibration table and the oracle values in `oracle_values.json`.
8. `.venv/bin/python -m benchmarks.tools.check_task dense-hydrogen-llt-crossover --cpuset <yours>`.

Other files:

- `analysis.py`: the reference protocol;
- `lmplog.py`: the LAMMPS log parser;
- `traj.py`: frame parsing and S(k).

## Archive files (Zenodo md5)

| file | bytes | md5 |
|---|---|---|
| N_128_T1000K.tar | 2456360960 | 50c4fde34416b06943b16e31f48ae9ed |
| N_128_T1100K.tar | 2448209920 | ee60cca9d81d03e491669bbb809b8493 |
| N_128_T1200K.tar | 868577280 | 6d6155e0edc84e06a0191fb069a5318d |
| N_128_T1300K.tar | 779008000 | 483b651a2c653339d7c5c3c1d619f22e |
| N_128_T1400K.tar | 776181760 | 788cf99436d90fb9811a4066abb3e9f0 |
| N_128_T900K.tar | 2195732480 | 6d3f79f06c50b1d44ff184f39d743ab9 |
| N_128_T950K.tar | 2612111360 | 1eed1b0681547fbccfbf287210776243 |
| N_2048_T1000K.tar | 7222497280 | f6290882f73409aabc0b88503f1034ae |
| N_2048_T1100K.tar | 4405657600 | 161bd2913dafac130e58280ecc386fcd |
| N_2048_T1200K.tar | 3425914880 | c61fa5feb52d258b19423f85a4bb7a4f |
| N_2048_T1400K.tar | 2930073600 | 09bf9bfa21965efd338796290c7d655e |
| N_2048_T1600K.tar | 3536885760 | dd3255eda3d98204ce6b4b0811b1d4dd |
| N_2048_T1800K.tar | 4511109120 | 968348bd5a73b3981188557a9ea15f61 |
| N_2048_T900K.tar | 5513656320 | 2ac2d3f56319a7b7f715ffffbf7e2441 |
| N_2048_T950K.tar | 6854748160 | 321e6d18e789ee3d6b09912afbd2657c |
| N_256_T1000K.tar | 3052728320 | 450604bbd1975a014549c4ea9c547175 |
| N_256_T1100K.tar | 2665809920 | a5c42dec35cf4760d43c0e8bdaf8dcf0 |
| N_256_T900K.tar | 4115251200 | cdea7ae5bc4bd930b04b97a94eca3e8e |
| N_256_T950K.tar | 2860984320 | 2bdfff54c5371ac455b5aee933e78ac5 |
| N_512_T1000K.tar | 3741696000 | e4da77a4f34a31c2ccc0b152e91f915b |
| N_512_T1100K.tar | 3736576000 | e7bad14c7987fb33b9ee673c530a23fc |
| N_512_T900K.tar | 4502804480 | 066185596491aae1c64485be2bb4f159 |
| N_512_T950K.tar | 4030689280 | 74b692d10016c57e2043f3dd8b14a953 |
| PIMD.tar | 831272960 | 2dc9ab6ea19227697d7822b8739c09ad |

## Findings in the archive that shaped the design

- **Restart chains.** Long runs are chains of 20–200 ps restart jobs, and `log.lammps` holds only the last job. The full series comes from stitching all `slurm-*.out` by step. Some jobs failed, for example a missing restart file, and these are listed in the served job history. Runs start from `init.lmp` (`read_data`).
- **Start-up transients.** The first few ps after `read_data` carry a drift that dominates var(E) for the short runs. The N = 2048 runs at 1200–1800 K last only 20–40 ps, and 900–1000 K runs started from a crystal melt within about 10–30 ps (S/N 0.1 → 0.01). Whole-run c_v is 2–3 times too large and moves the Widom line to lower rs by 0.010–0.020. This is dp1.
- **Loops.** The largest dP/drs (GPa per unit rs) along each isotherm:

  | T (K) | N = 128 | 256 | 512 | 2048 |
  |---|---|---|---|---|
  | 900 | +153 | +929 | +1045 | +1321 |
  | 950 | +219 | +650 | −152 | −179 |
  | 1000 | +143 | +190 | −253 | −241 |
  | 1100 | −219 | −221 | −239 | −237 |

  Small cells show loops up to 1000 K; at N = 2048 only 900 K does (dp2).
- **Structure.** max_k S(k)/N at N = 2048 is about 0.13 on the 900 K molecular side and about 0.008 in the liquids. N = 128 liquids give 0.05–0.3, which is dp4. The 900 K crystal still diffuses (MSD about 80 Å² over 100 ps), which is dp3.
- **PIMD.** All 55 files under `PIMD/900K` are byte-identical to files under `PIMD/1000K`, and their temperature column reads 1000 K, so no 900 K path-integral data exist (dp6). At 1000 K the PIMD EOS is smooth, with its flattest point near rs 1.488 and 138.8 GPa. The classical N = 2048 crossover is at rs 1.432 and 180.1 GPa.

## Deviations from the design (DETAILS.md section 47)

- **Widom line from the c_v maximum.** The paper's estimator is used instead of a compressibility maximum. It is graded at 1100, 1200 and 1400 K (±0.005 rs) rather than at "two hidden T". Nothing is withheld from purchase: the graded values are derived from the served data and are never served, as in the fluorite task.
- **Peak scaling.** "Saturating vs proportional" is not graded separately. The four per-temperature classifications carry the finite-size result. Peak heights at 950 K for N = 2048 are contaminated by melting transients, so a slope fit is not uniquely recoverable.
- **NQE** is graded as the shift of the crossover pressure at 1000 K (−41.3 GPa, ±8). At 900 K the answer must be `null`, because the 900 K PIMD data are duplicates. The design's "±50 % of the shift" is replaced by ±8 GPa, about 20 %: the defensible spread is −40.0 to −41.9 and the nearest wrong path is −9.5.
- **AIMD not used.** Each `aimd_*K.tar` is about 3 GB of QE output at a single density (rs 1.464, 2 ps). It cannot show an EOS kink, and streaming 13 GB was not justified. The "why short simulations suggest first order" point is a rubric item instead (rubric 4), supported by the served time series.
- **Molecular-fraction cutoff fork dropped.** The structure response serves the paper's 1.05 Å/80 fs definition. Its value does not enter the grade.

## Calibration

The full table is in `design.json` under `calibration`.

| Graded field | Reference | Defensible alternatives | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| Widom rs 1100 K | 1.4484 | 1.4475–1.4486 (t0 5/20 ps, 3-pt, PotEng c_v, grid max) | whole runs 1.4336; N = 128 1.4538 | ±0.005 |
| Widom rs 1200 K | 1.4644 | 1.4641–1.4653 | whole runs 1.4546 | ±0.005 |
| Widom rs 1400 K | 1.5019 | 1.5011–1.5050 | whole runs 1.4820 | ±0.005 |
| PIMD crossover rs 1000 K | 1.4880 | tanh fit 1.4840; no transient removal 1.4871 | classical 1.4319 | ±0.007 |
| NQE shift 1000 K (GPa) | −41.27 | −40.00 (tanh), −41.92 (N = 512 classical), −40.69 | fixed density −9.49; ignored 0 | ±8 |
| NQE shift 900 K | null | — | mislabelled files −55.0 | must be null |
| Classification 900/950/1000/1100 | FO/X/X/X | — | small cells FO/FO/FO/X; absolute S(k) FO/X/X/FO | exact |
| 900 K low-density phase | crystal | — | liquid (diffuses) | exact |

The tightest margins:

- **Widom 1400 K.** The grid maximum (1.505) is 0.0031 off, inside the tolerance with 0.0019 to spare.
- **Widom 1100 K.** The N = 128 isotherm (1.4538) fails by only 0.0004. It is not counted as an ablation on its own; the dp7 ablation also carries the small-cell classification.
- **PIMD rs.** The tanh fit passes with 0.003 to spare.

## Oracle

- **Purchases:** 66 evaluations for 117 of 135 units (87 %); the catalogue is worth 1134 units.
  - N = 2048 c_v isotherms at 1100, 1200 and 1400 K (8, 8 and 7 densities);
  - N = 2048 isotherms near the transition at 900, 950 and 1000 K (6, 5 and 7 densities);
  - N = 128 and 256 at 1000 K (4 densities each);
  - four structure purchases;
  - 11 PIMD densities at 1000 K and two 900K-labelled ones.
- **Gate:** passed (2026-10-26). All 8 ablations fail, all 8 robustness variants pass, and rounding to 4 significant digits passes. The gate container ran for 69 s.
- **Harbor oracle:** reward 1.0 (`runs/build-dense-hydrogen-llt-crossover`).
- **Recall:** the oracle's density windows were chosen with knowledge of the answers. A solver without recall can find them from the catalogue: the authors ran the long runs and the dense rs grids around the transition. Cheap N = 128 scans also locate the peaks within about 0.005 rs. No literature numbers enter any computation.

## Open risks

- **Memorization.** The paper's headline ("supercritical above melting, first order at 900 K in the crystal region") is public. A model that recalls it gets the classification and phase items, but not the Widom rs values, the PIMD shift or the 900 K-null trap.
- **Classification by loops.** The classification uses loops at N = 2048, and the 900 K loop is a crystal–liquid one. That is the intended reading: the question asks about the molecular-to-atomic transition, and the phase item separates the two.
- **Dev-easy risk.** Transient removal (dp1) and comparing crossovers (dp5) are standard. The difficulty rests on discovering the label duplication (dp6), the small-cell loops (dp2) and the absolute-S(k) confound (dp4) under the budget.
- **Peak disk:** about 2.2 GB under the cache (raw text members, PIMD tables and frame chunks). These were deleted after the build. The derived pools (14 MB) and all checksums remain.
