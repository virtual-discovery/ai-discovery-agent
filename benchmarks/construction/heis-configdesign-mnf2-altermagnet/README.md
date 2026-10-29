# Construction of heis-configdesign-mnf2-altermagnet

Benchmark task (see `benchmarks/tools/BUILDING.md`). Released start, from the same Materials Cloud record as `heis-configdesign-fe` (3p6fv-8hd59), but only the MnF2 archive is used. Nothing is imported from the Fe construction folder; only the general approach (pool extraction, `make_task.py` pattern) was adapted.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) ...`). Run from this folder.

1. `python3 fetch.py [CACHE]` downloads `MnF2.tar.gz` (51,865,693 bytes, md5 4311cefc120956738f6d26440354f3fd) into `~/.cache/ai-theorist/heis-configdesign-mnf2-altermagnet/`. It checks the md5 and unpacks the archive without POTCARs or `__pycache__`.
2. `python3 -B extract.py [CACHE]` writes `pool.json`: 101 spirals with their total energies, CSV energies, VASP warnings, INCARs, k-meshes, moments and SCF step counts. It asserts that the CSV energies are 1000 × ΔE of the 6-atom cell, and that the CSV moments equal the OSZICAR `M_int`.
3. `python3 -B reference.py` writes `j_reference.json`: a 10-shell fit to the 98 spirals with a complete basis.
4. Two Monte Carlo runs give the reference T_N. Each takes about 15–40 min on 4 logical CPUs; `mc_tn.c` is compiled automatically.
   - `python3 -B tn_mc.py j_reference.json tn_coarse.json --sizes 8 12 16 --temps-K 112 114 116 118 120 122 124 --neq 3000 --nmeas 20000`
   - `python3 -B tn_mc.py j_reference.json tn_reference.json --sizes 12 16 20 --temps-K 121.5 122.25 123 123.75 124.5 --neq 5000 --nmeas 40000`
5. `python3 -B make_task.py` writes:
   - the catalog, lab data and private targets;
   - `solution/ablation_values.json`;
   - the calibration numbers, printed to the console.
6. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task heis-configdesign-mnf2-altermagnet --cpuset <set>` writes `check.json`.

`model.py` holds the two-sublattice shell model and the wrong-convention variants used for the ablations.

## What the data showed (design inputs)

- **Normalization.** The archive README documents "meV / atom", but the CSV energies equal 1000 × (E_cell(q) − E_cell(0)) of the 6-atom cell to within 1.5 µeV. They are per cell (2 Mn), and the authors' J values correspond to that convention. This is decision point dp1; the label is quoted verbatim to the agent next to the raw total energy.
- **Uniform two-sublattice spirals.** The INCARs rotate Mn2 in the cell by −2π q·τ(Mn1). Together with VASP's +2π q·R rotation between cells, every Mn has the azimuth 2π q·r. This uniform model fits the 98 clean spirals to 0.083 meV/cell with 8–12 shells. The opposite rotation sense, and a lattice-phase-only model, fit no better than 19–28 meV/cell. q = 0 is the ferromagnet, and integer q with h + k + l odd is the Néel state, the minimum at −95.65 meV/cell.
- **Basis truncation (dp2).** Every OUTCAR prints "choose ENMAX > X". X follows (√ENINI + √(ħ²/2m)·|q|/2)², and with ENCUT = 1100 eV it exceeds ENCUT only for [214] at t = 0.90, 0.95 and 1.00 (X = 1121, 1147 and 1174 eV). Against the clean model these points deviate by +0.04, +0.69 and +1.45 meV/cell. The [214] t = 1 spiral, q = (2,1,4), is the Néel state, yet it lies 1.55 meV above the [001], [111] and [113] Néel energies, which agree with one another to 0.02 meV. A fit that keeps these points mispredicts the hidden Néel-state energies by 0.7–0.8 meV/cell. Over 300 randomized budgeted plans the effect is robust: the 5th percentile is 0.68 with the truncated points and the 95th percentile is 0.14 without them. The hidden-direction RMSE alone does not discriminate, because plan-to-plan variance (0.16–0.45) is larger than the contamination effect (+0.1). That is why the ground-state check exists.
- **Aliasing.** Rank of the 10-shell design: [001] alone has rank 4 (J2 aliases with the 8.0 Å shell, J2 = −0.60), [121] alone has rank 4, and [001] + [121] has rank 8. [214] or [113] alone is full rank. With [113] and [111] hidden, [214] is the only available direction that removes the aliasing, so the trap in dp2 cannot be avoided by skipping [214] altogether (dp3). Commensurate grids (t = 0.2, 0.4, …) alias long-range shells (dp4).
- **Model order (dp5).** The clean in-sample RMSE is 0.53 (3 shells), 0.40 (6), 0.083 (8) and 0.083 (10) meV/cell. The shells that matter beyond J4 are the 7.0 Å and 7.8 Å same-sublattice shells: J7 = −0.060 and J8 = +0.042 meV. These are exactly the altermagnetic-split shells.
- **Altermagnetism (dp8).** For the Mn1 site symmetry mmm with axes [110], [1-10] and [001], the first same-sublattice shell that splits is a⟨110⟩ (7.007 Å). Every inter-sublattice shell, and J1, J3, J4 and J6, are single classes. A uniform spiral gives both sublattices the same q, so the design columns for J(Mn1, a[110]) and J(Mn2, a[110]) are identical, and only their sum is identifiable.

## Build notes and deviations from the design

- **Hidden set.** All [111] and [113] spirals (40) are hidden. Γ, [001], [121] and [214] (61) are purchasable, with a budget of 24 (39 % of the catalog). The design suggested hiding [121]/[214], but [214] holds the basis-truncated points and must stay purchasable for dp2 to exist. [113] alone nearly determines the model, so hiding it makes the direction choice binding.
- **Energy units.** The design assumed the energies were per atom (6 atoms per cell). They are per cell. We grade per-cell predictions and per-Mn J values, and the mislabel became dp1. We present it as documented, not altered.
- **Grading tolerances.** Plan-to-plan variance of the hidden RMSE is large, so the design's RMSE < 0.3 meV per "atom" was replaced by two checks. RMSE ≤ 0.40 meV/cell catches model order, aliasing, normalization and phase errors. The ground-state energy check (mean error at the two hidden Néel-state spirals ≤ 0.35 meV/cell) catches basis-truncated data.
- **J tolerances.** The design's "J1–J3 within ±0.15 meV" became windows. J3 ≈ −0.02 meV, so ±0.15 meV would accept zero; the window is [−0.09, 0.03].
- **Altermagnetism is graded numerically.** The design had it as a rubric point only. Three checks apply: the first split-shell distance (7.007 ± 0.05 Å), its identifiable sublattice mean (window [−0.11, −0.025] meV), and the difference, which must be `null`.
- **T_N.** The reference is the classical Monte Carlo value for the reference model (unit vectors), with tolerance ±12 % instead of the design's 15 %; the nearest wrong path, mean field, is +34 %. The design's decision point on quantum versus classical spin length was dropped, because the prompt must fix the classical convention to make T_N unique, so that fork is not genuine. Experimental T_N (67 K) and neutron J values are never used.
- **Decision points by kind.** dp1 and dp2 are discovered from raw outputs; dp3 is a budgeted choice of discriminating purchases; dp4 is aliasing exposed only by an off-grid cross-check; dp5 and dp8 are regime judgements; dp6 and dp7 are conventions and estimators. Cascades: dp1 → J and T_N; over-reacting to dp2 by dropping all of [214] → dp3 aliasing.

## Oracle

The oracle buys Γ, [001] (t = 0.15, 0.35, 0.55, 0.75, 0.9, 1.0), [121] (t = 0.1–0.5, using the E(t) = E(1 − t) symmetry) and [214] (t = 0.05–0.85). That is 21 of 24 units.

- Before buying [214], it checks the documented normalization against total energies and calibrates the basis criterion on VASP's printed X values.
- It chooses 8 shells by leave-one-out CV.
- It runs a numpy heat-bath plus over-relaxation Monte Carlo with a Binder crossing for L = 8/12 (about 2–5 min on 2 CPUs; 123.2 K against the 122.9 K reference).
- It assigns the split shell by a site-symmetry orbit test.

Calibration table: `tasks/.../design.json` → `calibration`. Summary (reference → oracle | robustness | nearest wrong path → tolerance):

| field | reference | oracle | robustness | nearest wrong path | tolerance |
|---|---|---|---|---|---|
| hidden RMSE (meV/cell) | – | 0.165 | 0.176 (12 shells), 0.176 (all clean available) | 0.253 truncated kept (fails ground state), 0.76 (6 shells), 2.1 (commensurate grid) | ≤ 0.40 |
| ground-state error (meV/cell) | – | +0.108 | +0.063, +0.133; clean random plans 95th pct 0.135 | +0.689 / +0.829 truncated kept; contaminated random plans 5th pct 0.68 | ≤ 0.35 |
| J1 (meV) | −1.126 | −1.134 | −1.138, −1.126 | −0.61 (grid), −2.27 (per Mn) | [−1.25, −1.00] |
| J2 (meV) | −2.983 | −2.983 | −2.977, −2.982 | +2.98 (Néel reference), −0.60 ([001] only) | [−3.10, −2.86] |
| J3 (meV) | −0.027 | −0.013 | −0.017, −0.016; random plans 1–99 % [−0.026, 0.005] | −0.106 (grid), −0.372 (no [214]) | [−0.09, 0.03] |
| T_N (K) | 122.9 (C MC, L = 12/16/20) | 123.2 (numpy MC, L = 8/12) | 123.0 (coarse L = 8/12) | 164.6 (mean field) | ±12 % |
| split shell | 7.007 Å, mean −0.060 meV | 7.007 Å, −0.071 | −0.067, −0.070; random plans [−0.095, −0.064] | 4.955 Å; mean 0 (≤ 6 shells); difference 0.0 or 0.032 | ±0.05 Å; [−0.11, −0.025]; null |

The reference T_N runs gave Binder crossings at 123.1 K (L = 12/16) and 122.7 K (16/20), and at 123.0 K (8/12) and 123.6 K (12/16) in the coarse scan. T_N/T_MF = 0.75.

## Results

- **Gate (2026-10-19 18:02 UTC, after the raw downloads were deleted):** passed. 8 decision points; all 15 ablations fail; the 4 robustness variations and 4-significant-digit rounding pass. The oracle used 21/24 units in 184 s of container time.
- **Harbor oracle:** reward 1.0 (`runs/build-heis-configdesign-mnf2-altermagnet/oracle-r1-*/*/verifier/reward.txt`).
- **Disk:** peak use was about 0.38 GB (`MnF2.tar.gz` 50 MB plus 324 MB unpacked without POTCARs). The raw downloads were deleted after the gate passed; `fetch.py` refetches them, with the md5 pinned.

## Open risks

- Agents that sample q poorly and pick 11–12 shells by leave-one-out CV can exceed the 0.40 RMSE tolerance even without other errors. Of randomized reasonable plans, about 10 % exceed it. This is the intended model-order and aliasing judgement, but it adds variance.
- The `null` check for the altermagnetic difference is binary. The surrounding checks (shell distance, mean) and the rubric carry the reasoning.
- An agent with an accurate but less converged Monte Carlo (for example, a susceptibility peak in one small cell) may land 5–10 % off. The tolerance is 12 %.
- MnF2 altermagnetism and the first split shell (J7a/J7b at a√2) appear in recent literature, so the shell distance may be recalled. The numbers are protocol-specific (PBE, no U).
