# Construction of ti-codoping-exchange-tc

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). The source is Mozumder, Wasmer, Antognini Silva, Blügel and Rüßmann, *High-throughput magnetic co-doping and design of exchange interactions in topological insulators*, Materials Cloud record 6txnj-70460 (DOI 10.24435/materialscloud:b7-6k, CC BY 4.0), arXiv:2407.04413. Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

From this folder. All compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) -e HOME=/tmp -e TI_CACHE=/cache -v <cache>:/cache -v <repo>/benchmarks:/b -w /b/construction/ti-codoping-exchange-tc ...`).

1. `python3 -B fetch.py [CACHE]` downloads six small record files (1.9 MB) through the per-file API links and checks md5. The default CACHE is `~/.cache/ai-theorist/ti-codoping-exchange-tc`. The AiiDA exports (0.4–46 GB) are not needed.
2. Compile the lab kernel into the cache: `gcc -O2 -o /cache/mc_dilute /b/tasks/ti-codoping-exchange-tc/environment/lab/mc_dilute.c -lm`.
3. Run the reference Monte Carlo with the lab's own `dilute.py` and `mc_dilute.c`, on larger cells, with more placements and longer runs than the service allows. All runs together take about 4 h on 4 logical CPUs. The runs are listed in `REFERENCE_RUNS` below.
   - `python3 -B reference.py run TAG --comp ... --sizes ... --samples ... --temps ...`
   - `python3 -B reference.py summary TAG` prints the Binder ratios, the scaled magnetizations and their crossings with jackknife errors.
4. `python3 -B calibrate.py` writes:
   - `reference.json` (targets);
   - `ablations.json` (wrong-path and robustness values);
   - `calibration.json`;
   - the `calibration` block of the task's `design.json`.
5. `python3 -B make_task.py` writes:
   - the catalog;
   - the lab `data.json` (all 2077 archived rows, grouped by pair);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - both `public.json` budgets.

   It asserts the two archive checks described below.
6. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task ti-codoping-exchange-tc --cpuset <set>` writes `check.json`.

Modules:

- `model.py`: table loader, shells, single-impurity moments, the mean-field formulas.
- `reference.py`: reference runs and crossing analysis.
- `calibrate.py`: targets, ablations, robustness and the calibration table.
- `make_task.py`: staging.
- `environment/lab/dilute.py`: lattice and random placement. It is shared with the lab, so the reference uses the agent's code path.
- `environment/lab/mc_dilute.c`: heat-bath plus over-relaxation kernel.

`solution/solve.py` is a stand-alone client of the lab.

REFERENCE_RUNS: each line is tag, composition, sizes with placements, temperatures in K. All runs use `--neq 1500 --nmeas 5000` unless noted.

- V5Cr5 (`V:0.05,Cr:0.05`): sizes 3 4 5 6, placements 48 32 24 12; T = 18 16 15 14 13 12 11 9
- Cr5Mo5 (`Cr:0.05,Mo:0.05`): sizes 3 4 5 6, placements 48 32 24 12; T = 22 20 19 18 17 16 15 13
- Cr5Mn5 (`Cr:0.05,Mn:0.05`): sizes 3 4 5 6, placements 48 32 24 12; T = 13 11.5 10.5 10 9.5 9 8 6.5
- Cr5Nb5 (`Cr:0.05,Nb:0.05`): sizes 3 4 5, placements 48 32 24; T = 6 4 3 2.5 2 1.5 1.2 1
- Cr5Nb5_slow (same composition, independent placements): sizes 3 4 5, placements 32 24 16; T = 14 11 9 7.5 6.5 5.5 4.5 3.5 2.5 2; `--neq 3000 --nmeas 6000`
- V2Cr2 (`V:0.02,Cr:0.02`): sizes 3 4 5, placements 48 32 24; T = 4 3 2.5 2 1.5 1.2 1 0.8
- V5Cr5_coarse, Cr5Mo5_coarse, Cr5Mn5_coarse: size 3, 32 placements; T = 40 34 29 25 22 19 17 15 13 11 9 7 5 3 (the single-cell susceptibility peak, recorded but not an ablation)
- V5Cr5_near, Cr5Mo5_near (`--shells 0-6`): sizes 3 4 5, placements 32 24 16; T = 16 14 12 11 10 9 8 6 and 20 18 16 15 14 13 12 10
- Cr5 (`Cr:0.05`): sizes 3 4 5, placements 32 24 16; T = 8 6 5 4 3.5 3 2.5 2
- Cr5Co5_zero (`Cr:0.05,Co:0.05 --coco-zero`): sizes 3 4 5, placements 32 24 16; T = 16 13 11 9 7 5 3 2

`oracle_run.json` holds the gate oracle's submission; `calibrate.py` copies it into the calibration table.

## What the data showed (design inputs)

- **Convention verified.** The archived mean-field matrix `Tcs_Bi2Te3.txt` equals (1/3k_B) Σ_shells N J, with the authors' co-doping rule (T_AA + T_BB + T_AB)/3, to 1e-4 K per unit concentration. This confirms the paper's pair-once Hamiltonian, H = −½ Σ_{i≠j} J_ij ê_i·ê_j with unit vectors. The matrix's non-magnetic mask is permuted by a notebook bug: Ru (0.225 μB) is blanked and Zr (3e-8 μB) is kept. The matrix is not staged.
- **The 17 shells.** They match the Bi sublattice of the CIF: each distance belongs to one sublattice pair (Bi1–Bi1 or Bi1–Bi2), and the counts equal `num_neighbors`. Assigning couplings by distance on the real two-layer lattice is therefore unique.
- **`m2` column defect.** It is asserted in `make_task.py`. The notebook reads `magmoms[[0, 20]]`, and index 20 is impurity 2 only in the shell-0 cluster. For shells 1–16, `m2` ≈ 0.01 μB for every pair, including Mn (4.4 μB).
- **Missing calculations.** Co–Co is absent in shells 0–6 and 16. Its shell-9 run collapsed to a non-magnetic state (m1 = 0.005 μB). Scattered gaps elsewhere include Mn–Fe 6/16, Ti–Ti 2, Fe–Fe 8, Tc–Tc 9, Cr–Ni 11 and several Ru pairs. The graded ordered compositions use only complete pair tables.
- **Percolation-limited order.** At 2% + 2% the couplings strong enough to matter above 1 K do not percolate. The graph percolates only through |J| ≳ 0.1 meV, and nothing orders down to 0.8 K, although mean field gives 33 K. At 5% + 5% order appears at 11–20 K, three to eight times below mean field (60–97 K).
- **The far shells are essential.** With shells 0–6 only, V5Cr5 and Cr5Mo5 do not order at all in the probed range (down to 6 and 10 K). The far shells raise the mean degree from about 3 to about 10.
- **Frustration decides the co-dopant.** Nb has the largest mean-field T_C with Cr (97 K) because of its large Cr–Nb and Nb–Nb couplings. Its antiferromagnetic Nb–Nb shells at 6.10 Å (−2.07 meV) and 7.51 Å (−2.43 meV, 6 neighbours) frustrate the network. Cr5Nb5's Binder ratio stays at about 1.2 down to 1 K. Two independent reference sets give either no crossing or a weak, drifting one (10.6 → 7.7 K), all far below Cr5Mo5 (19.8 K). Mo's large antiferromagnetic shell (−22 meV at 6.49 Å) has only 3 neighbours, and its ferromagnetic Cr–Mo couplings carry the order.

## Deviations from the design (candidates.json / DETAILS.md §18)

- **Concentrations.** At the suggested 2% + 2% nothing orders above 0.8 K for V/Cr (the paper's showcase pair; mean field 33 K), so a T_C there is not gradable as a number. The graded ordered compositions use 5% + 5% (10% of Bi sites, within the experimental 10–20% range quoted by the paper). V2Cr2 is kept as the percolation case (graded ≤ 3 K, 0 expected).
- **No ground-state-type (FM/ferri/none) categorical.** Among pairs with complete exchange tables none is ferrimagnetic; the only clean ferrimagnet candidates (Nb–Tc, Cr–Fe) have missing archived dimers. Cross-species frustration enters instead through the co-dopant question (Nb vs Mo) and Cr5Mn5 (Mn–Mn antiferromagnetic, Cr–Mn ferromagnetic).
- **"Best pair ranking" became `best_codopant_for_Cr`** over {V, Nb, Mo, Mn} at 5% + 5%: Mo, with a margin (19.8 K vs ≤ 12.0 K upper bound for Nb, 10.8 K for Mn) beyond the 20% tolerance. Cr5Nb5 itself is not graded as a number because its weak, frustrated order is not uniquely recoverable at affordable sizes.
- **Non-magnetic dopants.** The design's "dopants with vanishing moments" fork is realised by the archive's own moment criterion (stated in the prompt for uniqueness) and by the Cr5Co5 missing-data case; Rh/Pd induced-moment couplings (Cr–Rh 4.8 meV at 6.49 Å) are real in the archive but a Cr+Rh composition reduces to Cr alone and was not added to keep the agent's MC load bounded.
- **Lab Monte Carlo service.** The design allowed "a lab MC service at higher cost". The lab runs the agent's own exchange model (so every modelling choice cascades), on sizes L = 2–6 with 1000 + 3000 sweeps; the budget binds the number of placements and sizes. Agents may also run their own MC.
- **Substitution-site fork dropped.** Bi1 and Bi2 are inversion-equivalent and the table lists every Bi1–Bi1/Bi1–Bi2 shell once; there is no separate Bi2 table to mix.
- **Tolerance.** 20% as designed. Defensible estimators (Binder or scaled-magnetization crossings, any size pair L = 3–6) span −12% to +6%.

## Decision points (7)

| dp | Fork | Ablation (graded result) |
|---|---|---|
| dp1 | mean field (two-sublattice or authors' average) | 60–82 K / 22–39 K vs 11–20 K; Nb best |
| dp2 | one random placement per size | V5Cr5 14.9 K (+39%), Cr5Mn5 no crossing |
| dp3 | co-dopant by mean field / cross-coupling | Nb instead of Mo |
| dp4 | T_C linear in concentration (V2Cr2) | 4.3 K (> 3 K) |
| dp5 | zero-fill the missing Co–Co dimers | a number instead of null |
| dp6 | near shells (0–6) only | no order for V5Cr5, Cr5Mo5 |
| dp7 | read m2 ≈ 0.01 μB as moment collapse; drop Mn/Mo | Cr 5% alone, 5.3 K |

## Calibration (full table in `design.json` → `calibration`, also `calibration.json`)

The reference T_C is the mean of the scaled-magnetization crossings, [<m²>] L^(1+η) with η = 0.0375, for L = 4/5 and L = 5/6. The runs use 48/32/24/12 placements for L = 3/4/5/6 and 1500 + 5000 sweeps.

| Field | Reference | Defensible spread (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| V5Cr5 | 10.71 K | 9.81–11.37 K (Binder or scaled crossings, L = 3/4 … 5/6); oracle 12.16 | single placement 14.9 (+39%); near shells 0 (no order); MF 81.8; authors' rule 39.2 | ±20% (8.56–12.85) |
| Cr5Mo5 | 19.79 K | 17.37–19.58 K; oracle 21.68 | m2 misread: Cr 5% alone 5.35; near shells 0; MF 76.7; authors' rule 34.0; single placement 17.6 (passes this field, fails on V5Cr5/Cr5Mn5) | ±20% (15.83–23.75) |
| Cr5Mn5 | 10.84 K | 9.89–11.96 K; oracle 11.55 | single placement: no crossing; m2 misread 5.35; MF 60.2; authors' rule 22.3 | ±20% (8.67–13.01) |
| V2Cr2 | no order down to 0.8 K | noise-level Binder crossing 2.04 K (L = 4/5, g ≈ 1.65); oracle 0 | linear in concentration 4.28; MF 32.7; authors' rule 15.7 | ≤ 3 K |
| Cr5Co5 | null | oracle null | zero-filled Co–Co: no order (0); MF 50.6 | must be null |
| best co-dopant | Mo (19.8 K) | oracle Mo | Nb (mean field 96.6 K, authors' rule); Nb Monte Carlo ≤ 12.0 K (largest crossing + error over two reference sets), V 10.7, Mn 10.8 | exact |

Tightest margins:

- **V5Cr5:** the oracle's L = 3/4 crossing is +13.6% against the 20% tolerance (lab runs are shorter and smaller than the reference; the reference's own L = 3/4 value is +6%). On the wrong side, the single-placement path gives +39%.
- **Cr5Mo5:** the Binder crossing at L = 5/6 is −12%.
- **V2Cr2:** a noise-level Binder crossing at 2.04 K passes the ≤ 3 K bound; the linear-in-concentration path at 4.28 K fails by 1.3 K.
- **Best co-dopant:** Mo at 19.8 K against at most 12.0 K for Nb.

Rounding to 4 significant digits passes.

## Oracle

The oracle uses 34 evaluations and 54 of 60 units (90%). The catalog holds 276 units of pair windows, plus Monte Carlo.

- **Pair purchases:** both windows of Cr–Cr, V–V, V–Cr, Nb–Nb, Cr–Nb, Mo–Mo, Cr–Mo, Mn–Mn and Cr–Mn (18 units), plus Co–Co shells 0–6 (1 unit). The Co–Co window returns 0 rows, so the answer is null.
- **Monte Carlo, per composition** (V5Cr5, Cr5Mo5, Cr5Mn5, Cr5Nb5, V2Cr2):
  1. a coarse L = 3 scan with 16 placements (1 unit);
  2. L = 3 with 32 placements (2 units) on a grid around the Binder ratio g = 1.2;
  3. L = 4 with 32 placements (4 units).
- **T_C:** the crossing of [<m²>] L^(1+η). It is accepted only if the larger cell's Binder ratio falls below 1.12 at the lowest temperature. Cr5Nb5 stays at about 1.2, so it is not ordered.
- **Runtime:** 1133 s of container time in the gate, on 4 CPUs.
- **Results:** 12.16 / 21.68 / 11.55 K; V2Cr2 0; Cr5Co5 null; best Mo.

## Open risks

- **Variance of agent Monte Carlo.** The lab's sweep counts (1000 + 3000) and the size caps make L = 3/4 crossings run about 5–14% high. An agent that uses L = 2/3 or few placements can exceed 20%. This is the intended disorder and finite-size judgement (dp2), but it adds variance. Larger sizes (5, 6) cost 4–8 units per 16 placements.
- **No-order thresholds.** The V2Cr2 check (≤ 3 K) tolerates noise-level Binder crossings near 2 K. An agent that trusts such a crossing at 2.5–3 K still passes, which is lenient by design.
- **Cr5Nb5 is deliberately not graded as a number.** Its frustrated, weak order is not uniquely recoverable at affordable sizes. Independent reference sets disagree: no crossing in one, 10.6 → 7.7 K drifting in the other. Only the ranking (Mo ≫ Nb) is graded.
- **Model-derived truth.** The targets are consequences of the archived table under the stated model, not archived numbers. The implementation is pinned by reproducing the authors' mean-field matrix exactly, by an exact dimer-energy check of the heat-bath kernel, and by an identity test of the one- and two-species code paths.
- **Dev-easy risk.** An agent that buys both windows for every pair, notices the Co–Co gap, and runs careful finite-size scaling for each candidate will pass. The decision points that can be discovered (m2, missing Co–Co, Nb frustration, percolation at 2% + 2%) are visible only in purchased rows or Monte Carlo output.

## Disk

- **Peak:** about 11 MB: 5.2 MB in `~/.cache/ai-theorist/ti-codoping-exchange-tc` (2.9 MB of archive files, 2.3 MB of reference and scan JSON, the compiled kernel) and 5.4 MB of scratch (the paper PDF and logs), plus 1.3 MB of gate run directory.
- **Task package:** 0.76 MB (gate count; 0.9 MB on disk).
- **After the build:** the raw downloads, the scan files, the compiled kernel, the paper PDF and the gate's lab image were deleted. `fetch.py` refetches the downloads with pinned md5s. The reference Monte Carlo JSON (2.0 MB, in the cache's `ref/` folder) was kept so that `calibrate.py` can be rerun without the 4 h of reference runs.

## Results

- **Gate (2026-10-21, cpuset 4,5,20,21):** passed. 7 decision points; all 8 ablations fail; the 5 robustness variants and 4-significant-digit rounding pass; 9 rubric items; 3 hidden values with no leaks; 0.76 MB package. The oracle made 34 evaluations, charged 54 of 60, with 0 failures, in 1133 s of container time. A final `--reuse-state` re-grade followed the last design.json edit.
- **Harbor oracle:** reward 1.0 (`runs/build-ti-codoping-exchange-tc/oracle-r1-ti-codoping-exchange-tc`, 28.6 min, charge 54). It reproduced the gate's submission exactly, because the lab's placements and RNG are seeded deterministically.
