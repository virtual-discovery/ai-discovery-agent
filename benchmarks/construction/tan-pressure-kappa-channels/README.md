# Construction of tan-pressure-kappa-channels

A released-start task built from A. Kundu et al., PRL 132, 116301 (2024), Zenodo record 10564835 (CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md` and `members.json`.

## Rebuild

CACHE = `~/.cache/ai-theorist/tan-pressure-kappa-channels`, CPUSET = your cpuset. Every step except 1, 4 and 6 runs in `ai-theorist-benchmark-runtime:20261007` with `--user $(id -u):$(id -g) -e HOME=/tmp -e OPENBLAS_NUM_THREADS=1`.

1. `python3 -B fetch.py` downloads the 25.5 MB tarball (md5 checked), extracts it to CACHE/raw and verifies the 43 members.
2. Four-phonon rates (about 25 min per run on one core; run one container per core):
   ```sh
   docker run --rm --cpuset-cpus <one cpu> ... -v $REPO/benchmarks/construction/tan-pressure-kappa-channels:/repo \
       -v CACHE/raw/IFCs-data-pressure-TaN:/raw -v CACHE/fourph:/s ai-theorist-benchmark-runtime:20261007 \
       python3 -B /repo/fourph.py P 16 0.1 8192 SEED T
   ```
   for (P, SEED, T) in (0,1,300) (0,2,300) (20,1,300) (40,1,300) (60,1,300) (80,1,300) (80,2,300) (100,1,300) (20,1,600) (60,1,600) (100,1,600). The sampler uses numpy's default_rng(seed), so the runs are reproducible.
3. `make_task.py` (about 11 min) with `-v $REPO:/repo -v CACHE/raw/IFCs-data-pressure-TaN:/raw -v CACHE/fourph:/s -w /repo/benchmarks/construction/tan-pressure-kappa-channels`: stage `data` writes the lab data, `tanbte.py` copy, catalog and both `public.json`; stage `reference` (`reference.py`) runs every reference, oracle-route, robustness and ablation evaluation through `task_lab._run` and writes `calibration_raw.json`.
4. `python3 -B finalize.py` writes `private/hidden.json`, `hidden_values.json`, `solution/ablation_values.json`, `calibration.json`.
5. Gate: `.venv/bin/python -B -m benchmarks.tools.check_task tan-pressure-kappa-channels --cpuset CPUSET`.
6. `python3 -B write_design.py` writes `design.json` (calibration and oracle length from `check.json`).

## The evaluator (tanbte.py)

Numpy implementation of the ShengBTE/FourPhonon protocol for the 2-atom cell: phonopy-format FC2 (atom-major, x-fastest supercell order, minimum-image weights), ShengBTE FC3, FourPhonon FC4, adaptive Gaussian broadening with the 2σ cut, three-phonon rates and ShengBTE's iterative Δ term, Tamura isotope rates, four-phonon rates (½, ½, 1/6 permutation factors, ħ²π/8 prefactor) and RTA/iterative κ with symmetry (12 point operations + time reversal: 270 irreducible points on 16³, as in the archive). A 16³ three-phonon evaluation takes ~25 s; a lab `kappa` call 30–60 s on 2 cores.

Checks against the archive: the irreducible list reproduces the archive's 270 points; the archived `BTE.sr_phel_metal` rows are branch-major over the ShengBTE irreducible list scanned with the third grid index fastest. Their frequency column matches our phonons to 1.1–1.3% median (max ~11%) at every pressure, with no pressure swapped; the residual is a different phonon calculation (the EPW phonon–electron workflow). The authors' κ outputs are not archived, so no direct κ check is possible; our released full-physics κ_xx (1081, 1509, 1904, 1812 W/mK at 0/20/60/100 GPa) reproduces the paper's qualitative non-monotonic shape (maximum near 60 GPa) at the published order of magnitude.

Four-phonon rates use unbiased Monte Carlo sampling of (q′, q″) pairs on the full 16³ grid (8192 pairs per irreducible q, all 6⁴ branch combinations exact). Seeds 1 and 2 at 80 GPa give κ_xx = 1879.3 and 1877.7 W/mK (0.09%); the lab serves the seed average where two exist.

## What the archive shows (and the task is built on)

- Four-phonon scattering is essential: with three-phonon + isotope only and the archived scalebroad = 0.1, hundreds of modes have no allowed process and the iterative solution diverges (κ_xx ~ 10⁴–10⁶ W/mK); with scalebroad 1.0 it converges at ~2200–3200 W/mK. Adding the four-phonon rates gives ~1100 (0 GPa) to ~1900 W/mK (60 GPa).
- Phonon–electron rates are negligible at 0–20 GPa (max 0.76 ps⁻¹) and jump at 40 GPa (max 14.5 ps⁻¹, optical branch means ×6); they make κ non-monotonic.
- The phonon–electron table's frequency column is from another phonon calculation and its row order is ShengBTE's, not the natural grid order: matching by frequency or reading rows in the catalog order are both wrong paths.
- Iterative κ is ~1.9× the RTA value (strong Normal scattering).
- The 4th-order constants do not satisfy the acoustic sum rule (residual 60–96 eV/Å⁴); as in FourPhonon, they are used as archived.

## Deviations from the design (candidates.json / DETAILS.md section 5)

- **Withheld set.** Design: 4th IFCs and ph-el tables hidden at 40 and 80 GPa. Built: 4th IFCs hidden at 40 and 80 GPa, ph-el hidden at 80 GPa only (40 GPa is released). Interpolating the ph-el rates across the 20→40 GPa onset is not uniquely recoverable (the rates jump by a factor ~6), so it would violate uniqueness; at 80 GPa the defensible spread is ±6% in κ.
- **Targets.** (1) κ_xx and κ_zz at 80 GPa (±10%) and κ_xx at 40 GPa (±5%; zz at 40 GPa dropped, it adds nothing beyond xx); (2) κ_xx at 60 GPa / 600 K (±4%) with the 60 GPa 600 K four-phonon rates withheld (600 K rates are served at 20 and 100 GPa); (3) the ph-el resistivity share at 80 GPa (±0.06 instead of ±0.07).
- **4ph lookups** are our own Monte Carlo FourPhonon-convention calculation (no Fortran FourPhonon was available) on the archived 16³ mesh only; no coarser-mesh 4ph lookups. Mesh convergence of 3ph remains explorable (meshes 8–20), but the graded protocol is the archived one, so mesh convergence is not a graded fork; the broadening setting is (dp6).
- **Units fork** kept (dp4) and a row-order fork added (dp5); RTA/iterative (dp1), channels (dp2, dp3), interpolation strategy (dp7), anisotropy (dp9) kept; temperature transfer of 4ph added (dp8).

## Calibration

| Field | Reference | Oracle route | Robustness | Nearest ablations | Tolerance |
|---|---|---|---|---|---|
| κ_xx 80 GPa 300 K | 1878.4 | −6.4% (linear in P) | geometric +5.5%, quadratic −3.9% | scalebroad 1.0 −17%; 2π −19%; row order −22%; no ph-el +31% | ±10% |
| κ_zz 80 GPa 300 K | 1536.8 | −4.2% | +3.7%, −1.8% | trace average +15%; scalebroad −16% | ±10% |
| κ_xx 40 GPa 300 K | 1835.3 | −1.4% | −0.6% | 2π −6.1%; κ-level interpolation −7.0%; scalebroad −14.5% | ±5% |
| κ_xx 60 GPa 600 K | 629.7 | +0.6% | −1.0%, +2.2% | 2π −5.0%; scalebroad −7.2%; T² scaling +21%; κ-level +51% | ±4% |
| ph-el share 80 GPa | 0.235 | +0.049 | −0.037, +0.035 | 2π +0.14; row order +0.17; no ph-el −0.235 | ±0.06 |

Tightest margins: oracle κ_xx(80) at −6.4% vs ±10%; ph-el share oracle +0.049 vs ±0.06; κ_xx(600 K) 2π path at −5.0% vs ±4% (the same path fails κ_80 by 19%); κ_xx(40) 2π path at −6.1% vs ±5%. The κ-level interpolation path lands within 1% on κ_80 and fails only on κ_40 (−7%) and κ_600K.

## Oracle and length

The oracle makes 17 evaluations (22 of 30 units) in ~200 s of container time: CONTROL, three ph-el tables (40/60/100), one harmonic set (60 GPa), five `rates` (four-phonon at 20/60/100 GPa and 300 K, 20/100 GPa at 600 K), and seven `kappa` (validation of the ph-el mapping at 60 GPa against the lab's own channel, then the targets). An agent needs roughly 45–75 tool calls. Gate passed 2026-10-26; Harbor oracle reward recorded in `runs/build-tan-pressure-kappa-channels`.

## Open risks

- The 80 GPa interpolation is the largest uncertainty; a different defensible transfer (e.g. interpolating in a descriptor) could land outside ±10%. The spread of the three tested routes is −6.4% to +5.5%.
- References are our implementation, not the authors' FourPhonon outputs (not archived); absolute values may differ from the paper by several percent.
- Difficulty: the row-order and units forks are discoverable but require comparing the archived frequency column with lab phonons; the lab's error message about the 4ph protocol (mesh 16, scalebroad 0.1) partly reveals dp6.
- θ-TaN is known (PRL 2024, Science 2025); graded values are protocol-specific new calculations at withheld pressures/temperatures.
