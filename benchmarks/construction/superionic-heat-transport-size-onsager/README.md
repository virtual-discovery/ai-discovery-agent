# Construction of superionic-heat-transport-size-onsager

Benchmark task, released start (see `benchmarks/tools/BUILDING.md`). Source: F. Grasselli, J. Chem. Phys. 156, 134705 (2022); Materials Cloud record yzhmw-vje16 (CC BY 4.0), the same record as the sibling `fluorite-superionic-finite-size`.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` (`docker run --rm --cpuset-cpus 6,7,22,23 --user $(id -u):$(id -g) -v $PWD:/c -v ~/.cache/ai-theorist/superionic-heat-transport-size-onsager:/cache ...`), except the network fetch.

1. `taskset -c 6,7,22,23 python3 -B fetch.py` reads 511 members (248 thermo logs + decks + scripts) of the heat folders out of `UO2_heat.zip`, `CaF2.zip` and `PbF2.zip` by HTTP range requests, checks CRC32 and the SHA256 in `members.json`, and stores each log as a float32 `.npy` of its first 12 columns (about 1.2 GB of .npy under the cache; the 248 logs are 2.3 GB of text, transferred as deflated byte ranges). Retries cover the archive's transient 503s and truncated range responses.
2. `extract.py` (container) writes `pool.npz` (int16 series, 0.08 ps, band-limited to 6.25 THz) and `pool_meta.json`.
3. `analysis.py` (container) writes `kappa_table.json`: per-run κ from the reference (multicomponent cepstral, full served band, AIC) and variants (3 THz low-pass, raw single-flux, direct multicomponent running integrals).
4. `make_task.py` (container, with the task folder mounted at `/task`) writes the catalog, `data.json`, `series.npz` (241 served runs), the private targets, `hidden_values.json` and `solution/ablation_values.json`, and prints the calibration.
5. `.venv/bin/python -m benchmarks.tools.check_task superionic-heat-transport-size-onsager --cpuset 6,7,22,23`.

`cepstral.py` is the estimator (numpy re-implementation of the thermocepstrum multicomponent algorithm). It was validated on synthetic AR processes with known S(0): single-current and two-current (Schur-reduced) estimates are unbiased within 1%. Two pitfalls found while building it, both relevant to agents: (i) removing the mean before the FFT zeroes the f = 0 periodogram bin and biases log-spectrum methods low by ~10%; (ii) the Schur-reduced periodogram of the component-averaged 2×2 matrix is S·χ²₄/6 (complex Wishart), not S·χ²₄/4, and the cross-spectrum must be kept complex. `explore.py` is the full-resolution exploration used to choose the design.

## Design as built

- **Request:** `md_run` with a run id `<material>-<NVE|NVT>-<n>x<n>x<n>-<T>K`. Cost 1/2/4 for 96/324/768 atoms; budget 48; catalog total 544 (241 runs).
- **Returned:** deck, Lx, mean/std temperature, mean pressure, total-energy drift, and six series (heat flux J·V of `compute heat/flux` and the anion velocity sum, 3 components each) as base64 int16 with scales.
- **Withheld states (no run served):** UO₂ 4×4×4 1700 K and 3400 K, CaF₂ 4×4×4 1600 K (NVE+NVT each), PbF₂ 4×4×4 1200 K (NVT only; PbF₂ has no 4×4×4 NVE run above 900 K). Target: mean multicomponent κ of the withheld runs.
- **Graded:** 4 predictions (±12%), κ of 3 served superionic 768-atom runs (±12%), thermostat verdict (false).

## Findings from the released data (multicomponent cepstral κ, W/m/K)

- Raw single-flux GK equals the multicomponent value in the crystal (UO₂ ≤ 2000 K: ratio 0.95–1.05) and is inflated once the anions diffuse: UO₂ 1.15–1.25 above 3000 K, CaF₂ 1.4–2.4 above 1400 K, PbF₂ 3–4.2 above 900 K. The raw value is gauge dependent in the superionic regime: adding c·V_anion to the flux (a change of the per-species energy zero) changes the raw UO₂ 3200 K value from 1.17 to 3.35 W/m/K for |c| ≤ 4.5 eV, while the multicomponent value is invariant.
- κ(T) of 768-atom cells falls far more slowly than 1/T and saturates: UO₂ 1.67 (1700 K) → ≈1.2 above ~2800 K; CaF₂ 2.1 (500 K) → ≈1.35 above ~1400 K; PbF₂ 1.15 (300 K) → ≈0.65 above ~900 K.
- Size: 96-atom cells are 15–30% low; 324 atoms are within 3–10% of 768 atoms except crystalline UO₂ (−15% at 1700 K).
- NVT (CSVR, τ = 0.2/0.4 ps) vs NVE: 116 pairs, κ_NVT/κ_NVE = 0.998 ± 0.064, consistent with the 3–5% per-run errors.

## Deviations from the design (and why)

- **(b) saturation temperature not graded.** With 3–5% noise per run the onset is located only to ~200–300 K; rubric item.
- **(c) mass-flux correction not graded as a number.** Asking for the raw-flux κ would name decision point 1. It is graded implicitly: κ of three served superionic runs, where the raw-flux path fails by +22%, +118% and +325%; the size of the correction is a rubric item.
- **(d) size-extrapolated κ not graded.** Three sizes, 3–5% noise, and 1/L vs 1/N forms differing by >20%: not uniquely recoverable. Graded states are 768-atom states; the rubric asks the report to say the limit is unsupported.
- **NVE vs NVT** is a fork about evidence, not a bias: the global CSVR thermostat does not bias κ in these runs, so the graded item is the verdict (and the NVT-only PbF₂ state).
- **CPU budget/compute:** all served data are precomputed lookups; no runtime dependencies beyond numpy/scipy.
- α-AgI heat runs (NVT only, no NVE) were not used.

## Calibration

Full table in `design.json` → `calibration`. Relative errors versus the hidden reference:

| Graded field | Reference | Oracle | Defensible variants | Nearest wrong paths | Tol |
|---|---|---|---|---|---|
| UO₂ 4×4×4 1700 K | 1.672 | −1.6% | −5.2 … +2.9% | 324-atom −15.2%; 96-atom −25.5% | ±12% |
| UO₂ 4×4×4 3400 K | 1.211 | +0.3% | +0.2 … +2.1% | 96-atom −15.1%; raw flux +21.1%; 1/T −27.9% | ±12% |
| CaF₂ 4×4×4 1600 K | 1.386 | +1.0% | −1.3 … +2.6% | 96-atom −11.8% (passes); 1/T −25%; raw +90% | ±12% |
| PbF₂ 4×4×4 1200 K | 0.689 | −4.6% | −9.2 … −1.3% | NVE 3×3×3 −11.8% (passes); 96-atom −19.1%; 1/T −43%; raw +270% | ±12% |
| κ UO₂-NVE-4x4x4-3500K | 1.196 | 0 | −1.6% | raw +22.5% | ±12% |
| κ CaF₂-NVE-4x4x4-1700K | 1.379 | 0 | −3.6% | raw +118% | ±12% |
| κ PbF₂-NVT-4x4x4-1300K | 0.626 | 0 | −1.3% | raw +325% | ±12% |
| thermostat verdict | false | false | — | true | exact |

Tightest margins: the 324-atom UO₂ ablation fails by 3.2 points; the raw-flux path on UO₂ 3400 K fails by 9 points (it is caught much more strongly on CaF₂/PbF₂); the PbF₂ single-neighbour variant passes with 2.8 points to spare. The direct multicomponent running-integral estimator scatters ±10% per run (it is ill-conditioned in the crystal, where the anion current is non-diffusive) and is not listed as a robustness variant; on the PbF₂ 1300 K run it gives +14%.

## Oracle

13 purchases, 43 of 48 units: UO₂ 4×4×4 NVE+NVT 1800 K, NVE 2600/3300/3500 K; CaF₂ 4×4×4 NVE 1500/1700 K; PbF₂ 4×4×4 NVT 1100/1300 K; UO₂ 2×2×2 and 3×3×3 NVE 1700 K; PbF₂ 3×3×3 NVE+NVT 1300 K. Gate container time 2 s. Gate (2026-10-26): passed, 6/6 ablations fail, 4/4 robustness variants pass, rounding passes. Harbor oracle: reward 1.0 (`runs/build-superionic-heat-transport-size-onsager`).

## Open risks

- **Dev-easy risk.** A model that knows the Bertossa–Grasselli multicomponent GK method and checks NVT against NVE will solve it; difficulty rests on discovering the raw-flux inflation (and that it is absent in the crystal, which invites the wrong generalisation), on buying 768-atom cells, and on using NVT where NVE is missing.
- **Noise.** Single-run κ carries 3–5% statistical error; the UO₂ 1700 K target is the mean of two runs (1.736 NVE, 1.608 NVT).
- **Estimator choice.** Agents implementing the reduced-spectrum statistics with the wrong χ² normalisation (S·χ²₄/4) would be 33% low in the superionic regime and fail; this is a genuine implementation pitfall, not counted as a decision point.
- **Response size.** Each purchase returns ~160 kB of JSON; the instruction tells the agent to redirect to a file.
- **Disk.** Peak use about 1.25 GB under the cache (raw logs + pool); raw logs and pool deleted after the build, refetch scripts and checksums kept.
