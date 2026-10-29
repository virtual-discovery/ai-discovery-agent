# Construction of skutterudite-wigner-coherence

This is a released-start task built from E. Di Lucente, M. Simoncelli and N. Marzari, Phys. Rev. Research 5, 033125 (2023). It uses Materials Cloud record 5jsfj-cdg70 (DOI 10.24435/materialscloud:xn-qr, CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md` and in `members.json`.

## Rebuild

1. Run `python3 -B fetch.py [DEST]`. It reads the zip's central directory and 87 single members (0.61 GB) with HTTP range requests (`remotezip.py`; no full download). It checks each member's CRC-32 and the sha256 in `members.json`. The default DEST is `~/.cache/ai-theorist/skutterudite-wigner-coherence/raw`.
2. Run `extract.py` with h5py, which is not in the runtime image. Use a throwaway virtualenv:
   ```sh
   python3 -m venv /tmp/h5env && /tmp/h5env/bin/pip install h5py numpy
   /tmp/h5env/bin/python -B extract.py
   ```
   It writes `pool/` next to `raw/`, about 15 MB:
   - per-run npz arrays;
   - `runs.json`, with the metadata, launch commands and parsed log tables;
   - `harmonic.json`, with the PDOS and MSD summed by element and the YbFe4Sb12 DFPT Γ frequencies.
3. Run `make_task.py` in the runtime image:
   ```sh
   docker run --rm --cpuset-cpus <cpus> --user $(id -u):$(id -g) -e HOME=/tmp \
       -v $PWD:/repo -v ~/.cache/ai-theorist/skutterudite-wigner-coherence:/cache \
       -w /repo/benchmarks/construction/skutterudite-wigner-coherence \
       ai-theorist-benchmark-runtime:20261007 python3 -B make_task.py /cache/pool
   ```
   It first asserts that `model.py` reproduces every archived `kappa_P_RTA` and `kappa_C` of the six phono3py 2.1.0 runs to 1e-6 (the observed agreement is 1e-8 to 1e-10). It then writes:
   - the lab data, which holds per-mode arrays and velocity operators but no archived κ arrays;
   - `catalog.json`;
   - the private reference;
   - `ablation_values.json`;
   - the budget in both `public.json` files.

   It prints the calibration numbers.
4. Run `.venv/bin/python -m benchmarks.tools.check_task skutterudite-wigner-coherence --cpuset <cpus>`. It writes `check.json`.
5. Run `python3 -B write_design.py`. It writes `design.json`, including the calibration table and the oracle length from `check.json`.

## What the task is

The agent buys four kinds of data:

- archived per-mode phono3py arrays per run and temperature (1 unit);
- velocity-operator products per run (4 units);
- archived output files of released runs (1 unit);
- harmonic summaries per compound (1 unit).

The budget is 32 units; the catalog is worth 202 units: 160 mode_data, 24 velocity operator, 11 files and 7 harmonic. The agent computes the Wigner conductivities itself.

The archived outputs of the Ir family (IrSb3, YbIr4Sb12) are withheld. The graded fields are:

- κ_P of YbIr4Sb12 at 300 and 700 K;
- κ_C of YbIr4Sb12 at 100, 300 and 700 K;
- T_half of YbIr4Sb12, where κ_C/κ_P = 0.5;
- the Yb rattler frequency (maximum of the Yb-projected DOS);
- κ_P of IrSb3 at 50 K and κ_TOT of IrSb3 at 300 K.

## What the archive shows (and the task is built on)

- **Formula.** Phono3py's Wigner κ_P sums each degenerate subspace as a block (off-diagonal velocity-operator elements within it included). κ_C is restricted to pairs with |Δf| ≥ 1e-4 THz, uses a Lorentzian of full width 2γ, and has a prefactor ħ × (κ_P units). These rules reproduce every archived table exactly.
- **Corrupted isotope rates (Yb-filled runs).** `gamma_isotope` has 17, 26 and 21 entries of 10³–10⁸ THz in YbCo4Sb12, YbFe4Sb12 and YbIr4Sb12. All of them sit at the irreducible q-point nearest Γ, (1/9, 0, 0), with weight 12; every other entry is below 0.01 THz.
  - In YbCo4Sb12 and YbIr4Sb12 the affected modes include the acoustic branches, so the archived κ_P is about 12% too low at 300–700 K and 22% too low at 50 K. In YbFe4Sb12 only optical modes are hit (0.2%).
  - With the corrupted entries kept, `model.py` reproduces the authors' κ_P(YbIr4Sb12, 300 K) = 0.816 W/mK. With them removed it gives 0.917 W/mK.
  - The other unfilled runs are clean. CoSb3 has 18 milder outliers (0.01–3.6 THz) at its first q-point, which do not matter (0.002%).
- **Inconsistent CoSb3 references.** In `v_operator_mean_free_path_isotope`, the phono3py 2.1.0 log gives K_P = 8.328 W/mK at 300 K, but `tc_data.dat` gives kP = 16.66. That file was written by `2_plot_MFP.py`, whose populations prefactor is 0.5 (instead of 0.25) and which divides C_s2 by ω_1.
- **Duplicate "isotope" files.** The phono3py 1.19 "isotope" folders (`CoSb3/tc_calc_with_isotope`, `YbFe4Sb12-NM/thirdorder_isotope`) hold kappa files byte-identical to the no-isotope runs. Only their logs carry the isotope results.
- **PDOS frequency axis.** All six phonopy PDOS runs used `--qe` (frequency factor 108.970772) on fc2.hdf5 files in eV/Å², where phono3py uses 15.633302. The `projected_dos.dat` axis is therefore inflated by 6.9704. The MSD runs used the default calculator and are correct: the Einstein frequency from the Yb MSD, 0.978 THz, matches the corrected Yb DOS peak at 0.968 THz.
- **YbFe4Sb12 fc3 source.** Finite differences (`thirdorder`, DFPT fc2) and hiPhive (hiPhive fc2 and fc3) differ by 2.4× in κ_P (0.568 vs 1.338 W/mK at 300 K), mostly through the linewidths (median ratio 2), but by under 10% in κ_C. Neither run is the paper's stated DFPT-fc2 + hiPhive-fc3 combination, so YbFe4Sb12 is not graded.

## Deviations from the design (candidates.json / DETAILS.md section 6)

- **No velocity-operator transfer.** The design withheld the velocity operator of the hidden compounds and graded κ_C from a model calibrated on analogues. Tested transfers to YbIr4Sb12 (frequency-binned |V|² maps, raw or normalised frequency) give κ_C ratios of:
  - YbCo4Sb12: 0.92–1.04;
  - CoSb3: 0.91–1.16;
  - IrSb3: 1.35–1.40;
  - YbFe4Sb12: 0.68–0.87;
  - FeSb3: 0.55.

  The "calibrate on CoSb3" wrong path is not wrong enough, and the defensible spread exceeds any useful tolerance, so the value is not uniquely recoverable. Linewidth transfer was tested too: κ_C is robust to it (0.95–0.97 from filled analogues), but κ_P is not (1.15–1.55). The velocity operators of the hidden runs are therefore served, and the archive's own defects (above) supply the difficulty.
- **Mesh decision point dropped.** The design assumed a CoSb3 9×9×9 run; the archive has CoSb3 only on 15×15×15.
- **Isotope decision point inverted.** The physical isotope effect is 0.6% for YbIr4Sb12 and 1.2% for IrSb3 at 300 K, so "omit isotopes, miss IrSb3 at 300 K" does not hold. The isotope fork is now:
  - the corrupted rates (dp1, dp6);
  - the need to add γ_iso where it is genuine, graded through IrSb3 κ_P at 50 K, where it is 6.3% (dp2).
- **Targets changed.**
  - The crossover is graded as T_half (κ_C/κ_P = 0.5) at ±25 K instead of ±75 K. The corrupted-isotope path lands 58 K low.
  - κ_C at 100 K was added to test degenerate handling.
  - IrSb3 κ_TOT at 300 K is kept as a sanity check.
  - IrSb3 κ_P at 50 K was added.
  - The rattler frequency was added for the PDOS-axis defect.
- **NM/AFM analogue choice dropped.** No κ data exist for AFM YbFe4Sb12, only DFPT Γ frequencies. Magnetism is a rubric point.
- **Lorentzian versus Gaussian, degeneracy and linewidth convention.** These are kept as the design's dp1, dp2 and dp6 (task dp4, dp5 and dp3).

## Decision points (7)

| dp | Fork | Kind |
|---|---|---|
| dp1 | Archived γ_iso of the Yb-filled runs is corrupted at the q-point nearest Γ; the archived YbCo4Sb12 table reproduces only with it | defect found in raw data; validation confirms the wrong path |
| dp2 | γ_iso must be added to γ (it matters for unfilled compounds at low T); the duplicate 1.19 "isotope" files suggest it does nothing | validation against logs; archive duplicate |
| dp3 | tc_data.dat (κ_P × 2) vs phono3py log for the lifetime convention | conflicting archived references; script bug |
| dp4 | Degenerate subspaces belong to the populations | validation (low-T κ_C) |
| dp5 | Lorentzian with full width 2γ vs Gaussian or half width | convention, validation |
| dp6 | Repair the corrupted entries vs discard the q-point | regime judgement (the q-point carries the long-wavelength acoustic heat) |
| dp7 | Units of the archived PDOS axis (phonopy `--qe` on eV/Å² constants) | metadata inconsistency; MSD cross-check |

## Calibration

The full table is in `design.json` under `calibration`.

| Field | Reference | Defensible variants | Ablations that move it | Tolerance |
|---|---|---|---|---|
| YbIr4Sb12 κ_P 300 K (W/mK) | 0.9167 | no γ_iso for the filled compound +0.6%; band-median repair −0.4%; diagonal populations −0.5% | archived γ_iso −11.0%; q-point discarded −14.6%; tc_data convention +100% | ±5% |
| YbIr4Sb12 κ_P 700 K | 0.3968 | +0.3%, −0.1%, −0.5% | archived γ_iso −10.8%; q-point discarded −14.3%; tc_data +100% | ±5% |
| YbIr4Sb12 κ_C 100 K | 0.1101 | no γ_iso −1.0%; threshold 0.1 cm⁻¹ −1.2% | degenerate pairs as coherences +10.3%; Gaussian −28.3%; half width −19.8% | ±6% |
| YbIr4Sb12 κ_C 300 K | 0.1784 | −0.4%, −0.3% | Gaussian −14.7%; half width −16.4%; degenerate +2.5% | ±6% |
| YbIr4Sb12 κ_C 700 K | 0.2094 | −0.1% | half width −10.3%; Gaussian −7.4% | ±6% |
| YbIr4Sb12 T_half (K) | 667.1 | −3.0 to +2.3 K (incl. separate interpolation +0.8) | archived γ_iso −58.2; q-point discarded −85.7; Gaussian +46.6; half width +65.8; tc_data +676 | ±25 K |
| YbIr4Sb12 rattler frequency (THz) | 0.968 | Yb-MSD Einstein +1.1%; lowest Γ optical mode −14.3% | PDOS axis read as THz +597% | ±20% |
| IrSb3 κ_P 50 K | 131.9 | diagonal populations −0.2% | isotope omitted +6.3%; tc_data +100% | ±4% |
| IrSb3 κ_TOT 300 K | 11.76 | −0.4% | tc_data +99%; isotope omitted +1.2% | ±5% |

The tightest margins are:

- isotope omitted against IrSb3 κ_P at 50 K (+6.3% against ±4%);
- degenerate pairs against κ_C at 100 K (+10.3% against ±6%);
- the Gaussian path against κ_C at 700 K (−7.4%), although the same ablation fails at 100 and 300 K by far more.

## Oracle and length

The oracle makes 18 purchases (27 of 32 units) in 22–40 s of container time. Its purchases are:

1. validation on YbCo4Sb12: log, mode data at 300 K and velocity operator;
2. the CoSb3 log and tc_data.dat;
3. the FeSb3 log and mode data at 100 K;
4. for YbIr4Sb12: the velocity operator and mode data at 100, 300, 600, 650 and 700 K;
5. for IrSb3: the velocity operator and mode data at 50 and 300 K;
6. harmonic data for YbIr4Sb12 and IrSb3.

It flags γ_iso entries above 1 THz (or above the mode frequency) and sets them to zero, and it rescales the PDOS axis by the ratio of the phonopy and phono3py conversion factors. An agent needs roughly 45–75 tool calls.

The gate passed on 2026-10-19 with all 8 ablations failing and all 7 robustness variants passing. The Harbor oracle scored reward 1.0 (`runs/build-skutterudite-wigner-coherence`).

## Open risks

- **The graded κ_P departs from the authors' archived value by 11%.** Graded κ_P values use the archived inputs with the corrupted isotope entries removed, not the authors' published κ_P (which includes them). The instruction defines the scattering physically (mass-variance isotope scattering), and any repair that keeps the q-point agrees within 0.4%. A reviewer could still argue that the authors' archived number is "the" answer.
- **Large responses.** The velocity-operator responses are 1.3 MB (5 MB for CoSb3) of JSON, and lab state.json grows to about 8 MB per episode.
- **Difficulty may be low.** Conventions (dp3–dp5) can be resolved by validating against released logs, which careful agents do. The discriminating forks are dp1, dp6 and dp7, which need physical judgement about archived numbers that reproduce the authors' tables.
