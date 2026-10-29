# Construction of cein3-kondo-rkky-magnons

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). The source is W. Simeth et al., Nat. Commun. 14, 8239 (2023). The data are in Zenodo record 10146787 (CC BY 4.0): 48 small `.dat` files, 1.0 MB in total.

## Rebuild

1. Run `python3 -B fetch.py`. It downloads the 48 files through the record's per-file API links into `~/.cache/ai-theorist/cein3-kondo-rkky-magnons/raw` and checks each md5 against `members.json`.
2. Run `calibrate.py` in the runtime image, on the builder's cpuset:
   ```sh
   docker run --rm --cpuset-cpus <cpus> --user $(id -u):$(id -g) -e CEIN3_RAW=/raw \
       -v <raw>:/raw:ro -v $PWD:/c -w /c ai-theorist-benchmark-runtime:20261007 python3 -B calibrate.py
   ```
   It takes about 6 minutes. It fits each model family to the purchasable data, predicts the withheld data with the same binning, and writes `calibration_fits.json`.
3. Run `python3 -B make_task.py` in the same image. It writes:
   - `catalog.json`;
   - `data.json`;
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - both `public.json` files;
   - `design.json`, including the calibration table and the oracle numbers from `oracle_result.json`.

   It also prints the calibration table.
4. Run `.venv/bin/python -m benchmarks.tools.check_task cein3-kondo-rkky-magnons --cpuset 4,5,20,21`. It writes `check.json`.

Supporting files:

- `cmodel.py`: the forward model. It is linear spin-wave theory for a simple-cubic antiferromagnet ordered at R, with real-space shells plus an optional long-range staggered term (a Gaussian minimum of J(q) at R). Each data bin is averaged over its momentum box, then energy-binned and broadened.
- `cfit.py`: the data sets, the variable-projection fits (per-set scale and offset solved linearly) and the metric that `task_lab.grade` uses.

## What the task is

The lab sells archived χ″ data (198 units in the catalogue; budget 60):

- single constant-Q spectra of the 12-meV maps (R–Γ, Γ–X, X–M, M–Γ), 1 unit each;
- single constant-q spectra of the high-resolution maps along R–Γ and R–X, 1 unit each;
- the four high-resolution constant-energy cuts, 2 units each;
- the band structure, 1 unit.

Withheld and graded:

- the high-resolution R–M map (35 × 51);
- the R–M cut at 1.4 meV (70 points);
- the Q-integrated χ″ (16 points).

Each prediction is scored by a reduced χ² after the best non-negative scale. The R–M map and cut also get a constant offset.

The submission also reports:

- W, η = v/(2πW) and J1 of the best nearest-neighbour model of the 12-meV data, each graded against a window;
- the exchange range, graded categorically.

## Deviations from the design (candidates.json / DETAILS.md section 20)

- **The lab serves data, not an LSWT calculator.** Simple-cubic LSWT is a few lines of numpy, so a calculator lab would not make the budget bind. Instead, the agent buys data. Choosing which spectra to buy (zone-boundary 12-meV spectra versus high-resolution data near R) is the main budget decision.
- **Graded targets are withheld archived data,** as the task notes ask. The design's targets were adapted as follows:
  - *Minimum energy near R within 0.15 meV and the local-minimum location within 0.03 r.l.u.:* not recoverable. The bins below about 0.25 meV are dominated by the subtracted elastic line. Free-gap fits give 0.01–0.43 meV, depending on the energy cut and the data set, and a fixed 0.35-meV gap changes no score by more than 0.02. These targets were replaced by χ² checks on the R–M map and the 1.4-meV cut.
  - *Q-integrated χ″ shape (correlation > 0.8):* kept, as a scale-only χ² check (≤ 5.5). It is a strong discriminator: every short-range LSWT model gives χ² ≥ 10, and long-range models give 0.4–3.3.
  - *Sign pattern of J1..J4:* dropped. It is not uniquely recoverable; two starts of a 10-shell fit give different signs and different bandwidths. That fit is now the dp7 ablation.
- **Extra graded scalars:** W, η, J1 (for S = 1/2) and the exchange-range category.
- **The R–X constant-energy cuts are served,** although the design withheld them. They are the only visible cuts besides R–Γ. The R–M 0.4-meV cut is withheld but not graded, because no model separates on it.
- **No damping decision point.** A resolution-only fit (no intrinsic Lorentzian) passes every graded check within the noise (for example η = 1.09–1.83, depending on the window reading), so it cannot be an ablation. Damping remains a rubric item: every good fit needs a HWHM of about 0.4–0.6 meV beyond the resolution.
- **Velocity is not graded directly.** It depends on the reading of the Methods phrase "within a distance of 0.08 r.l.u." (±0.08 gives v ≈ 22, a total width of 0.08 gives v ≈ 33 meV/r.l.u.). The paper's 38.6 meV/r.l.u. used cuts that are not archived. The η window [0.75, 3.0] accepts all of these readings, the crude cut-width estimate (0.83) and the paper's 2.23. It rejects:
  - J1 alone (0.577);
  - the overfit 10-shell model (0.61);
  - the model without the integration window (4.06).
- **Hidden-value scan.** The catalogue gives grids as formulas (`q = -0.33 + 0.02 * index`; `s = sqrt(3)/2 + index * 0.025`) and requests use integer indices. Listing several hundred decimal coordinates collided by chance with 4-significant-digit renderings of the 1,703 withheld χ″ values.

## Decision points (7)

| dp | Fork | Ablation result |
|---|---|---|
| dp1 | 12-meV data only, nearest-neighbour model | R–M map 2.05, cut 3.67, integrated 10.6, η 0.577 |
| dp2 | Short-range shells (J1–J3; 6 shells) | integrated 10.2 / 25.2; exchange_range short_range |
| dp3 | Perpendicular integration windows ignored | η 4.06, integrated 6.94 |
| dp4 | One J1 forced through all data | W 4.87, J1 1.62 |
| dp5 | Free-ion J = 5/2 as the effective spin | J1 0.198 |
| dp6 | Broad spectra read as an incoherent continuum | R–M map 2.05, cut 2.35 |
| dp7 | Ten free shells (overfit) | W 5.55, η 0.61 |

Where the difficulty comes from:

- dp1 is a budget decision: the 12-meV data alone are fitted well by J1.
- dp3 is evidence in the reduction metadata. The ±0.08 r.l.u. perpendicular window turns a steep cone into broad, plateau-like constant-q profiles that look like a continuum (dp6) or a gap.
- dp2 is model class plus transfer: short-range fits pass the R–M checks but mispredict the integrated χ″.

## Calibration

The full table is in `design.json` under `calibration`.

| Field | Oracle | Defensible spread | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| R–M map χ²_r | 1.123 | 1.06–1.10 (ref. 1.07) | continuum 2.048; J1 2.053 | ≤ 1.70 |
| R–M 1.4-meV cut χ²_r | 1.349 | 1.33–1.45 (ref. 1.31) | continuum 2.35; J1 3.67 | ≤ 1.90 |
| integrated χ²_r (scale only) | 0.399 | 0.43–3.26 (ref. 0.39) | no window 6.94; J1–J3 10.2; J1 10.6 | ≤ 5.5 |
| W (meV) | 2.644 | 2.61–2.98 (paper 2.75) | 6 shells 3.39; forced J1 4.87; 10 shells 5.55 | [2.35, 3.25] |
| η | 1.359 | 0.83–1.93 (paper 2.23) | J1 0.577; 10 shells 0.613; no window 4.06 | [0.75, 3.0] |
| J1 (meV, S = 1/2) | 0.981 | 0.917–0.992 | free ion 0.198; forced 1.62 | [0.8, 1.2] |

Tightest margins:

- **η lower edge:** the crude box-width estimate (0.834) passes by 0.08; the 10-shell overfit (0.613) fails by 0.14.
- **Integrated χ″:** the window-ignored model fails by 1.44.
- **R–M map:** the continuum fails by 0.35.
- **W:** the 6-shell fit fails by 0.14. It also fails on the integrated χ″.

The flat (continuum) integrated χ″ passes (2.52). It is caught by the R–M checks instead.

## Oracle

The oracle makes 49 evaluations for 53 of 60 units (88%):

- 20 12-meV spectra: 5 per path, from mid R–Γ to the zone boundary;
- 24 high-resolution spectra: 12 per line, within |q| ≤ 0.25;
- the 4 cuts;
- the band structure.

It fits:

- J1 to the 12-meV data, for J1;
- J1–J3 to all purchases;
- J1 plus a long-range term, and J1–J3 plus a long-range term, keeping the better of the two.

It then predicts the withheld data. Its scores are 1.123, 1.349 and 0.399, with W = 2.644, η = 1.359 and J1 = 0.981. The container time is 22 s.

Results:

- **Gate (2026-10-19):** passed; all 8 ablations fail, 5 robustness variants pass, and rounding to 4 significant digits passes.
- **Harbor oracle:** reward 1.0 (`runs/build-cein3-kondo-rkky-magnons`, 1 min 5 s, 53 units).

## Open risks

- **Memorisation.** The paper's W = 2.75 meV and η = 2.23 fall inside the graded windows, and its title gives away the verdict. The three predictions require modelling the data with the binning.
- **Model dependence of the χ² references.** All calibration models share one family: LSWT plus a Gaussian long-range term. An agent with a different long-range form (for example an RKKY J(r) with a fitted k_F) might give somewhat different χ² values. The thresholds sit 0.2–0.6 above the defensible spread.
- **Noise model of the map.** There are no archived errors. σ = 0.266 comes from the ±q symmetry and the map score uses it. The threshold was calibrated with the same σ.
- **Dev-easy risk.** An agent that fits all the purchasable high-resolution data with a flexible model and models the integration windows will pass. The difficulty rests on doing both under the budget, and on realising that short-range fits mispredict the integrated χ″.
- **Peak disk use:** about 3 MB in total (1.1 MB of raw files, 0.2 MB of calibration output, 1.6 MB of gate and Harbor run directories). The raw files were deleted after the build; `fetch.py` and `members.json` remain.
