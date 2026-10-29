# Construction of gan-hole-piezomobility

This is a released-start task. It is built from Chen, Leveillee, Van de Walle and Giustino, *Design of high-mobility p-type GaN via the piezomobility tensor* (submitted 2025), using Materials Cloud record qfdfa-cs517 (DOI 10.24435/materialscloud:zy-qw, CC BY 4.0).

## Rebuild

1. Run `bash fetch.sh [DEST]`. It downloads `strain_GaN.tar.gz` (12.9 MB, md5 9a8747242f4759d235c5552aaea0ba52) and `README.txt` through the archive's per-file links. It then unpacks the tarball (about 150 MB) and verifies the 148 members listed in `members.md5`. The default DEST is `~/.cache/ai-theorist/gan-hole-piezomobility`.
2. Run `python3 -B extract.py [DEST/x/cloud]`. It writes two files:
   - `pool.json` (260 KB): cells, positions, the valence-band maximum, the zero-field SERTA and IBTE tables at 100–500 K, and the solver logs of all 37 states;
   - the staged `environment/lab/bands.json` (3.5 MB): GW bands 47–56 on the 10×10×10 grid.

   Re-running it reproduces both files byte for byte.
3. Run `python3 -B make_task.py`. It writes:
   - `catalog.json`;
   - `lab/data.json`;
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files.

   It also prints the calibration numbers. `model.py` holds the shared estimators.
4. Run `python3 -B write_design.py`. It writes `design.json`, including the calibration table and the oracle runs from `oracle_run.json`.
5. Run `.venv/bin/python -m benchmarks.tools.check_task gan-hole-piezomobility --cpuset <cpus>`. It writes `check.json`.

## What the task is

The agent buys archived results for 37 states: unstrained GaN, and compress or expand 1, 2 and 3% along each of e1–e6. For each state it can buy:

- GW bands: 1 unit;
- SERTA mobility tensors: 1 unit;
- IBTE mobility tensors, which also return the SERTA starting point and the solver log: 3 units.

The budget is 50 units. Buying every band set and every IBTE costs 133.

The mobilities of five states are withheld; their bands are served. The agent must predict the full 300 K IBTE tensors of those five states. It must also report two design answers as Voigt strains, together with the mobilities:

- the 2% state with the largest μ_zz;
- the 2% state with the largest in-plane principal mobility.

## Deviations from the design (candidates.json / DETAILS.md section 9)

- **GW eigenvalues are on the 10×10×10 coarse grid, not a 1000-point path.** The lab serves the grid. It does not serve masses: band character has to be read from the dispersion, for example along Γ–A.
- **Linear piezomobility coefficients are not graded (design target 2).** The response is non-analytic at zero strain:
  - HH and LH are 7 meV apart at Γ, below kT;
  - the e3 response is one-sided (+15 cm²/Vs at −1%, +0.3 at +1%);
  - the e1 response is asymmetric even at ±1%.

  A "linear coefficient" is therefore not uniquely defined by the data. The failure of linearity is a decision point instead (dp3).
- **Withheld states changed.** I kept `compress_3%_e3` and `3% e4` (as `expand_3%_e4`).
  - `expand_3%_e1` became `expand_3%_e2`. The e2 series supports a defensible extrapolation, while the e1 tension series fails linear extrapolation by +12.7%.
  - `compress_2%_e2` was dropped. No route recovers it within 10%: the swap of compress_2%_e1 gives −13% on μ_xx, and interpolation of the e2 series gives −8.8% on μ_yy.
  - `expand_3%_e3` (one-sided regime) and `compress_3%_e6` (in-plane shear, off-diagonal sign) were added.
- **`compress_3%_e3` is graded with a regime window [150, 280] cm²/Vs (truth 213.6), not ±10%.** Its exact value is not recoverable from the purchasable data. Past the band crossing, the 2→3% increment lies between 0 (163.9) and the 1→2% increment (264.7). A two-population Boltzmann model gives about 178. The window rejects the pre-crossing and linear-response paths (70–94) and a quadratic extrapolation (351).
- **Categorical target kept as "best 2% strain for μ_zz", reported as a Voigt vector.** A second design answer was added: the best in-plane principal mobility among the 2% states.
- **The GW-vs-DFT fork was dropped.** Only GW energies are archived, and recomputing DFT bands for 37 SOC states was not justified. It is replaced by three forks found in the data:
  - the shear-label convention (dp7);
  - the SERTA survey trap (dp5);
  - the unconverged-flag judgement (dp8).

## Decision points (8) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | Exact mirror partner vs extrapolating the same-sign shear series; flip of the strained off-diagonal | symmetry to recognise; the partner is purchasable but not labelled as such |
| dp2 | e1↔e2 swap valid only at linear order (no 90° rotation in 6mm); expand_3%_e1 swapped is −18%/+17% off | textbook symmetry failing in the nonlinear regime; found by buying both |
| dp3 | Linear piezomobility tensor from ±1% fails everywhere | model class / regime |
| dp4 | c-axis band crossing between 1% and 2% (GW bands, 2% mobility) | regime; discriminating purchase |
| dp5 | SERTA survey ranks compress_2%_e4 above compress_2%_e3 for μ_zz; IBTE reverses it (IBTE/SERTA 2.1 vs 2.6) | budget allocation; consequence visible only after buying IBTE |
| dp6 | In-plane maximum needs the principal value (e6: about 205 along about 52°), not the diagonal (e1: 158) | off-diagonal elements |
| dp7 | Archive "X% e4/e5/e6" applies ε_ij = X%, so Voigt e = 2X% (visible in the cells only) | convention found in metadata |
| dp8 | 15 IBTE runs print "reached the maximum but did not converge"; residual ≤ 3.2e-4 cm²/Vs per iteration; drift < 0.06% after iteration 100 | data-quality judgement from raw logs |

## Calibration

Values are relative errors unless marked as values; the full table is in `design.json` under `calibration`.

| Field | Truth | Oracle | Defensible alternatives | Nearest wrong path | Tolerance |
|---|---|---|---|---|---|
| compress_3%_e3 μ_zz | 213.6 | 208.3 (value) | 163.9, 264.7 (values) | 93.9 (pre-crossing), 70.3 (linear) | [150, 280] |
| expand_3%_e3 μ_xx = μ_yy | 62.58 | −0.010 | 0.000 (flat from 1%) | +0.168 (scaled SERTA), +0.166 (linear) | ±10% |
| expand_3%_e3 μ_zz | 47.47 | +0.022 | +0.010 | −0.475 (linear) | ±10% |
| expand_3%_e2 μ_xx | 155.47 | −0.017 | +0.059 (linear series), −0.048 (ratio) | −0.176 (swap), −0.125 (linear) | ±10% |
| expand_3%_e2 μ_yy | 47.41 | −0.014 | −0.030, −0.035 | +0.168 (swap) | ±10% |
| expand_3%_e4 μ_xx / μ_yy / μ_zz | 58.26 / 134.13 / 185.81 | −0.008 / 0.000 / +0.010 | none | yy +0.187, zz +0.137 (same-sign series); yy +0.47 (scaled SERTA) | ±10% |
| expand_3%_e4 μ_yz | −49.06 | +0.059 | none | +0.21 (series), −2.06 (sign not flipped) | ±15% |
| compress_3%_e6 μ_xx / μ_yy / μ_zz | 117.45 / 166.28 / 52.23 | +0.014 / +0.011 / +0.013 | +0.074 / +0.059 / +0.035 (linear series) | linear tensor −0.53 / −0.68 | ±10% |
| compress_3%_e6 μ_xy | 81.48 | −0.001 | +0.094 | −2.00 (sign not flipped), +1.02 (linear) | ±15% |
| best_mu_zz_2pct | e3 = −0.02, μ_zz 163.88 | exact | none | SERTA picks e4 = −0.04 (148.0 in IBTE) | strain ±0.002; μ ±5% |
| best_inplane_2pct | e6 = ±0.04, 205.1 (203.5 / 206.7) | −0.008 | +0.008 (expand e6) | e1 = −0.02 at 158.3 (diagonal only); e6 = ±0.02 (label read as Voigt) | strain ±0.002; μ ±8% |

The 2% μ_zz ranking in IBTE is 163.9 (e3−), 148.0 (e4−), 144.1 (e4+) and 143.2 (e5−). In SERTA it is 71.0 (e4−), 68.6 (e4+), 64.7 (e5−) and 64.1 (e3−).

Mirror-pair noise in the archive (compress vs expand of e4, e5 and e6) is at most 1.8% on the diagonals at 3%, and at most 3.2% at 1–2%, except μ_zz at 1% e4 (6.0%). On the strained off-diagonal element it is 0.1–7.7% at 2–3%, and 15–32% at 1% e4/e5, where that element is only about 6–8 cm²/Vs.

## Oracle

The oracle makes 18 evaluations for 44 of 50 units (88%):

- 5 band purchases: unstrained, compress 1/2/3% e3 and expand 3% e3;
- 13 IBTE purchases:
  - unstrained;
  - compress 1 and 2% e3;
  - expand 1 and 2% e3;
  - expand 1 and 2% e2;
  - compress_3%_e4;
  - expand_3%_e6;
  - compress_2% e4, e5, e6 and e1.

It uses three estimators:

- saturating (geometric) extrapolation of the 0/1/2% series for expand_3%_e3, expand_3%_e2 and the non-graded elements of compress_3%_e3;
- the bracket midpoint for compress_3%_e3 μ_zz;
- mirror partners for the shear targets.

It computes the strains from the cells. The estimated agent length is 40–70 tool calls.

Results:

- **Gate (2026-10-19, cpuset 6,7,22,23):** passed. All 9 ablations fail, 7 robustness variants pass, and rounding to 4 significant digits passes. There are 8 decision points and 8 rubric items. The task is 3.6 MB.
- **Harbor oracle (`runs/build-gan-hole-piezomobility/oracle-r1-gan-hole-piezomobility`):** reward 1.0, runtime 47 s.

## Open risks

- **Dev-easy risk.** The c-axis crossing and its μ_zz gain are in the literature (crystal-field reversal in strained GaN), so dp4 and the μ_zz categorical may be recalled. Strong models probably also know mirror symmetry (dp1). Difficulty should come from four forks:
  - dp2: whether the agent validates the x↔y swap;
  - dp5: a SERTA survey under the budget;
  - dp6 and dp7: principal axes, and the factor 2 in the shear labels;
  - the dp1 temptation to extrapolate instead of using the partner.
- **Wide window.** The compress_3%_e3 window is wide by necessity: it checks the regime, not the value.
- **Inferred physics.** The e1/e2 inequivalence at 2–3% (up to 21% at expand 3%) is taken as physical. It is consistent with the GW splittings (44.5 vs 53.6 meV at expand 3%), but it could not be checked independently.
- **Thin margins.**
  - The W4 same-sign-series ablation fails μ_zz by only 3.7 points beyond the tolerance; μ_yy fails by 8.7.
  - The compress_3%_e6 linear-series robustness variant is 2.6 points inside the μ_xx tolerance.
- **Band character.** Band character is ambiguous on the coarse grid near the crossing (1% compression). The rubric only asks for "between 1% and 2%".
- **Size of results.** Each band purchase returns about 95 KB of JSON. The instruction asks the agent to redirect large results to a file.
