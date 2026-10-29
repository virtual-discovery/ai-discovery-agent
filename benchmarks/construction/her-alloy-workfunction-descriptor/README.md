# Construction of her-alloy-workfunction-descriptor

Released-start task built on Liasi et al., arXiv:2608.28347 (2026), Zenodo 22031607 (CC BY 4.0).

## Rebuild

1. `./fetch.sh` downloads the four tarballs (about 765 MB; md5 checked) into `~/.cache/ai-theorist/her-alloy-workfunction-descriptor/` and unpacks what is needed into `x/`.
2. `python3 -B extract.py` writes `pool.json.gz` (1,025 measurement areas with EDX composition, SECCM currents, surrogate histograms and work functions; DFT records; coverage paths).
3. `python3 -B make_task.py` writes the catalog, `descriptors.json`, the lab data, the private targets and the ablation values, and prints `calibration.json`.
4. `.venv/bin/python -m benchmarks.tools.check_task her-alloy-workfunction-descriptor --cpuset ...` writes `check.json`.

`model.py` holds the activity-model utilities (numpy/scipy; a distribution-averaged Sabatier rate with a Koutecky-Levich cap, library-specific k and j_lim, log10 Huber residuals).

## Validation against the archive

- The histogrammed descriptors (0.005 eV bins) reproduce the authors' production multi-potential fits. For Lib-1 model A with a fitted centre: j0 2.93, alpha_c 0.416, centre -0.460 eV, j_lim 1.947, R2_log 0.9759. `04-activity/runs/summary-all.csv` gives 2.92 / 0.416 / -0.460 / 1.947 / 0.976. Lib-2 and Lib-3 match to the same precision.
- The coverage omega reproduces `repulsion.json` (equimolar adsorption branch 0.2483 eV/ML) with a plain linear fit over all nine points. Note that the authors' code fits all nine points, although their README says 0.1–0.9 ML.

## Design and deviations

See `design.json` (`notes`, `decision_points`, `calibration`). In short:

- **Best composition not graded.** The measured Lib-3 optimum is broad and noisy. It is rubric item 8.
- **alpha_wf graded as a numeric window.** The conditional work-function coefficient is graded, not a significance verdict.
- **Surrogate output precomputed.** The surrogate inference is served as precomputed per-spot descriptors, so the lab needs no torch.
- **Coverage-shift target added.** It is a DFT route to the volcano centre that is independent of the activity fit.
- **Lib-3 calibration capped at 5 spots.** With 10–20 spots, a linear composition regression inside Lib-3 reaches RMSE about 0.10 without any physics.

## Calibration (Lib-3 RMSE in log10 units; tolerance 0.135)

| analysis | value |
|---|---|
| oracle (model A, coverage-corrected distribution, cap, library-specific k/j_lim, 5 FPS calibration spots) | 0.097 |
| combined model C / centre fixed at -0.461 eV / alternative FPS calibration design | 0.097 / 0.101 / 0.094 |
| model A with 5 random calibration spots: median, p90 (7% exceed 0.135) | 0.101, 0.128 |
| model C with 5 random calibration spots: median, p90 | 0.101, 0.114 |
| wrong: mean-energy descriptor / k-only calibration / no cap / WF-only / centre -0.24 / no calibration | 0.144 / 0.145 / 0.153 / 0.160 / 0.211 / 0.247 |
| wrong but failing elsewhere: linear-combination descriptor (centre -0.273) / centre from total isotherm slope (shift 0.70) | 0.125 / 0.121 |

Windows for the other graded fields:

- **Volcano centre:** [-0.52, -0.36]. Reference -0.443; nearest wrong values -0.273 (linear combination), -0.24 (literature) and -0.81 (mean energy).
- **Coverage shift:** [0.07, 0.45]. Reference 0.221; the 0.1–0.9 ML window gives 0.186 and the fit through the dilute point 0.274. Nearest wrong values 0.70 (total slope) and 0.77 (first-last differential).
- **alpha_wf:** [-0.06, 0.13]. Reference 0.065; row subsets give 0.064–0.100. Nearest wrong values 0.157 with dilute-centre binding (0.146 at minimum over row subsets), 0.178 for WF only, 0.279 for one shared prefactor.

## Oracle

The oracle buys:

- rows 3, 7, 12 and 16 of Lib-1 and Lib-2 at all four potentials (32 units);
- 5 Lib-3 spots chosen by farthest-point sampling in composition (10 units);
- the 10 equimolar coverage paths (10 units);
- pure Pd and Pt from DFT (2 units).

That is 54 of 60 units in 49 evaluations, and it runs in about 10 s.

- **Gate (2026-10-26):** passed. All 11 ablations fail and all 6 robustness variations pass.
- **Harbor oracle:** reward 1.0.

## Open risks

- **Lib-3 calibration-spot choice.** About 7% of random 5-spot choices with model A exceed the RMSE tolerance. The choice is a real design decision: spots must constrain both k and j_lim.
- **alpha_wf margin.** The upper margin is modest. The dilute-binding wrong path sits 0.016–0.03 above the window, depending on which rows were bought.
- **Variation between decorations.** The single-decoration coverage slopes span 0.01–0.65 eV/ML. Agents who buy three or fewer decorations may leave the shift window.
