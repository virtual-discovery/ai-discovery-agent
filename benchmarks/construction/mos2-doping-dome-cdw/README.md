# Construction of mos2-doping-dome-cdw

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** Girotto Erhardt, Berges, Poncé and Novko, npj 2D Mater. Appl. 9, 44 (2025), with Materials Cloud record 51yyy-ksf72 (DOI 10.24435/materialscloud:z1-aw, CC BY 4.0).
- **Provenance, checksums, transformations and what is hidden:** the task's `SOURCE.md`.

## Rebuild

Run compute in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) ...`), from this folder:

1. `bash fetch.sh [DEST]` downloads the six archives (11 MB) and `README.md` through the per-file API links, unpacks them and verifies the 93 members in `members.md5`. The default DEST is `~/.cache/ai-theorist/mos2-doping-dome-cdw`. `fig6.tar.gz` (model III) and `figS11.tar.gz` are not fetched.
2. `python3 -B extract.py [DEST]` writes `pool.json` (1.0 MB): every number used, unmodified.
3. `python3 -B validate_solver.py` (about 10 min) compares `eliashberg.py` with the 318 archived Tc values in `fig2e/tc_*.dat` and writes `solver_validation.json`. The maximum deviation is 0.21 % (one Einstein-model value at μ* = 0); the full-DOS and constant-DOS α²F values agree to ≤ 0.03 %.
4. `python3 -B make_task.py` (about 4 min) writes:
   - the catalog;
   - `lab/data.json` (0.85 MB) and `lab/eliashberg.py` (copied from here);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `reference.json` here.

   It asserts that the solver reproduces the archived Tc at the hidden dopings and that the high-doping CDW spectra map uniquely onto the authors' table.
5. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task mos2-doping-dome-cdw --cpuset <set>` writes `check.json`.
6. `python3 -B record_oracle.py` copies the gate's oracle run into `oracle_run.json`. `python3 -B write_design.py` then writes `design.json` with the calibration table.

## What the data showed (design inputs)

- **The solver can be validated exactly.** The authors archived Tc tables for six μ* values and seven methods, computed with their GPL solver ebmb 2.0. I implemented the paper's Eqs. (1)–(9) in numpy, using the numerical conventions of `tc.py`: cutoff 20 ω₂, μ* unscaled at ω₂ and rescaled at the Matsubara cutoff, and the spectrum rounded to 6 decimals. It reproduces the full-DOS (DA) and constant-DOS (CA) values to ≤ 0.03 %, so the lab runs a calculation that is equivalent to the authors'. ebmb was read for these conventions; no code was copied.
- **The archived fig5 Tc tables are Allen–Dynes values.**
  - The FA column of the 1×1 table equals Allen–Dynes (with strong-coupling and shape corrections) of the FA spectra to 0.1 %.
  - The high-doping CDW rows equal Allen–Dynes of the fig4 spectra to 0.2 %.
  - The A and NA columns of the 1×1 table do not match the archived A/NA spectra (computed Allen–Dynes 11.98 versus archived 5.07 K at 0.14, and 33.99 versus 29.65 K at 0.17). Their origin is unknown, so they are served as archived and flagged in the rubric only.
- **The high-doping CDW spectra carry no doping.** The files are `2x2cdw_dop1–3`. Allen–Dynes Tc of their FA (and A, NA) spectra matches the fig5 rows at 1.4, 1.5 and 1.6 electrons per 2×2 cell (0.35, 0.375 and 0.40 per 1×1 cell). The paper's Fig. 4 caption gives n1 = 3.0 × 10¹⁴ cm⁻²; 1.4 per 2×2 cell is 3.98 × 10¹⁴ cm⁻². n2 and n3 (4.3, 4.6) agree with the match.
- **λ jumps across the Q-valley/soft-mode regime.** λ_FA is 0.77 at 0.14, 1.14 at 0.15 and 1.41 at 0.16, and Tc rises steeply. Interpolating Tc linearly across 0.15 fails by 10.7 %. The authors' λ/ω_log table covers the withheld dopings.
- **The low-temperature chemical potential sits 43–66 meV above the DOS energy zero.** The zero is the Fermi level at the 800 K smearing. At Tc the chemical potential lies inside the Q-valley edge, so constant-DOS Eliashberg is 6–18 % low for 0.13–0.19.
- **Three different "critical dopings" in fig3:**
  - the PES energy crossing (E of the CDW endpoint minus E of 1×1 H) at 0.1890;
  - the 1×1 Γ* soft mode turning imaginary at 0.2029;
  - the lower edge of the CDW branch (hysteresis) at 0.1850.

  The PES file suffix is the 2×2-cell tot_charge (`pes_074` is 0.74 per 2×2 cell, 0.185 per 1×1 cell, 2.105 × 10¹⁴ cm⁻²). The fig3 axes are in 10¹⁴ cm⁻²; 1 electron per 1×1 cell (a = 3.1856 Å) is 11.3787 × 10¹⁴ cm⁻², which reproduces those axes exactly. The PES curves still descend slightly at the CDW endpoint; a cubic-fit well minimum moves x_c by less than 10⁻⁵.
- **SOC.** Spectra and DOS with SOC exist only at 0.18 and 0.19, from one Wannier Hamiltonian. SOC lowers Tc by 25–30 %.

## Deviations from the design (candidates.json / DETAILS.md section 2)

- **Graded fields changed.**
  - Design target 1 (Tc at 0.15 and 0.17, ±10 %) is kept with ±3 %. The defensible routes spread from −1.3 % to +0.3 %. The nearest wrong path is Allen–Dynes at 0.17 (−3.75 %); ±10 % would also have let constant-DOS and McMillan pass.
  - Design target 2 (2×2 CDW Tc at 0.77 and 1.4 per 2×2 cell) is replaced by Tc of the high-doping CDW at x = 0.35 and 0.40 (dop1 and dop3), ±4 %. No spectrum or DOS is archived for the low-doping CDW (0.74–0.77), whose Tc exists only as Allen–Dynes numbers. 1.4 per 2×2 cell is graded in its 1×1 form (0.35) and tests the label inference.
  - Design target 3 (dome maximum and its doping) is not graded. In the coexistence window the CDW Tc are archived Allen–Dynes values without spectra, slightly above the 1×1 values, so the maximum is not uniquely defined by the data. The paper's headline 34.6 K (1×1 at 0.20, Allen–Dynes) is within 0.5 % of the scalar-relativistic Tc at x_c and would pass by recall. It is replaced by the highest Tc of the 1×1 ground-state branch with SOC (27.04 K, ±3 %), which is not in the paper (the paper quotes 21.6 K with SOC at 0.18, Allen–Dynes).
  - Design target 4 (x_c from the PES) is kept, with ±0.002 instead of ±0.005. ±0.005 would accept the CDW-branch edge at 0.185.
- **No Wannier Hamiltonians in the lab.** The archived DOS files are served at the released dopings, and the `eliashberg_custom` request allows rigid-band use of a released DOS with any electron count. The Wannier files are 3.6 MB per doping, and the DOS at the hidden dopings would make the withheld Tc directly computable.
- **Decision points.** The design's forks on units, the adiabatic divergence and flavour, the DOS treatment, formula vs Eliashberg and phase selection are kept (here dp6, dp2, dp3, dp4 and dp5). The design's dp6 (tuning μ* to experiment) is a rubric item only, because μ* = 0.13 is stated. I added three evidence-based forks:
  - the λ jump at the withheld doping (dp1);
  - the unlabelled high-doping spectra and the inconsistent caption (dp8);
  - the archived tables being Allen–Dynes (dp9).
- **Not staged:** phonon dispersions, phonon spectral function, band structures, the substrate data (figS11) and the elphmod model (fig6).

## Decision points (9)

| dp | Fork | Kind |
|---|---|---|
| dp1 | Hidden-doping Tc across the λ jump: neighbour interpolation fails at 0.15; the tabulated λ/ω_log or a rescaled spectrum with a rigid-band DOS is needed | regime found in the data; discriminating purchase |
| dp2 | A / NA / FA spectra (the paper's Fig. 2d is NA) | convention (stated) |
| dp3 | Constant vs energy-dependent DOS | model class |
| dp4 | Allen–Dynes / McMillan vs Eliashberg at λ > 1.5 | model class |
| dp5 | PES energy crossing vs soft-mode instability vs CDW-branch edge; cascades into the SOC maximum | regime judgement, coexistence and hysteresis in the data |
| dp6 | 2×2-cell charges and 10¹⁴ cm⁻² axes vs electrons per 1×1 cell | convention found in metadata |
| dp7 | SOC for the ground-state maximum | data coverage |
| dp8 | Dopings of `2x2cdw_dop1–3`: the caption (n1 = 3.0 × 10¹⁴) and table order are wrong; Allen–Dynes matching against the archived table resolves them | label inconsistency found by a cross-check |
| dp9 | The archived Tc tables are Allen–Dynes (header says only "Tc"); their 1×1 A/NA columns do not match the spectra | metadata inconsistency found by a cross-check |

dp1, dp5, dp8 and dp9 rest on evidence the agent has to find and pay for. dp5 cascades into the SOC field, and dp8 into the CDW field. dp2–dp4 are stated protocol, but their wrong paths remain realistic.

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Defensible spread | Nearest wrong paths | Tolerance |
|---|---|---|---|---|
| Tc 1×1, 0.15 (K) | 15.991 (archived DA) | −1.3 % to −1.0 % (moments × Eliashberg/AD ratio; rescaled two-neighbour spectrum, either rigid-band DOS; λ-weighted mix) | AD −1.8 %, CA −8.3 %, unscaled average spectrum −8.1 %, linear Tc interpolation −10.7 %, archived table −12.3 %, NA +26 % | ±3 % |
| Tc 1×1, 0.17 (K) | 26.802 (archived DA) | −1.1 % to +0.3 % | AD −3.75 %, archived table −5.1 %, CA −6.2 %, McMillan −16 %, NA +21 % | ±3 % |
| x_c (e per 1×1 cell) | 0.18895 | +0.0007 (quadratic crossing), 0.0000 (well minimum) | CDW-branch edge −0.0040, soft mode +0.0139, 2×2 units 0.756, 10¹⁴ cm⁻² 2.15 | ±0.002 |
| Tc 1×1 with SOC at x_c (K) | 27.042 | +1.4 % (quadratic x_c), +0.1 % (ratio method) | AD −4.0 %, CA −6.2 %, CDW edge −7.9 %, no SOC +27.6 %, soft mode +27.8 % | ±3 % |
| Tc high-doping CDW, 0.35 / 0.40 (K) | 12.930 / 6.333 | direct lab output | McMillan −1.7 % / +6.1 %, AD or archived table +6.3 % / +10.7 %, caption mapping −23 % / +6.7 %, table order −31 % / −41 % | ±4 % |

Gray zone, not graded as a wrong path: one neighbour's spectral shape rescaled to λ(x) ignores the tabulated ω_log. It gives −0.2 % / −1.6 % at 0.15 and +2.9 % / −3.5 % at 0.17.

Rounding every submitted number to 4 significant digits passes.

## Oracle

The oracle makes 20 evaluations for 29 of 34 units (85 %):

- moments table: 1 unit;
- Eliashberg (full DOS, FA) at 0.14, 0.16 and 0.18: 6 units;
- FA spectra at 0.14 and 0.16, plus one `eliashberg_custom` cross-check at 0.15: 4 units;
- four PES: 4 units;
- both soft-mode branches and the CDW Mo displacement: 3 units;
- SOC Eliashberg at 0.18 and 0.19: 4 units;
- the archived CDW table: 1 unit;
- constant-DOS Eliashberg for dop1–3: 6 units.

Its graded errors are:

- Tc at 0.15 and 0.17: −1.09 % and −0.01 %;
- x_c: 0;
- SOC maximum: +0.001 %;
- CDW: 0.

Each lab Eliashberg run takes 0.1–8 s. The catalogue costs 58 units for retrievals alone, and every Eliashberg run costs 2 more. The expected agent length is 45–75 tool calls.

## Results

- **Gate (2026-10-19, cpuset 2,3,18,19):** passed.
  - 9 decision points; all 15 ablations fail.
  - 7 robustness variants and 4-significant-digit rounding pass.
  - 8 rubric items, 6 hidden values, no leaks.
  - Package 0.89 MB; the oracle took 27 s in the gate container.
- **Harbor oracle (`runs/build-mos2-doping-dome-cdw/oracle-r1-mos2-doping-dome-cdw`):** reward 1.0, solver charge 29, runtime 69 s.

## Open risks

- **Dev-easy risk.** A careful agent buys the moments table, transfers Eliashberg/Allen–Dynes ratios, buys the PES, divides by 4, uses the SOC runs and cross-checks the CDW labels against the table. The traps are:
  - the λ jump (linear interpolation −10.7 %);
  - the soft-mode criterion;
  - the caption, which places dop1 at 0.264;
  - Allen–Dynes, which misses by 3.75–4 % on two fields.
- **Thin margins.**
  - Allen–Dynes at 0.17 fails by 0.75 points beyond the 3 % tolerance, and on the SOC field by 0.97 points. The same path also fails the CDW field by 2.3 points.
  - x_c from the CDW-branch edge fails by 0.002 beyond the tolerance.
- **Gray-zone estimators.** Single-neighbour spectral shapes at 0.17 sit near ±3 % (above). An agent that validates its estimator on a released doping will reject them, but an agent that does not may pass or fail by luck.
- **Inferred, not archived.**
  - The dopings of `2x2cdw_dop1–3` come from the Allen–Dynes match (0.2 %), not from an archive label.
  - The PES path endpoints are read from the paper's Fig. 3(c) (first point 1×1 H, last point 2×2 CDW).
  - The SOC maximum assumes that the scalar-relativistic PES boundary also holds with SOC; the rubric asks for this caveat.
- **Electronic temperature.** All spectra and the boundary are at T_scf = 800 K. The paper notes that x_c moves between 1.5 and 2.4 × 10¹⁴ cm⁻² with T_scf. The graded x_c is the archive's value, not a physical prediction.
- **Excluded-study separation.** No `fig6` (elphmod model III) data are used, and no Schobert (2024) quantity is reproduced. The task shares only the research group and the Eliashberg methodology with that study.

## Disk

- **Peak:** 134 MB in `~/.cache/ai-theorist/mos2-doping-dome-cdw`: the 11 MB of archives and the unpacked tree, which includes the unused Wannier files and figS11. The scratch subfolder held 5 MB (the paper, and the ebmb and elphmod sources read for conventions).
- **After the gate:** both were deleted. `pool.json` (1.0 MB) is kept here, and `fetch.sh` refetches the rest.
