# Construction of hzo-vacancy-local-environment

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). The source is Taper et al., PRB 113, 134114 (2026); data DOI 10.17863/CAM.121316, CC BY 4.0. Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. `python3 benchmarks/construction/hzo-vacancy-local-environment/fetch.py` downloads `Data.zip` (12.8 MB, md5 70721b2c…) and `Readme.md` to `~/.cache/ai-theorist/hzo-vacancy-local-environment/`. It verifies md5 and sha256 and extracts the nine members used.
2. `python3 benchmarks/construction/hzo-vacancy-local-environment/extract.py` writes `pool.json` (134 KB) and asserts every archive relation used.
3. `cd benchmarks/construction/hzo-vacancy-local-environment && python3 -B make_task.py` writes the catalog, structures, x = 0.50 table, lab data, private targets and ablation and robustness values, and prints the calibration table.
4. `python3 -B calibrate.py` prints the pass rates of the correct method and of the wrong paths under other purchase plans (400 designs each).
5. `.venv/bin/python -m benchmarks.tools.check_task hzo-vacancy-local-environment --cpuset 8,9,24,25` writes `check.json`.

`model.py` holds the environment classification (2.3 Å first shell), the model fits and the random-alloy population code shared by steps 3–4. `solution/solve.py` re-implements the same oracle logic stand-alone; the gate run reproduces `oracle_expected` in `ablation_values.json` exactly.

## What the task tests

The released data show four things:

- O3 vacancies cost about 0.28 eV more than O4 vacancies at equal Zr count.
- Each Zr first-shell neighbour lowers E_f by 0.08–0.10 eV, approximately additively.
- A composition shift of about −0.03 eV per 0.25 in x appears at fixed environment. It is below the 0.02–0.035 eV within-class scatter, and it is confounded with the Zr slope in pooled fits because each cell's environment mix follows its composition.
- The CASTEP cells use their own decorations, per functional.

The graded quantities are:

- held-out site energies at x = 0.25 and 0.75;
- the random-alloy vacancy population at 700 K (the n_Zr distribution and the site fraction);
- the common-slope LDA description.

Decision points (`design.json`):

| DP | Kind | Wrong path |
|---|---|---|
| dp1 | sublattice-resolved cation model | pooled n_Zr regression; geometry descriptors |
| dp2 | discovered confounder | composition-blind pooled model |
| dp3 | discriminating purchases | transfer of the free x = 0.50 model |
| dp4 | regime judgement | 96-atom cell statistics; fixed-composition sampling |
| dp5 | degeneracy × Boltzmann | no degeneracy; lowest class only |
| dp6 | cascade | unpurchased rare classes dropped (the x = 0.25 4-Zr site is withheld) |
| dp7 | metadata inconsistency | CASTEP sites given the VASP cell's environments |
| dp8 | functional transfer | PBE description reported for LDA |

dp2, dp3, dp4, dp6 and dp7 rest on evidence that has to be found in the data or bought. dp1 and dp5 are standard.

## Build notes and deviations from the design

- **Chemical potential.** The design wanted a composition-dependent O-poor limit (oxide equilibrium). The archive has no O2, metal or oxide reference energies, and the authors used μ_O = −10.09 / −7.33 / −4.57 eV for every composition. The prompt therefore states μ_O = −10.09 eV, and that fork was dropped.
- **Graded targets changed:**
  - **Per-site energies:** the design graded each site within 0.12 eV, plus Spearman ≥ 0.8. A 0.12 eV per-site tolerance cannot see the composition effect (about 0.03 eV), and Spearman over sites spanning two sublattices is trivial. The checks are now per-composition RMSE ≤ 0.035 eV and |mean error| ≤ 0.02 eV over 10 + 10 withheld sites.
  - **"Lowest-E_f class at x = 0.75":** dropped; it is 4-Zr at every composition.
  - **O3 fraction at x = 0.50:** dropped. It is below 1 %, depends exponentially on the 0.28 eV offset, and does not separate the degeneracy paths (no-degeneracy gives 0.0062 against 0.0084).
  - **Replacements:** the n_Zr distribution of vacancies, which tests random-alloy statistics, degeneracy and rare classes; the site fraction, whose x = 0.75 value is sensitive to the composition shift; and the LDA description.
- **x = 0.50 data.** The design's "seed set of about 10 x = 0.50 energies" became the authors' full x = 0.50 table (FIGS12), staged free as released. Its `E_vac` is the O-rich formation energy, which the release does not state. One vacancy plus pristine purchase pins it (μ = −4.57 eV), and the oracle does this.
- **CASTEP data.** The design assumed "CASTEP energies for the reference x = 0.50 cell", but the README of the archive is wrong on this point. The released FIGS3 labels and coordinates show that every functional's calculation used its own decoration: 7–24 of 64 labels match any staged cell by index. The labels, not the staged cells' environments, fit the energies (residuals 0.02–0.04 eV against 0.10–0.16 eV). The labels are returned only with a purchase, which makes this inconsistency the dp7 trap. The reference cell `Fig2a_Fig3_HZO.cif` is not staged, because it is CompMat-Bench's cell.
- **Rubric direction.** The design's rubric said O3 vacancies are cheaper ("fewer bonds broken"). The data say the opposite: O3 vacancies cost about 0.28 eV more, with shorter and stronger bonds. The rubric follows the data.
- **Geometry model.** The released geometric descriptors (FIGS12) are computed on the ideal cells, where they take one value per sublattice, and the archive has no relaxed geometries. The geometry explanation is therefore a dp1 wrong path and a rubric item, not a separate graded target.
- **Charged vacancies excluded.** The `Bulk VBM` column is, according to CompMat-Bench task 028's notes, the pristine Fermi energy and not the true VBM (0.92 eV apart for the reference cell). Charged-defect thermodynamics is also CompMat-Bench's task 028. The +2 energies also show no monotonic Zr-count trend: class means lie within about 0.1 eV of each other, with within-class spreads of 0.01–0.08 eV.
- **Dropped decision point.** "Underdetermined per-composition fits" was dropped: with 9 purchased sites per composition, the independent fit passes (it is now a robustness variant).
- **CompMat-Bench separation.** See `SOURCE.md`. No graded value involves the Fig2a_Fig3 cell, any of the charged, HSE or COHP data, or any CompMat-Bench target. Nothing from CompMat-Bench is staged.
- **Leak fix.** FIGS3 rows are sorted by enthalpy in the archive, so the catalog lists the CASTEP sites by index; the archive's order would leak the energy ranking.

## Calibration (full table in `design.json` → `calibration`)

| Check | Oracle | Robustness range | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| mean error x=0.25 / 0.75 (eV) | +0.005 / −0.000 | ≤ 0.016 | −0.026 / +0.030 (x0.50 transfer) | 0.02 |
| RMSE x=0.25 / 0.75 (eV) | 0.022 / 0.025 | ≤ 0.032 | 0.034 / 0.040 (transfer), 0.125 (pooled) | 0.035 |
| mean n_Zr x=0.25 (ref 2.239) | 2.255 | 2.20–2.36 | 2.785 (pooled), 3.071 (cell statistics) | ±0.25 |
| 4-Zr share x=0.25 (ref 0.103) | 0.101 | 0.092–0.121 | 0.056 (fixed-composition cells), 0 (dropped class) | ±0.03 |
| mean n_Zr / 4-Zr share x=0.75 (ref 3.693 / 0.727) | 3.678 / 0.715 | 3.66–3.75 / 0.68–0.77 | 4.0 / 1.0 (lowest class) | ±0.15 / ±0.08 |
| log10 site fraction x=0.25 (ref −5.866) | −0.060 | −0.098 to −0.034 | +0.212 (cell statistics) | ±0.20 |
| log10 site fraction x=0.75 (ref −4.402) | −0.018 | −0.049 to +0.041 | −0.190 (composition-blind) | ±0.15 |
| LDA offset Δ (ref 0.120 eV) | 0.127 | 0.164 (other 12 sites) | 0.226 (CASTEP-PBE) | ±0.075 |
| LDA slope b (ref −0.111 eV) | −0.121 | −0.106 | +0.045 (VASP-cell environments) | ±0.05 |

Pass rates from `calibrate.py`:

| Path | Purchase plan | All checks pass |
|---|---|---|
| Correct method | class-targeted 9 + 9 | 0.975 |
| Correct method | random 9 + 9 | 0.983 |
| Correct method | random 6 + 6 | 0.965 |
| Correct method | random 12 + 12 | 0.970 |
| Composition-blind model | any of these | 0.000 |

All failures of the correct method come from noise in the LDA estimate: the LDA check passes 0.93, 0.977 and 0.992 of the time with 8, 10 and 12 random LDA sites. With the VASP cell's environments it passes 0.03–0.05.

Rounding every submitted number to 4 significant digits passes.

## Oracle and length

- **Oracle:** 32 evaluations, 32 of 36 units (89 %):
  - 3 pristine cells;
  - 1 x = 0.50 vacancy to pin the table's μ_O;
  - 9 + 9 class-targeted vacancies at x = 0.25 and 0.75;
  - 10 CASTEP-LDA vacancies.
- **Runtime:** 3 s in the gate container; the Harbor oracle took 47 s.
- **Catalog:** 713 units, so the budget covers 5 % of it.
- **Expected agent length:** 45–75 tool calls (structure parsing, table conversion, about 30 purchases, model fitting, population code).

## Open risks

- **Dev-easy risk.** A careful agent that fits per-composition offsets, uses binomial statistics and reads the CASTEP labels passes. The traps are a small composition shift, the cells' unrepresentative rare environments, an unpurchasable rare class and per-functional decorations; strong models may find all of them.
- **Population reference.** The population truth averages Boltzmann factors over each class's archived sites (within-class scatter included). Model-based estimates without a scatter correction sit 0.03–0.10 low in log10 at x = 0.25, inside the ±0.20 tolerance.
- **Paper figures.** The paper plots per-site E_f for the three cells (figure level, readable to about 0.02 eV). The populations and the LDA fit are not in the paper.
- **CASTEP labels.** The FIGS3 labels are trusted because they are internally consistent with the energies. The archive README's claim that they belong to `Fig2a_Fig3_HZO.cif` is contradicted by the data.
- **Fixed lattice.** The composition shift cannot be split into chemistry and fixed-volume strain with these data. It is a rubric point only.

## Disk

- **Peak:** about 130 MB in `~/.cache/ai-theorist/hzo-vacancy-local-environment/`, made of `Data.zip` and a full unzip during exploration (113 MB). The raw downloads were deleted after the gate and the Harbor run.
- **Kept:** `pool.json` (134 KB) and the refetch scripts with checksums.
- **Task package:** 0.26 MB.

## Gate

The gate passed on 2026-10-19. `check.json` reports 8 decision points, all 12 ablations failing, all 7 robustness variants passing, rounding passing and 40 hidden values with no leaks. The Harbor oracle run scored reward 1.0 (`runs/build-hzo-vacancy-local-environment/oracle-r1-hzo-vacancy-local-environment`).
