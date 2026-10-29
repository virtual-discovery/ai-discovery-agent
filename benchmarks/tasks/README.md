# AI Theorist physics benchmark: 60 long tasks

Status: **built and frozen, 2026-10-27.** 60 tasks are frozen in `evidence/release.json`: 50 released starts and 10 fresh starts. Every task passes `check_task` and scores reward 1.0 with the Harbor oracle.

Designs came from a 61-candidate literature survey (minus Na3RuO4), kept in git history. The build spec is `benchmarks/tools/BUILDING.md`. Each task's construction scripts, calibration table and build notes are in `benchmarks/construction/<task>/README.md`.

## What every task guarantees

- **Labs and data:**
  - The lab is a budgeted service with private data.
  - The agent image has no Quantum ESPRESSO, so every calculation goes through the budget.
  - Released-start tasks serve author-archived data. Many also run a small numpy evaluator that is validated against the archive.
  - Fresh-start tasks run new pw.x/ph.x calculations whose cost scales with wall time.
- **Decision points:** at least 6 per task, median 8. Each comes with ablations computed from the real data, and every ablation fails grading.
- **Tolerances:** set between the spread of defensible methods and the nearest wrong path. Each task records this in `design.json` → `calibration`. Defensible variants and 4-significant-digit rounding must pass.
- **No leaks:** the graded hidden values are scanned and absent from everything the agent can read.
- **No recall:** graded values must be recoverable without the paper's numbers, because famous headline values are avoided or graded at protocol-specific settings. In MnWO4 this rule removed the fitted-exchange targets.
- **Archive defects:** many decision points rest on verified defects in the archives, for example:
  - mislabelled units, runs and folders;
  - duplicated or corrupted files;
  - wrong signs in published formulas;
  - unconverged or unrelaxed runs left in place.

## Development pilots (Sonnet 5 and GPT-5.6 Terra, one run each)

| | Sonnet 5 | GPT-5.6 Terra |
|---|---|---|
| Tasks passed | 13/60 | 12/60 |
| Mean fraction of checks passed | 0.71 | 0.75 |
| Median API calls | 58 | 25 |
| Median wall time | 23 min | 5 min |
| Median GPT-6-astra-equivalent cost | $10.0 | $2.1 |

- **Outcomes by task:** both models pass 6 tasks (dev-easy), one model passes 13, and neither passes the other 41.
  - The dev-easy tasks are Fe exchange, CuO exchange, hydrogen plasma-facing, InP, Cr2O3 and kagome 135.
- **By start type:** released starts pass 19 of 100 runs, fresh starts 6 of 20.
- **Pilot spend:** $194 in total. The ledger is `benchmarks/construction/pilots.json`.
- **Why runs fail:** failures typically miss one to four checks. Each was traced to a designed decision point, not a task bug.
- **Shared answers investigated:** three times both models converged on the same wrong answer, and each case was diagnosed before freezing.
  - **Cu2WSe4:** the orbital-moment operator was ambiguous. The wording was fixed and the models' path added as an ablation.
  - **MnTe magnons:** two genuine stacked errors. Added as an ablation, and the slices are now sold in tiles.
  - **HgTe:** a ridge misassignment that the raw data rule out. Added as an ablation.

## Harness settings

- **Limits:** $20 cap, 200 calls, 32k output tokens, up to 3 retries of a truncated turn, and transient 52x errors retried (harness 0.8.0).
- **Time limit:** 3 h per agent episode. Three fresh-start tasks get 4 h because their oracles need 60–100 min of QE: `fresh-li-bcc-soft-shear`, `fresh-gete-ferroelectric-well` and `fresh-srtio3-soft-modes`.
- **QE tasks need the CPU-pinning daemon:** run `python -m benchmarks.cpu_pinning`.

To run a model on the frozen set:

```sh
source .venv/bin/activate
set -a; . ~/.env; set +a
python -m benchmarks.run --policy baselines --suite all \
    --campaign <name> --concurrency 2 --wait --execute
```

## Open review items

Physical-reasoning rubrics still need physicist review, as they did in v1. These choices should be reviewed before results are published:

- **Truth that departs from the paper because of a verified correction:**
  - **skutterudite:** corrupted isotope rates are repaired.
  - **1D wires:** the paper's C_V formula has a sign error.
  - **ZnO:** the authors' band filling is corrected.
  - **half-Heusler:** the hinge termination is set by reversed hopping vectors in the authors' code.
- **Truth computed by our validated evaluators or re-implementations rather than archived numbers:**
  - Cu2WSe4, TaN (Monte Carlo four-phonon rates), MPEA (barrier model), ZrSiX (Chambers solver), Pd–H, kagome 135, Ca3Ti2O7, Ge polymorphs.
  - The fresh-start tasks, whose truth comes from our own QE runs.
- **Synthetic data:** `raman-polymorph-orientation` uses synthetic measurements from a documented forward model.
- **Hand-traced references:** `tmd-moire-chern-twist` (C3 labels from a figure) and the HgTe SdH ridges.
- **Not peer reviewed:** the sources of `hydride-coulomb-mustar` and `mg2x-midir-forbidden-transition`.
- **Thin margins and other flags:** recorded per task in each construction README under "Risks".
