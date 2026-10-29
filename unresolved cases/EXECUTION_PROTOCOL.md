# Execution protocol for an agent taking one of these cards

## 1. Pick-up

1. Read, in this order: `README.md` (this folder), `LI_ANALYSIS.md`, `SHARED_RUNBOOK.md`, then ONLY the card folder you take (`goal.md`,
   `acceptance_draft.json`). The cards' "Read this and nothing more" sections are binding.
2. Open a branch `case/<nn>-<slug>` from `main`. Work only inside `unresolved cases/<nn>_<slug>/`. Commit the frozen card first.
3. Confirm the environment (SHARED_RUNBOOK §2). If a required S0 input is marked `TO_FILL` and you cannot obtain it (paywalled paper),
   stop at S0, write `BLOCKED_ON_INPUT.md` naming the exact numbers needed, commit, and report. Do not guess them.
4. Budget: state the planned GPU hours and $ in `results/budget.json` before renting anything; the cap is 60 % of the remaining credit;
   speed over cost within the cap (round-8 instruction). Rented resources are destroyed before the turn ends (RUNBOOK §6).

## 2. Running

- Stages are those in the card's §8 table; each writes a state file and is resumable. Run the LITCHECK stage (G0) before anything heavier
  and stop if it fails — that is a legitimate terminal outcome (`LITCHECK_FAIL`), not a problem to engineer around.
- Write the scorer BEFORE production and validate it on a synthetic or literature input so it can both pass and fail.
- Do not add stages, k-points, parameters, or "quick checks" beyond the card. If something looks necessary, finish the card as written, then
  propose a successor card in `RESULT.md` §"Next".

## 3. Reporting

- `RESULT.md` opens with a verdict box: status string, one sentence, the gate table. Then figures with captions, provenance, compute log,
  limitations, "Next".
- Every number in the verdict traces to a file in `results/` and a stage directory; cite the file next to the number.
- Report failures as failures. A gate that could not be scored is "not scored (reason)", never "pass".
- Keep raw heavy data out of git; list them in `LARGE_FILES.tsv` with size, sha256 and where they are.

## 4. Coverage audit (mandatory before the final commit)

Enumerate the card's stage table and gate table and check each against files on disk and keys in `results/*.json`. State coverage as
"N of M planned stages ran; gates scored: …; gates not scored: … (reason)". Missing work hides as confident prose assembled from neighbouring
numbers — this exact failure cost us a week once.

## 5. Hand-back

Push the branch, open a pull request titled `case <nn>: <STATUS> — <one line>`, and paste the verdict box in the PR body. Do not merge.
