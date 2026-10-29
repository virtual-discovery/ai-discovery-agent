# AI Theorist

A benchmark of 60 long physics model-building tasks for AI agents, and the harness to run it.

**Status (2026-10-28):** 60 tasks built and frozen (50 released starts, 10 fresh starts). Development pilots are done; Opus 5.5 and GPT-6-astra baselines are next. The improvement goal is in [benchmarks/tasks/RSI_GOAL.md](benchmarks/tasks/RSI_GOAL.md).

## Benchmark

Each task asks an agent to build and test a physics model from a budgeted lab, then submit graded numbers and a report. Tasks are modelled on GeneBench-Pro: an attempt is meant to cost a few dollars for GPT-6-astra, typically 25–80 calls and 0.1–4 h.

- **Released starts (50).** The lab serves author-archived data from 49 open-licence studies (Materials Cloud, Zenodo). Many labs also run a small numpy evaluator validated against the archive: Wannier/Berry, Eliashberg, BTE, Chambers, spin waves, kMC, and ports of ML potentials.
- **Fresh starts (10).** The lab runs new Quantum ESPRESSO pw.x/ph.x calculations; the charge scales with wall time.
- **Budget binds.** The agent image has no Quantum ESPRESSO, so every calculation goes through the lab budget. Budgets cover roughly 15–60% of each catalogue.
- **Graded on hidden numbers.** A pure `grade()` checks the hidden numbers. Graded values cannot be recalled from the literature: they are protocol-specific or withheld archive values.
- **Explanations reviewed separately.** Each task's report is reviewed against its `rubric.md`.
- **Decision points.** Each task has 6 or more, median 8. Many rest on verified archive defects that must be discovered in the data, for example mislabelled units or folders, duplicated or corrupted files, sign errors in published formulas, and unrelaxed or unconverged runs. Each decision point has ablations computed from the real data, and every ablation must fail grading.
- **Build gate.** Every task passes `check_task` (reference accepted, ablations fail, robustness variants and 4-significant-digit rounding pass, no hidden-value leaks, deterministic grading) and scores 1.0 with the Harbor oracle.

**Development pilots** (one run per task and model):

| | Sonnet 5 | GPT-5.6 Terra |
|---|---|---|
| Tasks passed | 13/60 | 12/60 |
| Mean fraction of checks passed | 0.71 | 0.75 |
| Median API calls | 58 | 25 |
| Median GPT-6-astra-equivalent cost | $10.0 | $2.1 |

- **By task:** both models pass 6 tasks, one model passes 13, and neither passes the other 41.
- **Cost:** pilots cost $194 in total.
- **More detail:** per-task results, harness settings and open review items are in [benchmarks/tasks/README.md](benchmarks/tasks/README.md).

## Setup

You need:
- Linux;
- Python 3.12+ with `venv`;
- Docker Engine with Compose;
- about 8 cores, 16 GB RAM and 30 GB of free disk.

Authenticate to GitHub to clone this private repository.

```bash
git clone https://github.com/exon-research/AI-Theorist.git
cd AI-Theorist
bash scripts/setup.sh
source .venv/bin/activate
python -m benchmarks.verify          # 60 tasks verified (60 frozen)
```

For paid runs, create `.env` (`cp .env.example .env; chmod 600 .env`) and set `OPENAI_API_KEY`, `ANTHROPIC_API_KEY` and `GEMINI_API_KEY`. Keep credentials on the host, outside task containers.

## Running

Commands preview the workload unless `--execute` is given.

```bash
# Fresh-start (QE) tasks need dedicated cores; keep this running during campaigns.
python -m benchmarks.cpu_pinning &

# Reference solutions through Harbor.
python -m benchmarks.run --policy oracle --suite all --campaign oracle --execute

# Baselines: Opus 5.5 and GPT-6-astra (or one model: --policy opus55 / astra).
set -a; . ./.env; set +a
python -m benchmarks.run --policy baselines --suite all --campaign baselines \
    --concurrency 2 --wait --execute
python -m benchmarks.tools.pilot_report baselines       # length, cost, checks passed
```

Notes:
- **Repeats and selection:** use `--task <id>` (repeatable) to pick tasks and `--replicate N` for repeats.
- **Reruns:** terminal jobs are skipped on rerun.
- **Paid campaigns run one at a time:** a budget lock serializes them, and `--wait` queues behind a running one.
- **Budget ceiling:** spend is checked against the $2,000 ceiling, and `--other-spend-usd` accounts for costs outside the runner.
- **Output locations:** runs go in the ignored `runs/` directory.

## Files

- **`benchmarks/tasks/`:** the 60 Harbor task packages (instruction, lab, private data, reference solution, ablations, rubric, provenance).
- **`benchmarks/construction/<task>/`:** scripts that refetch the pinned archives and rebuild each task, plus its calibration table, gate result and build notes. `pilots.json` holds the development-pilot ledger.
- **`benchmarks/tools/`:**
  - `BUILDING.md` (build spec);
  - `labkit/` (shared lab service, client and verifier);
  - `new_task.py`, `check_task.py` (gate), `release.py` (verify/freeze) and `pilot_report.py`.
- **`benchmarks/api_agent.py`, `durable_shell.py`:** the provider harness and durable command execution.
- **`benchmarks/run.py`:** the campaign runner.
- **`benchmarks/cpu_pinning.py`:** dedicated cores for QE labs.
- **`evidence/release.json`:** frozen task checksums and gate summaries.
- **`benchmarks/costs.json`:** API pricing and recorded spend.
