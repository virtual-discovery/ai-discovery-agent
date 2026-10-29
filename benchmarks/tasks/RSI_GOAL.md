# Goal: solve as close to 100% of the benchmark as possible

Benchmark has 60 physics model-building tasks. Every task is solvable: its reference solution scores reward 1.0 through the same harness an agent uses. **The goal is a self-improving agent system that solves as close to 100% of the 60 tasks as possible.** A task counts as solved at reward 1.0.

## 1. Run the benchmark

```bash
git clone https://github.com/exon-research/AI-Theorist.git && cd AI-Theorist
bash scripts/setup.sh                         # venv, Harbor, runtime Docker image
source .venv/bin/activate
python -m benchmarks.verify                   # 60 tasks verified (60 frozen)
python -m benchmarks.cpu_pinning &            # needed for the 10 fresh-start (QE) tasks
set -a; . ./.env; set +a                      # OPENAI_API_KEY, ANTHROPIC_API_KEY

# Reference solutions (sanity check, no API cost):
python -m benchmarks.run --policy oracle --suite all --campaign oracle --execute

# A model on all 60 tasks:
python -m benchmarks.run --policy <model> --suite all --campaign <name> --concurrency 2 --wait --execute
python -m benchmarks.tools.pilot_report <name>   # reward and checks passed per task
```

Each episode runs in a Harbor container. The agent works in a shell, buys calculations from the task's lab, and submits `submission.json` and `report.md`. The grade is numerical: reward 1.0 means every check passed. Use `--task <id>` to run selected tasks, `--replicate N` for repeats, and `--policy oracle` to see a reference solution in action.

To evaluate your own agent, use one of these:
- **As a Harbor agent:** pass it the same way `benchmarks.api_agent:PaperBaselineAgent` is passed in `benchmarks/run.py`.
- **As a model:** add it to `MODELS` in `run.py` and to the rate tables in `api_agent.py`.

## 2. First, benchmark GPT-6-astra and Opus 5.5

```bash
python -m benchmarks.run --policy baselines --suite all --campaign baselines --concurrency 2 --wait --execute
```

`baselines` runs `astra` (`openai/gpt-6-astra`) and `opus55` (`anthropic/claude-opus-5-5`). Confirm the GPT-6-astra API model id before the first run. Record their scores; they are the starting point the loop must beat.

For reference, the development pilots solved 13/60 (Sonnet 5) and 12/60 (GPT-5.6 Terra).

## 3. Improve

Run the self-improvement loop on the tasks and iterate until the solve rate is as close to 100% as you can get. Report the solve rate of each iteration and of the final system, alongside the two baselines.

One rule keeps the score meaningful: the agent and the loop may use only what the agent sees in a task:
- the instruction;
- the lab's catalog and purchased results;
- its own transcripts and grades.

They must never read the answer files: `environment/lab/private/`, `solution/`, `design.json`, `rubric.md`, `SOURCE.md` and `benchmarks/construction/`.
