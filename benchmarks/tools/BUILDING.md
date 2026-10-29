# Building a benchmark task

This is the build spec for one Experiment 1 task. The goal is a long agentic physics task, similar in length to a GeneBench-Pro problem: a few dollars per GPT-6-astra attempt, which is about 40–80 agent tool calls and 1–3 hours.

- **Designs:** the survey records `benchmarks/expansion/candidates.json` and `DETAILS.md`, now removed from the tree (see git history at commit 7d0bc84).
- **Worked example:** `benchmarks/tasks/heis-configdesign-fe`, with its construction scripts in `benchmarks/construction/heis-configdesign-fe`. Read every file of it before you start.

## What you deliver

1. **The task package:** `benchmarks/tasks/<id>/`, created by `new_task.py`, which must pass `check_task`.
2. **The construction folder:** `benchmarks/construction/<id>/`. It holds:
   - scripts that refetch the pinned source files (URL plus checksum);
   - scripts that regenerate every staged, private and ablation file;
   - `check.json`, written by the gate;
   - `README.md`, with build notes: deviations from the design and why, the calibration table, measured oracle length, and open risks.
3. **A Harbor oracle run with reward 1.0** (command below).
4. **A short final report** to the coordinator, in your last message:
   - task id;
   - whether the gate passed;
   - the oracle reward;
   - budget and the oracle's evaluations;
   - the number of decision points;
   - the tightest pass/fail margins;
   - peak disk used;
   - deviations from the design;
   - risks;
   - anything you could not do.

If the design is infeasible with the released data, adapt it. Keep the physics question and at least 5 genuine decision points, and record the change. If no honest version is possible, stop and report why. Never fake data, targets or provenance.

## Hard rules

- **Write only inside your own paths:**
  - `benchmarks/tasks/<id>/`
  - `benchmarks/construction/<id>/`
  - `runs/check/<id>/`
  - `runs/build-<id>/`
  - `~/.cache/ai-theorist/<id>/`
  - a subfolder of your scratchpad named after your task id (the scratchpad is shared with other builders)
- **Do not touch shared files.** That means `benchmarks/tools/*`, `run.py`, `api_agent.py`, the v1 tasks and other builders' tasks. If shared tooling blocks you, say so in your report.
- **No git operations** (commit, push, branch, stash, checkout). The coordinator commits.
- **No `__pycache__` in the task package.** Run Python with `-B` or `PYTHONDONTWRITEBYTECODE=1`, and never import task modules in place; the release check rejects caches.
- **No paid model APIs.** Do not run the `sonnet`, `terra`, `opus`, `gpt` or `gemini` policies; the coordinator runs the pilots. Never read or print `~/.env`.
- **Process safety:**
  - Never use `pkill -f` or `killall`. Kill only PIDs you started, found with `pgrep` and killed in a separate command.
  - Never stop other containers.
  - Never wait with `until ! pgrep -f <pattern>` loops: the loop matches its own command line and never ends. Wait on a PID (`while kill -0 $pid`), or use `docker wait <container>`.
  - Never run `docker system prune` or `docker image prune -a`.
- **Disk is scarce: about 26 GB free on the whole host, shared by six builders.**
  - Keep at most 3 GB of downloads at a time, under `~/.cache/ai-theorist/<id>/`.
  - Prefer partial fetches, such as HTTP range reads of a remote zip's central directory and of single members, or the archive API's per-file links.
  - When the package passes, delete the raw downloads. Keep the refetch scripts and checksums.
  - Never write large files into the repository. The task package must be at most 60 MB, and no file may exceed 40 MB.
- **CPU:**
  - Use only the cpuset you were given:
    - `check_task --cpuset <yours>`;
    - `docker run --cpuset-cpus <yours>` for any construction compute.
  - Run containers with `--user $(id -u):$(id -g)` so they do not leave root-owned files.
  - Construction compute runs in the runtime image `ai-theorist-benchmark-runtime:20261007`, not on the host.

## Runtime facts

- **Lab and verifier image:** `ai-theorist-benchmark-runtime:20261007`, Ubuntu 22.04.
  - Python 3.10, numpy 1.21.5, scipy 1.8.0, requests.
  - Quantum ESPRESSO 6.7, all `*.x` executables, including ph.x, matdyn.x, q2r.x, epw.x, projwfc.x, pp.x and dos.x.
  - OpenMPI `mpirun`.
  - No ASE, spglib, phonopy, pymatgen, Wannier90, LAMMPS or torch.
- **Version compatibility.** `task_lab.py` and `solution/*.py` must run under these versions. The gate grades inside the lab image, so host-only numpy features will fail there.
- **Agent image.** Same base with Quantum ESPRESSO removed: python3, numpy and scipy only, and no network. All calculations the agent pays for go through the lab service, so the budget is binding.
- **Extra lab dependencies.** Pinned pip or apt installs in the task's `Dockerfile.lab` are allowed only if unavoidable. Prefer precomputing an evaluator's outputs at construction time and serving them as lookups. Record any extra dependency, with versions and image size, in `design.json` under `runtime_dependencies`. Keep the extra layer under 1 GB.
- **Async jobs.** `labkit` runs each request in a thread. Calls return within 20 s with either the result or a running job id, and the agent polls `job/<id>`. A computing lab should keep each evaluation to 10 minutes or less on the task's `cpus`; use `mpirun -np <cpus>`. Set `OMPI_MCA_btl_vader_single_copy_mechanism=none` (the template already does).

## Design standard (GeneBench-style)

### Prompt (`instruction.md`)

- **Include:** the scientific question; the lab interface (request kinds, costs, budget); any conventions needed to make the graded numbers unique (units, sign and counting conventions, definitions); the deliverable JSON with placeholder values; and what the report must explain.
- **Leave out:** the workflow, the traps, and hints that name a decision point. Follow the Fe example's structure and length.

### Decision points

- **At least 5, ideally 6–8,** each a real inferential fork that can be resolved from evidence the agent can see. Typical forks:
  - conventions and units;
  - unconverged or aliased data;
  - choice of model class;
  - underdetermined fits;
  - finite-size or k-mesh convergence;
  - symmetry;
  - which calculation discriminates between competing explanations;
  - transfer to a regime where the model breaks down.
- **Ablations.** Each decision point needs at least one ablation in `solution/ablations.py`: the realistic wrong path, computed from the real data by construction code (store the values in `solution/ablation_values.json`). Arbitrary perturbations do not count. Every ablation must fail grading.
- **Robustness.** List at least 2 defensible alternative analyses in `ROBUSTNESS`, for example a different shell count or fit window, or another valid estimator; they must pass.

### Grading

- **Tolerances** sit between the spread of defensible methods and the nearest wrong path.
- **Calibration table.** Record it in `design.json` under `"calibration"`: for each graded field, the reference value, the robustness values, the nearest ablation, and the tolerance.
- **Rounding.** Rounding every number to 4 significant digits must still pass.

### Hidden targets

- **Recoverable.** Uniquely recoverable from data the agent can buy within budget, and never served by the lab.
- **Not leaked.** Not present in any agent-visible file. The gate scans `hidden_values.json` at 4–6 significant digits, so list every graded hidden number there. Also avoid derivable leaks, such as the same value under another key or a trivial transform of it.
- **Not memorized.** Not a famous headline number. Grade protocol-specific values instead.
- **Recoverable without recall.** If your oracle starts from or uses literature numbers (e.g. a paper's fitted parameters as a fit's starting point), show that a route without them reaches every graded window with compute feasible in the agent container; otherwise grade only what that route recovers (lesson from MnWO4).
- **Fresh starts.** Truth comes from your own reference calculations, run with the same lab code path the agent uses, and stored only under `private/`.
- **Released starts.** Truth comes from author-archived values that the lab withholds.

### Difficulty (lesson from the Fe pilots)

Sonnet 5 and GPT-5.6 Terra both solved the Fe task on their first attempt (Sonnet: 55 calls; Terra: 24 calls and 6 minutes). They handled the cone factor, 2π, shell aliasing and the separate large-angle model without help. Current models know textbook conventions and standard fitting practice, so those make poor decision points.

Build difficulty from these instead:

- **Evidence that has to be discovered, not recalled.** For example:
  - released calculations that are inconsistent in a way visible only in their metadata or raw outputs (smearing, k-mesh, cutoff, volume, unconverged runs);
  - a confounder that mimics the target effect;
  - aliasing or finite-size artefacts that only a cross-check exposes.
- **Discriminating experiments under a binding budget.** The agent must decide which calculation separates the competing explanations. The consequence of a choice should be visible only after spending budget on it.
- **Cascades.** An early choice silently biases later graded quantities.
- **Regime judgement.** Which data are converged or asymptotic, where a model stops applying, and what the data cannot support.

Two limits apply:

- **Uniqueness.** Never create difficulty through ambiguity: the graded answer must stay uniquely recoverable.
- **Tolerances.** Never tighten them below the robustness spread to manufacture failures, or loosen them past the calibration rule.

The coordinator pilots every task with both development models. A task both models solve is flagged *dev-easy*. If too many tasks are dev-easy, the easiest ones get a hardening pass.

### Budget and length

- **Targets:** 15–60 lab evaluations and 40–80 agent tool calls.
- **Budget:** it must force choices, so buying everything should be impossible (about half the catalog or less). The oracle should use at most about 90% of it.
- **Recording:** put the oracle's evaluation count and an estimated tool-call count in `design.json` under `length_estimate`.

### Grading function and rubric

- **Grading.** `grade(payload, trace, validation)` is pure and deterministic. It requires `evidence_ids` to be a non-empty subset of completed job ids, and it returns a per-check report.
- **Rubric.** `rubric.md` has at least 4 numbered, checkable statements, covering the competing explanations in the design.

### Provenance

`SOURCE.md` records:

- the citation and DOI;
- the archive URL and license (open licenses only; CC BY-SA must be stated per file);
- the files used, with checksums;
- every transformation;
- what is hidden.

## Procedure

1. **Read.** Read the Fe example (all files), `labkit/labkit.py`, `check_task.py` and your candidate's record and DETAILS section. For computing labs, the fresh-start tasks show Quantum ESPRESSO usage, for example `benchmarks/tasks/fresh-inp-polar-phonons/environment/lab/`.
2. **Fetch and explore.** Fetch the minimal files, explore the data, and confirm feasibility before writing code.
3. **Scaffold** (task id = catalog id = candidate id):
   ```sh
   cd ~/Desktop/workspace/AI-Theorist
   .venv/bin/python -m benchmarks.tools.new_task <id> --start released|fresh --catalog-id <id> \
       --source-group <source_group> --title "<short title>" [--cpus N --memory-mb M --agent-timeout S]
   ```
4. **Construct.** Write the construction scripts, then `task_lab.py`, `solution/solve.py`, `solution/ablations.py`, `instruction.md`, `rubric.md`, `SOURCE.md` and `design.json` (fill every field). Make sure `public.json` and `tests/public.json` agree on `solver_budget`.
5. **Gate.** Iterate until it passes. Never weaken a check or claim an exception you cannot justify in `design.json`:
   ```sh
   .venv/bin/python -m benchmarks.tools.check_task <id> --cpuset <yours>
   ```
6. **Harbor oracle.** The reward must be 1.0:
   ```sh
   source .venv/bin/activate && python -m benchmarks.run --policy oracle --development \
       --task <id> --suite all --campaign build-<id> --execute
   ```
   Read `runs/build-<id>/oracle-r1-<id>/*/verifier/reward.txt`. For tasks whose lab runs `pw.x`, the runner needs the coordinator's CPU-pinning daemon. If it refuses, report that instead of starting the daemon yourself.
7. **Clean up and report.** Delete raw downloads and `docker image rm` the per-trial Harbor images named `<id>__*`. Write the build notes, then report.
