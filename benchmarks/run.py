"""Preview or execute bounded Harbor campaigns from a fresh repository checkout."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import threading
from benchmarks.verify import ROOT
from benchmarks.tools.release import task_rows, verify_release

MODELS = {'opus55': 'anthropic/claude-opus-5-5', 'astra': 'openai/gpt-6-astra',
          'sonnet': 'anthropic/claude-sonnet-5', 'terra': 'openai/gpt-5.6-terra',
          'opus': 'anthropic/claude-opus-4-6', 'gpt': 'openai/gpt-5.3-codex', 'gemini': 'google/gemini-3.1-pro-preview'}
KEYS = {name: ('ANTHROPIC_API_KEY' if model.startswith('anthropic') else 'GEMINI_API_KEY' if model.startswith('google')
               else 'OPENAI_API_KEY') for name, model in MODELS.items()}
GROUPS = {'baselines': ['opus55', 'astra'], 'dev': ['sonnet', 'terra']}
TASKS = 'benchmarks/tasks'
CAP, CALLS = 20., 200   # per episode, the same for every model
SMOKE = ['heis-configdesign-fe', 'fresh-h-chain-peierls-vs-afm']


def select(rows, suite, tasks=None):
    if tasks:
        unknown = set(tasks) - {r['task_id'] for r in rows}
        if unknown:
            raise ValueError(f'Unknown tasks: {sorted(unknown)}')
        return [r for r in rows if r['task_id'] in tasks]
    return rows if suite == 'all' else [r for r in rows if r['task_id'] in SMOKE]


def completed_job(folder):
    result = folder/'result.json'
    return result.exists() and bool(json.loads(result.read_text()).get('finished_at'))


def cost_in_runs(folder):
    total = 0.
    for path in folder.glob('*/*/agent/usage.json'):
        total += float(json.loads(path.read_text()).get('cost_usd') or 0)
    return total


def build_command(row, policy, output, replicate, cap):
    task = row['task_id']
    name = f'{policy}-r{replicate}-{task}'
    cmd = ['harbor', 'run', '--path', TASKS, '--include-task-name', task,
           '--n-concurrent', '1', '--max-retries', '0', '--jobs-dir', str(output), '--job-name', name]
    if policy in MODELS:
        cmd += ['--agent', 'benchmarks.api_agent:PaperBaselineAgent', '--model', MODELS[policy],
                '--ak', f'max_calls={CALLS}', '--ak', f'max_cost_usd={cap}', '--ak', 'agent_qe=false',
                '--ak', 'max_output_tokens=32000', '--ak', 'truncation_retries=3']
    else:
        cmd += ['--agent', 'oracle']
    return name, cmd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--policy', required=True, choices=[*MODELS, *GROUPS, 'oracle'])
    p.add_argument('--development', action='store_true', help='Allow tasks not yet frozen in evidence/release.json')
    p.add_argument('--suite', choices=['smoke', 'all'], default='smoke')
    p.add_argument('--task', action='append')
    p.add_argument('--replicate', type=int, default=1)
    p.add_argument('--max-cost-per-run', type=float, default=CAP)
    p.add_argument('--other-spend-usd', type=float, default=0., help='Additional incurred API, rental, training, storage and transfer charges not already in this repository ledger')
    p.add_argument('--campaign', default='evaluation')
    p.add_argument('--concurrency', type=int, default=1, help='Independent Harbor jobs run at once within this campaign')
    p.add_argument('--wait', action='store_true', help='Wait for a running paid campaign instead of exiting')
    p.add_argument('--execute', action='store_true', help='Actually run; default only prints the plan')
    args = p.parse_args()
    if args.concurrency < 1 or args.replicate < 1 or not 0 < args.max_cost_per_run <= 2000 or not 0 <= args.other_spend_usd <= 2000:
        p.error('Invalid replicate or dollar limit')
    if not args.campaign or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_' for c in args.campaign):
        p.error('Campaign names may contain letters, digits, hyphens and underscores')
    rows = task_rows(require_frozen=not args.development)
    policies = GROUPS.get(args.policy, [args.policy])
    output = ROOT/'runs'/args.campaign
    plans = [(policy, row, *build_command(row, policy, output, args.replicate, args.max_cost_per_run))
             for policy in policies for row in select(rows, args.suite, args.task)]
    if not plans:
        p.error('No matching tasks; use --suite all or --task')
    for policy, row, name, cmd in plans:
        print(('SKIP completed ' if completed_job(output/name) else '') + shlex.join(cmd))
    pending = [(policy, row, name, cmd) for policy, row, name, cmd in plans if not completed_job(output/name)]
    reserve = sum(args.max_cost_per_run for policy, *_ in pending if policy in MODELS)
    print(f'{len(plans)} planned episodes; {len(pending)} pending; API cap reservation ${reserve:.2f} (not a token forecast).')
    if not args.execute:
        return
    if not pending:
        print('All selected jobs are terminal; nothing to rerun.')
        return
    verify_release(require_frozen=not args.development,
              only={row['task_id'] for _, row, *_ in pending} if args.development else None)
    if not shutil.which('harbor'):
        p.error('Activate .venv first')
    # Colima/Docker Desktop may not share an arbitrary temporary checkout path.
    # Test a file bind before creating any jobs; native Linux also checks the mount.
    probe = subprocess.run(['docker', 'run', '--rm', '--network', 'none', '--mount',
        f'type=bind,src={ROOT / "requirements.txt"},dst=/probe,readonly',
        'ai-theorist-benchmark-runtime:20261007', 'python3', '-c',
        'import hashlib; print(hashlib.sha256(open("/probe","rb").read()).hexdigest())'],
        capture_output=True, text=True)
    if probe.returncode or probe.stdout.strip() != hashlib.sha256((ROOT/'requirements.txt').read_bytes()).hexdigest():
        p.error('Docker cannot read this checkout or the runtime image is missing. Build the image and use a Docker-shared path; no jobs started.')
    from benchmarks.cpu_pinning import HEARTBEAT, qe_tasks
    if {row['task_id'] for _, row, *_ in pending} & set(qe_tasks()):
        import time
        fresh = HEARTBEAT.exists() and time.time() - json.loads(HEARTBEAT.read_text())['epoch'] < 15
        if not fresh:
            p.error('QE tasks need dedicated cores: start `python -m benchmarks.cpu_pinning` first; no jobs started.')
    (ROOT/'runs').mkdir(exist_ok=True)
    # The lock serializes paid campaigns for the budget check; unpaid controls skip it.
    lock_path = ROOT/'runs/.campaign.lock' if reserve else os.devnull
    with open(lock_path, 'w') as lock:
        if reserve:
            fcntl.flock(lock, fcntl.LOCK_EX | (0 if args.wait else fcntl.LOCK_NB))
        ledger = json.loads((ROOT/'benchmarks/costs.json').read_text())
        spent = sum(float(ledger[k]) for k in ['measured_development_cost_usd', 'prior_credential_check_usd', 'prior_rental_estimate_usd'])
        spent += args.other_spend_usd + sum(cost_in_runs(d) for d in (ROOT/'runs').iterdir() if d.is_dir())
        if spent + reserve > float(ledger['total_budget_usd']):
            p.error(f'Would exceed $2000: recorded/additional spend ${spent:.2f} + reserved API cap ${reserve:.2f}')
        # Validate credentials before starting any paid episode.
        from benchmarks.api_agent import credentials
        for policy in policies:
            if policy in MODELS:
                credentials(KEYS[policy])
        env = dict(os.environ, PYTHONPATH=str(ROOT))
        for policy, row, name, cmd in pending:
            if (output/name).exists():
                p.error(f'Unfinished job exists: {output/name}. Inspect it and use harbor job resume; do not restart a live calculation.')

        stop = threading.Event()

        def launch(job):
            policy, row, name, cmd = job
            if stop.is_set():
                return None
            result = subprocess.run(cmd, cwd=ROOT, env=env)
            trials = [json.loads(path.read_text()) for path in (output/name).glob('*/result.json')]
            if result.returncode or not trials or any(trial.get('exception_info') or not trial.get('finished_at') for trial in trials):
                stop.set()
                return name
        # After an execution error no new job starts; jobs already running finish.
        with ThreadPoolExecutor(args.concurrency) as pool:
            errors = [name for name in pool.map(launch, pending) if name]
        if errors:
            raise RuntimeError(f'Execution errors in {errors}; inspect saved results before continuing the campaign')
    print('Summary: python -m benchmarks.tools.pilot_report ' + args.campaign)


if __name__ == '__main__':
    main()
