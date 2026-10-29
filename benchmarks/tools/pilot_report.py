"""Summarize development pilots of the benchmark tasks: length, actual cost and GPT-6-astra-equivalent cost.

    python -m benchmarks.tools.pilot_report CAMPAIGN [CAMPAIGN ...] [--json OUT]

Token counts come from the pilot model's tokenizer and are a proxy for GPT-6-astra.
Astra-equivalent cost reprices them at $10 per M uncached input, $1 per M cached
input and $50 per M output (reasoning included).
"""
import argparse
from datetime import datetime
import json
from benchmarks.verify import ROOT

ASTRA = (10, 1, 50)


def load(path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def episode(result_path):
    trial = result_path.parent
    result = load(result_path)
    usage = load(trial/'agent/usage.json')
    meta = usage.get('metadata', {})
    state = load(trial/'artifacts/lab-results/state.json')
    details = load(trial/'verifier/details.json')
    final = state.get('submission') or {}
    started, ended = result.get('started_at'), result.get('finished_at')
    minutes = ((datetime.fromisoformat(ended) - datetime.fromisoformat(started)).total_seconds() / 60
               if started and ended else None)
    n_in, n_cache, n_out = (usage.get(k) or 0 for k in ('n_input_tokens', 'n_cache_tokens', 'n_output_tokens'))
    astra = ((n_in - n_cache) * ASTRA[0] + n_cache * ASTRA[1] + n_out * ASTRA[2]) / 1e6
    report = final.get('report') or {}
    return {'job': trial.parent.name, 'task': (result.get('task_name') or '').split('/')[-1],
            'model': meta.get('model') or (result.get('agent_info') or {}).get('name'),
            'reward': ((result.get('verifier_result') or {}).get('rewards') or {}).get('reward'),
            'submitted': bool(final), 'accepted': final.get('success'),
            'checks': report.get('checks'),
            'checks_passed': (f"{sum(map(bool, report['checks'].values()))}/{len(report['checks'])}"
                              if isinstance(report.get('checks'), dict) else None),
            'api_calls': meta.get('api_calls'), 'shell_commands': len(list((trial/'agent').glob('command-*-launch.json'))),
            'lab_evaluations': len(state.get('trace', [])), 'charge_used': state.get('used'), 'budget': state.get('budget'),
            'minutes': round(minutes, 1) if minutes is not None else None,
            'input_tokens': n_in, 'cached_tokens': n_cache, 'output_tokens': n_out,
            'cost_usd': round(usage.get('cost_usd') or 0, 3), 'astra_equiv_usd': round(astra, 2),
            'termination': meta.get('termination_reason'),
            'exception': (result.get('exception_info') or {}).get('exception_type'),
            'verifier_error': details.get('error')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaigns', nargs='+')
    p.add_argument('--json')
    a = p.parse_args()
    rows = [episode(r) for c in a.campaigns for r in sorted((ROOT/'runs'/c).glob('*/*/result.json'))]
    cols = ['task', 'model', 'reward', 'checks_passed', 'api_calls', 'lab_evaluations', 'charge_used', 'budget', 'minutes',
            'cost_usd', 'astra_equiv_usd', 'termination', 'exception']
    print('| ' + ' | '.join(cols) + ' |\n|' + '---|' * len(cols))
    for r in rows:
        print('| ' + ' | '.join('' if r[c] is None else str(r[c]) for c in cols) + ' |')
    for r in rows:
        failed = [k for k, v in (r['checks'] or {}).items() if not v]
        if failed:
            print(r['task'], r['model'], 'failed checks:', failed)
    if a.json:
        (ROOT/a.json).write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main()
