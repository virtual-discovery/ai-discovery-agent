"""Numerical/evidence verifier for the benchmark (vendored unchanged into each task).

Reward 1 requires a terminal, graded, accepted submission whose saved file equals
the executed payload, within budget, with a report. Physical reasoning is reviewed
separately against rubric.md.
"""
import json
from pathlib import Path

out = Path('/logs/verifier')
out.mkdir(parents=True, exist_ok=True)
checks, score = {}, 0.
try:
    public = json.loads(Path('/tests/public.json').read_text())
    state = json.loads(Path('/lab-results/state.json').read_text())
    submitted = json.loads(Path('/app/submission.json').read_text())
    final = state.get('submission') or {}
    report = Path('/app/report.md')
    checks = {'task_matches': public['task_id'] == state['task_id'],
              'finished': bool(state['finished']),
              'within_budget': state['used'] <= public['solver_budget'],
              'executed_submission': final.get('payload') == submitted,
              'graded': final.get('status') == 'complete',
              'scientific_numerical_checks': bool(final.get('success')),
              'report_present': report.exists() and len(report.read_text().strip()) >= public.get('min_report_chars', 200)}
    score = float(all(checks.values()))
    checks.update(solver_charge=state['used'], details=final.get('report'),
                  interpretation_review='Physical reasoning is reviewed separately against rubric.md')
except Exception as exc:
    checks['error'] = type(exc).__name__
(out/'reward.txt').write_text(str(score))
(out/'details.json').write_text(json.dumps(checks, indent=2) + '\n')
print(json.dumps(checks))
