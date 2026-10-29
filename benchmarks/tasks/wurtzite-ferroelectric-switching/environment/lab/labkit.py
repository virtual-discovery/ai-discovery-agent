"""Shared calculation service for the benchmark tasks. Vendored unchanged into each task.

The task's lab directory provides `task_lab.py` with:

    validate(request) -> int
        Check a request and return its charge (positive integer). Raise ValueError
        for an invalid request; invalid requests are not charged.
    run(request, workdir) -> dict
        Execute an accepted request and return the result shown to the agent.
        It must never contain hidden targets. Exceptions mark the job failed;
        failed jobs stay charged. A ValueError message is shown to the agent.
    grade(payload, trace, validation) -> (bool, dict)
        Grade the terminal submission against hidden data. `trace` holds the
        executed records (request, status, result); `validation` holds the
        results of validation_requests. Must be a pure function of its inputs
        and the task's private data so that ablations can be graded offline.
    validation_requests(payload, trace) -> list   (optional)
    run_private(request, workdir) -> dict          (optional)
        Hidden calculations executed after submission. Prefer precomputed hidden
        references in private data; use these only when the check depends on the
        submission itself.

and `public.json`, served at GET /task, with at least `task_id` and `solver_budget`
(optional `reserved_validation_charge`). Private data lives under `private/`.
"""
import copy
import json
import math
import os
from pathlib import Path
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = Path(os.environ.get('LAB_ROOT', '/lab'))
OUT = Path(os.environ.get('LAB_RESULTS', '/lab-results'))
SYNC_WAIT = float(os.environ.get('LAB_SYNC_WAIT', '20'))
MAX_REQUEST_BYTES = 65536


def finite_json(value):
    """Round-trip through strict JSON so NaN/inf or non-JSON types fail loudly."""
    return json.loads(json.dumps(value, allow_nan=False))


class Lab:
    def __init__(self, task):
        self.task = task
        self.public = json.loads((ROOT/'public.json').read_text())
        self.budget = int(self.public['solver_budget'])
        self.reserved = int(self.public.get('reserved_validation_charge', 0))
        self.lock = threading.RLock()
        self.state = {'task_id': self.public['task_id'], 'budget': self.budget, 'used': 0, 'finished': False,
                      'trace': [], 'submission': None, 'validation': []}
        OUT.mkdir(parents=True, exist_ok=True)
        self.save()

    def save(self):
        tmp = OUT/'state.tmp'
        tmp.write_text(json.dumps(self.state, indent=1, allow_nan=False))
        tmp.replace(OUT/'state.json')

    def busy(self):
        return any(r['status'] == 'running' for r in self.state['trace'])

    def record_view(self, record):
        view = {k: record[k] for k in ('id', 'status', 'charge', 'request') if k in record}
        if record['status'] == 'complete':
            view['result'] = record['result']
        elif record['status'] == 'failed':
            view['error'] = record['error']
        else:
            view['instruction'] = 'Running: poll job/%d; do not restart the calculation.' % record['id']
        view.update(cumulative_charge=self.state['used'], budget=self.budget)
        return view

    def evaluate(self, request):
        with self.lock:
            if self.state['finished']:
                raise ValueError('Episode submitted; no further calculations')
            if self.busy():
                raise ValueError('A calculation is running; poll its job and do not restart it')
            charge = self.task.validate(request)
            if type(charge) is not int or charge < 1:
                raise RuntimeError('task_lab.validate must return a positive integer charge')
            if self.state['used'] + charge + self.reserved > self.budget:
                raise ValueError('Budget exhausted (remaining %d, request costs %d)'
                                 % (self.budget - self.reserved - self.state['used'], charge))
            self.state['used'] += charge
            record = {'id': len(self.state['trace']) + 1, 'request': finite_json(request), 'charge': charge,
                      'status': 'running', 'started_epoch': time.time()}
            self.state['trace'].append(record)
            self.save()
        worker = threading.Thread(target=self.execute, args=(record,), daemon=True)
        worker.start()
        worker.join(SYNC_WAIT)
        with self.lock:
            return self.record_view(record)

    def execute(self, record):
        workdir = OUT/('job-%03d' % record['id'])
        workdir.mkdir(parents=True, exist_ok=True)
        try:
            result = finite_json(self.task.run(copy.deepcopy(record['request']), workdir))
            if not isinstance(result, dict):
                raise RuntimeError('task_lab.run must return a JSON object')
            update = {'status': 'complete', 'result': result}
        except Exception as exc:  # the job fails; the episode continues
            (workdir/'error.log').write_text(traceback.format_exc())
            message = str(exc)[:500] if isinstance(exc, ValueError) else 'calculation failed (%s)' % type(exc).__name__
            update = {'status': 'failed', 'error': message}
        with self.lock:
            record.update(update, ended_epoch=time.time())
            self.save()

    def job(self, job_id):
        with self.lock:
            for record in self.state['trace']:
                if record['id'] == job_id:
                    return self.record_view(record)
        raise KeyError('Unknown job')

    def submit(self, payload):
        with self.lock:
            if self.state['finished']:
                raise ValueError('Episode already submitted')
            if self.busy():
                raise ValueError('Wait for the running calculation before submitting')
            self.state['finished'] = True
            self.state['submission'] = {'payload': finite_json(payload), 'status': 'validating',
                                        'submitted_epoch': time.time()}
            self.save()
        worker = threading.Thread(target=self.validate_submission, args=(payload,), daemon=True)
        worker.start()
        worker.join(SYNC_WAIT)
        return self.submission_view()

    def validate_submission(self, payload):
        try:
            trace = copy.deepcopy([r for r in self.state['trace'] if r['status'] != 'running'])
            requests = getattr(self.task, 'validation_requests', lambda p, t: [])(copy.deepcopy(payload), trace)
            results = []
            for i, request in enumerate(requests):
                workdir = OUT/('validation-%02d' % i)
                workdir.mkdir(parents=True, exist_ok=True)
                results.append(finite_json(self.task.run_private(request, workdir)))
                with self.lock:
                    self.state['validation'].append({'request': request, 'result': results[-1]})
                    self.save()
            success, report = self.task.grade(copy.deepcopy(payload), trace, results)
            update = {'status': 'complete', 'success': bool(success), 'report': finite_json(report)}
        except Exception as exc:
            (OUT/'grading-error.log').write_text(traceback.format_exc())
            update = {'status': 'failed', 'success': False, 'report': {'grading_failed': type(exc).__name__}}
        with self.lock:
            self.state['submission'].update(update, graded_epoch=time.time())
            self.save()

    def submission_view(self):
        with self.lock:
            sub = self.state['submission']
            if sub is None:
                raise KeyError('Nothing submitted')
            view = {'episode_finished': True, 'status': sub['status']}
            if sub['status'] == 'validating':
                view['instruction'] = 'Validation running: poll submission until it completes.'
            else:
                view['accepted'] = sub['success']
            return view

    def status(self):
        with self.lock:
            return {'used': self.state['used'], 'budget': self.budget, 'reserved': self.reserved,
                    'finished': self.state['finished'], 'running': self.busy()}


def handler_for(lab):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def respond(self, body, status=200):
            data = json.dumps(body, allow_nan=False).encode()
            self.send_response(status)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            try:
                if self.path == '/task':
                    self.respond(lab.public)
                elif self.path == '/status':
                    self.respond(lab.status())
                elif self.path.startswith('/job/'):
                    self.respond(lab.job(int(self.path.rsplit('/', 1)[1])))
                elif self.path == '/submission':
                    self.respond(lab.submission_view())
                else:
                    self.respond({'error': 'Unknown endpoint'}, 404)
            except (KeyError, ValueError) as exc:
                self.respond({'error': str(exc).strip("'")}, 404)

        def do_POST(self):
            try:
                size = int(self.headers.get('Content-Length', 0))
                if not 0 < size <= MAX_REQUEST_BYTES:
                    raise ValueError('Invalid request size')
                body = json.loads(self.rfile.read(size))
                if not isinstance(body, dict):
                    raise ValueError('JSON object required')
                if self.path == '/evaluate':
                    self.respond(lab.evaluate(body))
                elif self.path == '/submit':
                    self.respond(lab.submit(body))
                else:
                    raise ValueError('Unknown endpoint')
            except (ValueError, KeyError, TypeError) as exc:
                self.respond({'error': str(exc).strip("'")}, 400)
    return Handler


def main():
    sys.path.insert(0, str(ROOT))
    import task_lab
    lab = Lab(task_lab)
    ThreadingHTTPServer(('0.0.0.0', int(os.environ.get('LAB_PORT', '8000'))), handler_for(lab)).serve_forever()


if __name__ == '__main__':
    main()
