"""Command-line client for the task's calculation service (the benchmark).

  python3 /app/lab_client.py task                    # task description and budget
  python3 /app/lab_client.py status                  # budget used, running job
  python3 /app/lab_client.py evaluate '<json>'       # request a calculation
  python3 /app/lab_client.py job/<id>                # poll a calculation
  python3 /app/lab_client.py submit "$(cat /app/submission.json)"   # terminal
  python3 /app/lab_client.py submission              # poll grading of the submission
"""
import json
import os
import sys
import urllib.error
import urllib.request


def call(endpoint, payload=None):
    url = os.environ.get('LAB_URL', 'http://lab:8000') + '/' + endpoint
    data = None if payload is None else json.dumps(payload).encode()
    request = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        raise ValueError(json.load(exc).get('error', 'calculation service error')) from None


if __name__ == '__main__':
    try:
        print(json.dumps(call(sys.argv[1], json.loads(sys.argv[2]) if len(sys.argv) > 2 else None), allow_nan=False))
    except (ValueError, IndexError) as exc:
        print(json.dumps({'error': str(exc) or 'usage: lab_client.py ENDPOINT [JSON]'}))
        raise SystemExit(1)
