"""Check a checkout against the frozen the benchmark manifest (evidence/release.json)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def verify():
    from benchmarks.tools.release import verify_release
    verify_release()


if __name__ == '__main__':
    verify()
