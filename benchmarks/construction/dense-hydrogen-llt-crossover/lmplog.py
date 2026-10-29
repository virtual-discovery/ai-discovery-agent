"""Parse LAMMPS thermo output (log.lammps / slurm-*.out) of the Tenti et al. runs."""
import re
from pathlib import Path
import numpy as np

COLS = ['Step', 'Temp', 'PotEng', 'KinEng', 'TotEng', 'Press', 'Volume']


def blocks(text):
    """Yield (header, rows) for every thermo block in a LAMMPS output."""
    lines = text.split('\n')
    i = 0
    while i < len(lines):
        if lines[i].split()[:2] == ['Step', 'Temp']:
            head = lines[i].split()
            rows = []
            i += 1
            while i < len(lines):
                parts = lines[i].split()
                if lines[i].startswith('WARNING'):
                    i += 1
                    continue
                if len(parts) != len(head):
                    break
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    break
                i += 1
            yield head, np.array(rows)
        i += 1


def run_series(files):
    """Stitch every thermo block of a run's output files: returns rows sorted by step (later jobs win)."""
    def jobid(p):
        m = re.search(r'slurm-(\d+)', p.name)
        return int(m.group(1)) if m else 10**12   # log.lammps duplicates the last job
    segs = []
    for p in sorted(files, key=jobid):
        for head, rows in blocks(Path(p).read_text(errors='replace')):
            if len(rows):
                idx = [head.index(c) for c in COLS]
                segs.append((jobid(p), rows[:, idx]))
    return segs
