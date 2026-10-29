"""Run the lab kinetic Monte Carlo engine (kmc.c, compiled at image build) for an agent-specified barrier model.

The engine is a rejection-free single-vacancy kMC on an L x L x L conventional bcc supercell (N = 2 L^3 sites,
N - 1 atoms and one vacancy). See kmc.c for the barrier expression. Trajectories are split over worker processes.
"""
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).parent
EXE = str(HERE/'kmc')
ELEMENTS = ['V', 'Cr', 'Nb', 'Mo', 'Ta', 'W']
SPECIES7 = ELEMENTS + ['Va']


def counts_for(fractions, L):
    """Atom counts filling N - 1 sites (largest-remainder rounding, ties broken by element order)."""
    n = 2 * L ** 3 - 1
    x = [fractions.get(e, 0.0) for e in ELEMENTS]
    tot = sum(x)
    raw = [xi / tot * n for xi in x]
    c = [int(math.floor(r)) for r in raw]
    order = sorted(range(6), key=lambda i: (-(raw[i] - c[i]), i))
    for i in order[:n - sum(c)]:
        c[i] += 1
    return c


def build_input(req, orbits, shells, ntraj, seed):
    m = req['model']
    kra, es = m['kra'], m.get('end_state')
    n_orb = kra.get('orbits', 0)
    K0 = [float(kra['K0'][e]) for e in ELEMENTS]
    lines = [f"{req['cell']} {req['steps']} {ntraj} {seed} {float(req['temperature_K'])} 1.0 "
             f"{req.get('equilibration_sweeps', 0)} {req.get('checkpoints', 0)} {0 if es is None else 1}",
             ' '.join(str(c) for c in counts_for(req['composition'], req['cell'])), ' '.join(repr(v) for v in K0), str(n_orb)]
    for o in orbits[:n_orb]:
        lines.append(f"{len(o)} " + ' '.join(str(x) for site in o for x in site))
    W = kra.get('W', {})
    lines.append(' '.join(repr(float(W[e][o][k])) if n_orb else '0' for e in ELEMENTS for o in range(n_orb) for k in range(6)))
    Q = kra.get('Q')
    lines.append(str(3 if Q else 0))
    lines.append(' '.join(repr(float(Q[e][i][j])) if Q else '0' for e in ELEMENTS for i in range(6) for j in range(6)))
    n_sh = 0 if es is None else es['shells']
    lines.append(str(n_sh))
    for s in shells[:n_sh]:
        lines.append(f"{len(s)} " + ' '.join(str(x) for site in s for x in site))
    lines.append(' '.join(repr(float(es['V'][s][a][b])) for s in range(n_sh) for a in range(7) for b in range(7)) if n_sh else '')
    return '\n'.join(lines) + '\n'


def run(req, orbits, shells, workers):
    ntraj = req['trajectories']
    chunks = [ntraj // workers + (1 if i < ntraj % workers else 0) for i in range(workers)]
    chunks = [c for c in chunks if c > 0]

    def one(i):
        text = build_input(req, orbits, shells, chunks[i], 1000003 * (req.get('seed', 0) + 1) + 7919 * i)
        p = subprocess.run([EXE], input=text, capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError(p.stderr[-300:])
        return p.stdout
    with ThreadPoolExecutor(len(chunks)) as pool:
        outs = list(pool.map(one, range(len(chunks))))
    traj, ck, sro = [], {}, []
    offset = 0
    for i, out in enumerate(outs):
        for line in out.splitlines():
            f = line.split()
            if f[0] == 'T':
                t = int(f[1]) + offset
                traj.append({'trajectory': t, 'time_ps': float(f[2]), 'vacancy_r2_a2': float(f[3]) / 4,
                             'sum_r2_a2': {e: float(f[4 + k]) / 4 for k, e in enumerate(ELEMENTS)},
                             'hops': {e: int(f[10 + k]) for k, e in enumerate(ELEMENTS)},
                             'mean_encountered_barrier_eV': float(f[16]), 'mean_accepted_barrier_eV': float(f[17]),
                             'equilibration_acceptance': float(f[18])})
            elif f[0] == 'C':
                t = int(f[1]) + offset
                ck.setdefault(t, []).append({'step': int(f[3]), 'time_ps': float(f[4]), 'vacancy_r2_a2': float(f[5]) / 4,
                                             'sum_r2_a2': {e: float(f[6 + k]) / 4 for k, e in enumerate(ELEMENTS)}})
            elif f[0] == 'S':
                sro.append([int(x) for x in f[2:]])
        offset += chunks[i]
    for t in traj:
        if t['trajectory'] in ck:
            t['checkpoints'] = ck[t['trajectory']]
    pairs = [[sum(s[6 * a + b] for s in sro) for b in range(6)] for a in range(6)]
    return traj, pairs
