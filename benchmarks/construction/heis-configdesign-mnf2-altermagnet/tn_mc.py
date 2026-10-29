"""Reference classical Monte Carlo Neel temperature for an MnF2 exchange model (runs mc_tn.c).

    python3 -B tn_mc.py J.json OUT.json [--sizes 8 12 16] [--temps-K 110 ...] [--neq N] [--nmeas N] [--procs 4]

J.json holds {"J_meV": [J1, J2, ...]} for the shells in model.SHELLS order (sublattice-symmetrised: every bond
of a shell carries the same J). T_N is located from the crossings of the Binder ratio
g = <m^4>/<m^2>^2 of the staggered magnetisation between successive sizes.
"""
import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
KB = 0.08617333262


def bonds_file(J, path):
    lines = []
    for Js, s in zip(J, M.SHELLS):
        for d in s['vectors']:
            for site in (0, 1):
                if s['kind'] == 'same':
                    lines.append((site, site, *np.round(d).astype(int), Js))
                elif site == 0:     # A at tau_A: neighbour B at d, cell = d + tau_A - tau_B
                    lines.append((0, 1, *np.round(d + M.TAU['A'] - M.TAU['B']).astype(int), Js))
                else:               # B at origin: neighbour A at -d (inversion-symmetric set), cell = -d - tau_A + tau_B
                    lines.append((1, 0, *np.round(-d - M.TAU['A'] + M.TAU['B']).astype(int), Js))
    path.write_text(''.join('%d %d %d %d %d %.10f\n' % l for l in lines))


def compile_mc(workdir):
    exe = workdir/'mc_tn'
    subprocess.run(['gcc', '-O3', '-march=x86-64', '-o', str(exe), str(HERE/'mc_tn.c'), '-lm'], check=True)
    return exe


def run_point(exe, bonds, L, T_K, neq, nmeas, seed):
    out = subprocess.run([str(exe), str(bonds), str(L), '%.8f' % (T_K * KB), str(neq), str(nmeas), str(seed)],
                         capture_output=True, text=True, check=True).stdout
    blocks = np.array([[float(x) for x in line.split()] for line in out.strip().splitlines()])
    g_blocks = blocks[:, 1] / blocks[:, 0] ** 2
    m2, m4 = blocks[:, 0].mean(), blocks[:, 1].mean()
    # jackknife error of g over blocks
    n = len(blocks)
    jk = np.array([np.delete(blocks[:, 1], k).mean() / np.delete(blocks[:, 0], k).mean() ** 2 for k in range(n)])
    g = m4 / m2 ** 2
    return {'L': L, 'T_K': T_K, 'g': g, 'g_err': float(np.sqrt((n - 1) * np.var(jk))), 'm2': m2,
            'E': blocks[:, 2].mean(), 'blocks': n}


def crossings(points, sizes):
    out = []
    for La, Lb in zip(sizes[:-1], sizes[1:]):
        pa = sorted([p for p in points if p['L'] == La], key=lambda p: p['T_K'])
        pb = sorted([p for p in points if p['L'] == Lb], key=lambda p: p['T_K'])
        T = np.array([p['T_K'] for p in pa])
        diff = np.array([a['g'] - b['g'] for a, b in zip(pa, pb)])
        err = np.array([np.hypot(a['g_err'], b['g_err']) for a, b in zip(pa, pb)])
        k = int(np.argmax(np.sign(diff[:-1]) != np.sign(diff[1:]))) if np.any(np.sign(diff[:-1]) != np.sign(diff[1:])) else None
        if k is None:
            out.append({'sizes': [La, Lb], 'T_cross_K': None})
            continue
        # weighted linear fit of the difference over the two bracketing points and their neighbours
        sel = slice(max(k - 1, 0), min(k + 3, len(T)))
        c = np.polyfit(T[sel], diff[sel], 1, w=1 / err[sel])
        lin = float(T[k] + (T[k + 1] - T[k]) * diff[k] / (diff[k] - diff[k + 1]))
        out.append({'sizes': [La, Lb], 'T_cross_K': float(-c[1] / c[0]), 'T_cross_bracket_K': lin,
                    'bracket': [float(T[k]), float(T[k + 1])]})
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('jfile')
    p.add_argument('out')
    p.add_argument('--sizes', type=int, nargs='+', default=[8, 12, 16])
    p.add_argument('--temps-K', type=float, nargs='+', required=True)
    p.add_argument('--neq', type=int, default=4000)
    p.add_argument('--nmeas', type=int, default=20000)
    p.add_argument('--procs', type=int, default=4)
    p.add_argument('--workdir', default='/tmp/mc')
    a = p.parse_args()
    J = json.loads(Path(a.jfile).read_text())['J_meV']
    work = Path(a.workdir)
    work.mkdir(parents=True, exist_ok=True)
    exe = compile_mc(work)
    bonds = work/'bonds.txt'
    bonds_file(J, bonds)
    jobs = [(L, T) for L in a.sizes for T in a.temps_K]
    with ThreadPoolExecutor(a.procs) as pool:
        points = list(pool.map(lambda job: run_point(exe, bonds, job[0], job[1], a.neq, a.nmeas,
                                                     1000 * job[0] + int(job[1] * 10)), jobs))
    result = {'J_meV': J, 'mean_field_K': M.neel_mean_field_K(J), 'points': points,
              'crossings': crossings(points, a.sizes)}
    Path(a.out).write_text(json.dumps(result, indent=1))
    print(json.dumps({'mean_field_K': result['mean_field_K'], 'crossings': result['crossings']}, indent=1))
    for pt in points:
        print(pt['L'], pt['T_K'], '%.4f +- %.4f' % (pt['g'], pt['g_err']))


if __name__ == '__main__':
    main()
