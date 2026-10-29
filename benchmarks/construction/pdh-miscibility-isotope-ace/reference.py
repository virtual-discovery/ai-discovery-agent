"""Reference calculations for pdh-miscibility-isotope-ace with the lab's own code path (ace_numpy + relax).

Run inside the runtime image (numpy 1.21 / scipy 1.8), e.g.
  docker run --rm --cpuset-cpus 6,7,22,23 --user $(id -u):$(id -g) -e OMP_NUM_THREADS=1 \
     -v $PWD:/c -v $CACHE/work:/w ai-theorist-benchmark-runtime:20261007 python3 -B /c/reference.py base
Parts: base | dilute | grids | random SEED NCELL "n1,n2,..." REPS | ordered
Each part writes ref/<part>[...].json in this folder.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ace_numpy import ACE, scan_atom   # noqa: E402
from relax import relax                 # noqa: E402

POT = sys.argv[sys.argv.index('--pot') + 1] if '--pot' in sys.argv else '/w/PdH_ace_ctilde.json'
ace = ACE(POT)
OUT = HERE/'ref'
OUT.mkdir(exist_ok=True)
BASE = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])


def fcc(a, n):
    pts = np.array([b + [i, j, k] for i in range(n) for j in range(n) for k in range(n) for b in BASE]) * a
    return np.eye(3) * a * n, pts


def base():
    r = relax(ace, np.eye(3) * 3.9, ['Pd'] * 4, BASE * 3.9, mode='cell', fmax=1e-4)
    a0 = float(np.cbrt(abs(np.linalg.det(r['cell']))))
    cell, pos = fcc(a0, 1)
    e0 = ace.compute(cell, ['Pd'] * 4, pos, forces=False)['energy']
    V0 = a0 ** 3

    def E(eps):
        F = np.eye(3) + eps
        return ace.compute(cell @ F.T, ['Pd'] * 4, pos @ F.T, forces=False)['energy']
    res = {'a0': a0, 'E_Pd_atom': e0 / 4}
    for d in (0.005, 0.01):
        c11 = (E(np.diag([d, 0, 0])) + E(np.diag([-d, 0, 0])) - 2 * e0) / d ** 2 / V0
        s = (E(np.diag([d, d, 0])) + E(np.diag([-d, -d, 0])) - 2 * e0) / d ** 2 / V0
        c12 = (s - 2 * c11) / 2
        sh = np.zeros((3, 3))
        sh[0, 1] = sh[1, 0] = d / 2
        c44 = (E(sh) + E(-sh) - 2 * e0) / d ** 2 / V0
        res[f'Cij_eV_A3_d{d}'] = [c11, c12, c44]
    h2 = relax(ace, np.eye(3) * 12, ['H', 'H'], [[0, 0, 0], [0.75, 0, 0]], fmax=1e-5)
    res['E_H2'] = h2['energy']
    res['d_H2'] = float(np.linalg.norm(h2['positions'][1] - h2['positions'][0]))
    # harmonic H2 frequency of the potential (for the record; the task supplies experimental gas data)
    dd = res['d_H2']
    es = [ace.compute(np.eye(3) * 12, ['H', 'H'], [[0, 0, 0], [dd + s, 0, 0]], forces=False)['energy'] for s in (-0.005, 0, 0.005)]
    kk = (es[0] + es[2] - 2 * es[1]) / 0.005 ** 2
    res['H2_harmonic_zpe_eV'] = 0.5 * 0.0646541 * np.sqrt(kk / (1.00782503 / 2))
    json.dump(res, open(OUT/'base.json', 'w'), indent=1)
    print(res)


def dilute():
    b = json.load(open(OUT/'base.json'))
    a0, EH2 = b['a0'], b['E_H2']
    res = {}
    for n in (1, 2, 3, 4):
        cell, pos = fcc(a0, n)
        N = len(pos)
        EPd = ace.compute(cell, ['Pd'] * N, pos, forces=False)['energy']
        for site, frac in (('O', [0.5, 0, 0]), ('T', [0.25, 0.25, 0.25])):
            P = np.vstack([pos, np.array(frac) * a0])
            for mode in ('positions', 'volume'):
                t = time.time()
                r = relax(ace, cell, ['Pd'] * N + ['H'], P, mode=mode)
                key = f'{n}_{site}_{mode}'
                res[key] = {'N': N, 'E_sol': r['energy'] - EPd - EH2 / 2,
                            'dV': abs(np.linalg.det(r['cell'])) - abs(np.linalg.det(cell)),
                            'p_GPa': -np.trace(r['stress']) / 3 * 160.21766208, 'max_force': r['max_force']}
                if n == 3:
                    res[key]['structure'] = {'cell': r['cell'].tolist(), 'positions': r['positions'].tolist()}
                print(key, {k: v for k, v in res[key].items() if k != 'structure'}, round(time.time() - t), 's', flush=True)
    json.dump(res, open(OUT/'dilute.json', 'w'))


def grids():
    d = json.load(open(OUT/'dilute.json'))
    out = {}
    for site, half, h in (('O', 0.84, 0.06), ('T', 0.45, 0.045)):
        s = d[f'3_{site}_volume']['structure']
        cell, pos = np.array(s['cell']), np.array(s['positions'])
        sym = ['Pd'] * (len(pos) - 1) + ['H']
        c = pos[-1]
        g = np.arange(-half, half + 1e-9, h)
        X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
        pts = c + np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
        t = time.time()
        e = scan_atom(ace, cell, sym, pos, len(pos) - 1, np.vstack([c, pts]))
        print(site, len(pts), round(time.time() - t), 's', flush=True)
        np.save(OUT/f'grid_{site}.npy', (e[1:] - e[0]).reshape(X.shape))
        dd = 0.01
        hp = [c + s * dd * np.eye(3)[ax] for ax in range(3) for s in (1, -1)]
        eh = scan_atom(ace, cell, sym, pos, len(pos) - 1, np.array(hp))
        out[site] = {'h': h, 'half': half, 'k_diag': [(eh[2 * a] + eh[2 * a + 1] - 2 * e[0]) / dd ** 2 for a in range(3)]}
    json.dump(out, open(OUT/'grids.json', 'w'), indent=1)


def random_part(seed, ncell, nhs, reps):
    b = json.load(open(OUT/'base.json'))
    a0, EH2 = b['a0'], b['E_H2']
    cell, pos = fcc(a0, ncell)
    N = len(pos)
    octa = pos + np.array([0.5, 0, 0]) * a0
    EPd = N * b['E_Pd_atom']
    rng = np.random.default_rng(seed)
    res = []
    for nh in nhs:
        for rep in range(reps):
            sel = rng.choice(N, nh, replace=False)
            P = np.vstack([pos, octa[sel]])
            sym = ['Pd'] * N + ['H'] * nh
            t = time.time()
            rc = relax(ace, cell, sym, P, mode='cell', max_steps=1500)
            rv = relax(ace, cell, sym, P, mode='volume', max_steps=1500)
            rf = relax(ace, cell, sym, P, mode='positions', max_steps=1500)
            row = {'N': N, 'nh': nh, 'x': nh / N, 'e_cell': (rc['energy'] - EPd - nh * EH2 / 2) / N,
                   'e_volume': (rv['energy'] - EPd - nh * EH2 / 2) / N, 'e_fixed': (rf['energy'] - EPd - nh * EH2 / 2) / N,
                   'v_cell': abs(np.linalg.det(rc['cell'])) / N, 'conv': [rc['converged'], rv['converged'], rf['converged']]}
            res.append(row)
            print(row, round(time.time() - t), 's', flush=True)
    json.dump(res, open(OUT/f'random_{seed}_{ncell}.json', 'w'), indent=1)


def ordered():
    """Small ordered cells (4-Pd conventional cell): the 'small ordered supercell' sampling of the mixing energy."""
    b = json.load(open(OUT/'base.json'))
    a0, EH2 = b['a0'], b['E_H2']
    cell, pos = fcc(a0, 1)
    octa = pos + np.array([0.5, 0, 0]) * a0
    EPd = 4 * b['E_Pd_atom']
    res = []
    for occ in ([0], [0, 1], [0, 1, 2], [0, 1, 2, 3]):
        P = np.vstack([pos, octa[occ]])
        rc = relax(ace, cell, ['Pd'] * 4 + ['H'] * len(occ), P, mode='cell')
        res.append({'x': len(occ) / 4, 'e_cell': (rc['energy'] - EPd - len(occ) * EH2 / 2) / 4})
        print(res[-1], flush=True)
    json.dump(res, open(OUT/'ordered.json', 'w'), indent=1)


if __name__ == '__main__':
    part = sys.argv[1]
    if part == 'base':
        base()
    elif part == 'dilute':
        dilute()
    elif part == 'grids':
        grids()
    elif part == 'random':
        random_part(int(sys.argv[2]), int(sys.argv[3]), [int(v) for v in sys.argv[4].split(',')], int(sys.argv[5]))
    elif part == 'ordered':
        ordered()
