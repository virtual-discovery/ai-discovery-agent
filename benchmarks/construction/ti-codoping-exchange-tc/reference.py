"""Reference Monte Carlo for ti-codoping-exchange-tc, with the lab's own code path (environment/lab/dilute.py and
mc_dilute.c), larger supercells, more disorder samples and longer runs than the lab service allows.

Run inside the runtime image (numpy 1.21 / scipy 1.8), from this folder:
    python3 -B reference.py run TAG --comp V:0.05,Cr:0.05 --sizes 3 4 5 6 --samples 64 48 32 16 \
        --temps 20 18 16 15 14 13 12 10 [--neq 1500 --nmeas 6000] [--shells 0-16] [--coco-zero] [--double]
    python3 -B reference.py summary TAG [TAG ...]
Raw results go to $TI_CACHE/ref/TAG_L{size}.json (one file per size).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]/'tasks/ti-codoping-exchange-tc/environment/lab'))
os.environ.setdefault('MC_BIN', str(Path(os.environ.get('TI_CACHE', Path.home()/'.cache/ai-theorist/ti-codoping-exchange-tc'))/'mc_dilute'))
import dilute  # noqa: E402
import model as M  # noqa: E402

OUT = M.CACHE/'ref'


def couplings(names, shells=(0, 16), coco_zero=False, scale=1.0):
    J = {}
    for a in names:
        for b in names:
            j = M.pair_J(a, b)
            if any(x is None for x in j):
                if coco_zero:
                    j = [0.0 if x is None else x for x in j]
                else:
                    raise ValueError(f'missing archived couplings for {a}-{b}')
            J[(a, b)] = [scale * x if shells[0] <= k <= shells[1] else 0.0 for k, x in enumerate(j)]
    return J


def run(a):
    comp = {k: float(v) for k, v in (x.split(':') for x in a.comp.split(','))}
    names = sorted(comp, key=lambda s: M.ZOF[s])
    conc = {n: comp[n] for n in names}
    lo, hi = (int(x) for x in a.shells.split('-'))
    J = couplings(names, (lo, hi), a.coco_zero, 2.0 if a.double else 1.0)
    OUT.mkdir(parents=True, exist_ok=True)
    for size, ns in zip(a.sizes, a.samples):
        path = OUT/f'{a.tag}_L{size}.json'
        if path.exists():
            print('exists', path)
            continue
        seeds = [10**7 + 1000 * size + i for i in range(ns)]
        t0 = time.time()
        res, temps = dilute.run_many(size, conc, J, seeds, a.temps, a.neq, a.nmeas, 2, nworkers=a.workers)
        path.write_text(json.dumps({'tag': a.tag, 'comp': conc, 'size': size, 'temps': temps, 'neq': a.neq,
                                    'nmeas': a.nmeas, 'shells': a.shells, 'coco_zero': a.coco_zero, 'double': a.double,
                                    'seconds': time.time() - t0,
                                    'samples': [{'seed': s, 'n_spins': r['n_spins'], 'n_bonds': r['n_bonds'],
                                                 'data': r['data'].tolist()} for s, r in zip(seeds, res)]}))
        print('wrote', path, round(time.time() - t0), 's', flush=True)


def load(tag):
    out = {}
    for p in sorted(OUT.glob(f'{tag}_L*.json')):
        d = json.loads(p.read_text())
        out[d['size']] = d
    return out


def binder(d, q=0):
    arr = np.array([s['data'] for s in d['samples']])
    m2, m4 = arr[:, :, 2 + 2 * q], arr[:, :, 3 + 2 * q]
    g = m4.mean(0) / m2.mean(0) ** 2
    n = len(arr)
    gj = np.array([(m4.sum(0) - m4[i]) / (n - 1) / ((m2.sum(0) - m2[i]) / (n - 1)) ** 2 for i in range(n)])
    return np.array(d['temps']), g, np.sqrt((n - 1) / n * ((gj - gj.mean(0)) ** 2).sum(0)), m2.mean(0), gj


def crossing(T, ga, gb, tmin=0.0):
    """Highest-temperature crossing of two Binder curves (ga: smaller size). None if the larger size never drops below."""
    order = np.argsort(T)[::-1]
    T, d = T[order], (ga - gb)[order]
    for i in range(len(T) - 1):
        if d[i] <= 0 < d[i + 1] and T[i + 1] >= tmin:
            return float(T[i] + (T[i + 1] - T[i]) * d[i] / (d[i] - d[i + 1]))
    return None


ETA = 0.0375   # 3D Heisenberg anomalous dimension: <m^2> ~ L^-(1+eta) at T_C (dilution is irrelevant, Harris)


def scaled_m2(d, q=0):
    """<m^2> L^(1+eta) per sample (box length proportional to L), for the scaling crossing."""
    arr = np.array([s['data'] for s in d['samples']])
    return arr[:, :, 2 + 2 * q] * d['size'] ** (1 + ETA)


def scaled_crossings(tag, q=0):
    runs = load(tag)
    sizes = sorted(runs)
    out = {}
    for a, b in zip(sizes[:-1], sizes[1:]):
        T = np.array(runs[a]['temps'])
        xa, xb = scaled_m2(runs[a], q), scaled_m2(runs[b], q)
        # ordered side: the larger size has the larger scaled <m^2>; "crossing" uses smaller - larger
        x = crossing(T, xb.mean(0), xa.mean(0))
        na, nb = len(xa), len(xb)
        ja = [crossing(T, xb.mean(0), (xa.sum(0) - xa[i]) / (na - 1)) for i in range(na)]
        jb = [crossing(T, (xb.sum(0) - xb[i]) / (nb - 1), xa.mean(0)) for i in range(nb)]
        def jk(xs):
            xs = [v for v in xs if v is not None]
            n = len(xs)
            return float(np.sqrt((n - 1) / n * np.sum((np.array(xs) - np.mean(xs)) ** 2))) if n > 1 else None
        e1, e2 = jk(ja), jk(jb)
        out[f'{a}/{b}'] = {'T': x, 'err': None if (e1 is None or e2 is None) else float(np.hypot(e1, e2))}
    return out


def crossings(tag, q=0):
    runs = load(tag)
    sizes = sorted(runs)
    out = {}
    for a, b in zip(sizes[:-1], sizes[1:]):
        Ta, ga, ea, _, gja = binder(runs[a], q)
        Tb, gb, eb, _, gjb = binder(runs[b], q)
        assert np.allclose(Ta, Tb)
        x = crossing(Ta, ga, gb)
        # jackknife over the larger set's samples (and the smaller set's), combined in quadrature
        xs_a = [crossing(Ta, gj, gb) for gj in gja]
        xs_b = [crossing(Ta, ga, gj) for gj in gjb]
        def jk(xs):
            xs = [v for v in xs if v is not None]
            n = len(xs)
            return float(np.sqrt((n - 1) / n * np.sum((np.array(xs) - np.mean(xs)) ** 2))) if n > 1 else None
        ea_, eb_ = jk(xs_a), jk(xs_b)
        err = None if (ea_ is None or eb_ is None) else float(np.hypot(ea_, eb_))
        out[f'{a}/{b}'] = {'T': x, 'err': err}
    return out


def summary(tags):
    for tag in tags:
        runs = load(tag)
        print('==', tag, {L: (len(d['samples']), round(np.mean([s['n_spins'] for s in d['samples']])), round(d['seconds']))
                         for L, d in runs.items()})
        sizes = sorted(runs)
        T = np.array(runs[sizes[0]]['temps'])
        rows = {L: binder(runs[L]) for L in sizes}
        for i, t in enumerate(T):
            print(f'{t:6.2f} ' + '  '.join(f'L{L} g={rows[L][1][i]:.4f}({rows[L][2][i]:.4f}) m2={rows[L][3][i]:.3f}' for L in sizes))
        print('binder crossings', json.dumps(crossings(tag)))
        print('scaled <m^2> L^(1+eta) crossings', json.dumps(scaled_crossings(tag)))
        for L in sizes:
            print(f'  L{L} m2*L^(1+eta):', ' '.join(f'{v:.3f}' for v in scaled_m2(runs[L]).mean(0)))


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('tag')
    r.add_argument('--comp', required=True)
    r.add_argument('--sizes', type=int, nargs='+', required=True)
    r.add_argument('--samples', type=int, nargs='+', required=True)
    r.add_argument('--temps', type=float, nargs='+', required=True)
    r.add_argument('--neq', type=int, default=1500)
    r.add_argument('--nmeas', type=int, default=6000)
    r.add_argument('--shells', default='0-16')
    r.add_argument('--coco-zero', action='store_true')
    r.add_argument('--double', action='store_true')
    r.add_argument('--workers', type=int, default=4)
    s = sub.add_parser('summary')
    s.add_argument('tags', nargs='+')
    a = p.parse_args()
    if a.cmd == 'run':
        run(a)
    else:
        summary(a.tags)


if __name__ == '__main__':
    main()
