#!/usr/bin/env python3
"""Linear (two-branch) quasiparticle corrections fitted on a BerkeleyGW eqp1.dat - the eigenvalue-update step of the evGW ladder.
Extracted verbatim (logic) from runs/116_lisoc2/lisoc2_chain.py:read_eqp/fit_scissors/scissor_lines (SOC) and
runs/138_complete/nosoc_chain.py:fit_scissors_1spin (no-SOC, spin 1).
  valence   n <= VBM : E_QP = E_MF + evs + evdel (E_MF - ev0),  ev0 = max E_MF over the fit window (VBM)
  conduction n >= CBM: E_QP = E_MF + ecs + ecdel (E_MF - ec0),  ec0 = min E_MF over the fit window (CBM)
Fit windows: SOC valence 145-148, conduction 149-152 (VBM 148 / CBM 149); no-SOC spin 1 valence 71-74, conduction 75-78 (VBM 74 / CBM 75).
Usage: fit_scissors.py soc|nosoc eqp1.dat [--lines] [--compare epsilon.inp|sigma.inp|absorption.inp]
  --lines    print the 6 BerkeleyGW keywords (evs ev0 evdel ecs ec0 ecdel) for the next epsilon.inp / sigma.inp / absorption.inp
  --compare  print the difference to the keywords present in an existing input (the archived inputs carry the original fits)"""
import json, sys
import numpy as np
W = {'soc': dict(fitv=(145, 148), fitc=(149, 152), vbm=148, cbm=149), 'nosoc': dict(fitv=(71, 74), fitc=(75, 78), vbm=74, cbm=75)}
def read_eqp(path, spin=1):
    out, k = {}, None
    for line in open(path):
        p = line.split()
        if len(p) == 4 and '.' in p[0]: k = tuple(round(float(x), 9) for x in p[:3]); out[k] = {}
        elif len(p) >= 4 and k is not None and int(p[0]) == spin: out[k][int(p[1])] = (float(p[2]), float(p[3]))
    return out
def fit(path, chain):
    w = W[chain]; eqp = read_eqp(path)
    val = [(mf, qp) for k in eqp for b, (mf, qp) in eqp[k].items() if w['fitv'][0] <= b <= w['fitv'][1]]
    con = [(mf, qp) for k in eqp for b, (mf, qp) in eqp[k].items() if w['fitc'][0] <= b <= w['fitc'][1]]
    ev0 = max(mf for mf, _ in val); ec0 = min(mf for mf, _ in con)
    def linfit(rows, x0):
        x = np.array([mf - x0 for mf, _ in rows]); y = np.array([qp - mf for mf, qp in rows]); A = np.stack([np.ones_like(x), x], axis=1)
        (s, d), *_ = np.linalg.lstsq(A, y, rcond=None); return float(s), float(d), float(np.sqrt(np.mean((A @ np.array([s, d]) - y) ** 2)))
    evs, evdel, rv = linfit(val, ev0); ecs, ecdel, rc = linfit(con, ec0)
    gap = min(eqp[k][w['cbm']][1] for k in eqp) - max(eqp[k][w['vbm']][1] for k in eqp); dgap = min(eqp[k][w['cbm']][1] - eqp[k][w['vbm']][1] for k in eqp)
    return dict(evs=evs, ev0=float(ev0), evdel=evdel, ecs=ecs, ec0=float(ec0), ecdel=ecdel, rms_v=rv, rms_c=rc, gap_qp=gap, gap_direct=dgap)
def lines(sc): return ''.join('%s %.6f\n' % (k, sc[k]) for k in ('evs', 'ev0', 'evdel', 'ecs', 'ec0', 'ecdel'))
if __name__ == '__main__':
    chain, path = sys.argv[1], sys.argv[2]; sc = fit(path, chain)
    if '--lines' in sys.argv: sys.stdout.write(lines(sc))
    else: print(json.dumps({k: round(v, 6) for k, v in sc.items()}))
    if '--compare' in sys.argv:
        f = sys.argv[sys.argv.index('--compare') + 1]; old = {l.split()[0]: float(l.split()[1]) for l in open(f) if l.split() and l.split()[0] in ('evs', 'ev0', 'evdel', 'ecs', 'ec0', 'ecdel')}
        print('vs %s: ' % f + ', '.join('%s %+.2e' % (k, round(sc[k], 6) - old[k]) for k in old))
