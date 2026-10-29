"""Reference solution for ge-pressure-polymorph-pathways; runs as the Harbor oracle and in check_task.

Route (no literature values): buy the archived relaxation series of each phase; keep near-equilibrium frames
(hydrostatic, symmetric stress; rattled training frames have shear stress and 1-3 eV/A forces); notice that the
compressed ST12/BC8/R8 frames carry residual forces (the ISIF = 6 runs relax the cell only, so the internal
coordinates stay at their zero-pressure values) and correct each such frame's DFT energy by the released
potential's fixed-cell ionic relaxation energy; fit third-order Birch-Murnaghan E(V) per phase; equate enthalpies.
The printed stresses of intermediate relaxation steps are not used (they scatter by up to ~4 GPa).
"""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, curve_fit

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
EV_A3_GPA = 160.21766208
PHASES = ['CD', 'bSn', 'ST12', 'BC8', 'R8']
METASTABLE = ['ST12', 'BC8', 'R8']
PLAN = {'CD': [35, 36, 37, 38, 98], 'bSn': [25, 26, 28, 30, 61], 'ST12': [1, 48, 49, 50, 51],
        'BC8': [30, 31, 32, 33, 61], 'R8': [29, 30, 31, 32]}


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(req):
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    return r


def bm3(V, E0, V0, B0, Bp):
    x = (V0 / V) ** (2 / 3)
    return E0 + 9 * V0 * B0 / 16 * ((x - 1) ** 3 * Bp + (x - 1) ** 2 * (6 - 4 * x))


def p_bm3(V, E0, V0, B0, Bp):
    x = (V0 / V) ** (2 / 3)
    return 1.5 * B0 * (x ** 3.5 - x ** 2.5) * (1 + 0.75 * (Bp - 4) * (x - 1))


def fit(V, E):
    V, E = np.asarray(V, float), np.asarray(E, float)
    i = int(np.argmin(E))
    return curve_fit(bm3, V, E, p0=[E[i], V[i], 0.4, 4.5], maxfev=20000)[0]


def enthalpy(p, P):
    Pe = P / EV_A3_GPA
    V = brentq(lambda v: p_bm3(v, *p) - Pe, 0.6 * p[1], 1.4 * p[1])
    return bm3(V, *p) + Pe * V


def crossing(pa, pb, lo=-2.0, hi=20.0):
    def g(P):
        try:
            return enthalpy(pa, P) - enthalpy(pb, P)
        except ValueError:
            return np.nan
    Ps = np.linspace(lo, hi, 221)
    vals = [g(P) for P in Ps]
    for i in range(len(Ps) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            return float(brentq(g, Ps[i], Ps[i + 1]))
    return float('nan')


def main():
    evidence, frames = [], {g: [] for g in PHASES}
    for g, cfgs in PLAN.items():
        for c in cfgs:
            r = evaluate({'kind': 'record', 'id': '%s-%03d' % (g, c)})
            evidence.append(r['id'])
            res = r['result']
            n = res['natoms']
            for st in res['ionic_steps']:
                s = np.array(st['stress_kB']) / 10
                p = s[:3].mean()
                F = np.linalg.norm(np.array(st['forces_eV_A']), axis=1).max()
                if np.abs(s[:3] - p).max() < 0.5 and np.abs(s[3:]).max() < 0.05 and F < 0.6:
                    frames[g].append({'V': st['volume_A3'] / n, 'E': st['energy_sigma0_eV'] / n, 'n': n, 'F': F,
                                      'cell': st['cell_A'], 'positions': st['positions_A']})
    # internal-coordinate relaxation of the metastable frames with the released potential
    todo = [(g, k) for g in METASTABLE for k in range(len(frames[g]))]
    for i in range(0, len(todo), 20):
        chunk = todo[i:i + 20]
        r = evaluate({'kind': 'dp', 'structures': [{'cell': frames[g][k]['cell'], 'positions': frames[g][k]['positions'],
                                                    'mode': 'relax_ions'} for g, k in chunk]})
        evidence.append(r['id'])
        for (g, k), out in zip(chunk, r['result']['results']):
            frames[g][k]['dE'] = (out['energy_eV'] - out['energy_initial_eV']) / frames[g][k]['n']
    fits = {g: fit([o['V'] for o in frames[g]], [o['E'] + o.get('dE', 0.0) for o in frames[g]]) for g in PHASES}
    raw = {g: fit([o['V'] for o in frames[g]], [o['E'] for o in frames[g]]) for g in METASTABLE}
    tp = {'CD_to_bSn': crossing(fits['CD'], fits['bSn']), 'bSn_to_ST12': crossing(fits['bSn'], fits['ST12']),
          'bSn_to_BC8': crossing(fits['bSn'], fits['BC8'])}
    de = {'ST12_minus_BC8': (fits['ST12'][0] - fits['BC8'][0]) * 1e3, 'R8_minus_BC8': (fits['R8'][0] - fits['BC8'][0]) * 1e3}
    order = sorted(METASTABLE, key=lambda g: fits[g][0])
    submission = {'transition_pressures_GPa': tp, 'zero_pressure_energy_differences_meV_per_atom': de,
                  'metastable_order_at_0GPa': order, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    raw_tp = {g: crossing(fits['bSn'], raw[g]) for g in METASTABLE}
    r8 = crossing(fits['bSn'], fits['R8'])
    (APP/'report.md').write_text(f"""# Ge under pressure (reference solution)

Data: the archived ISIF = 3 (0 GPa) and ISIF = 6 (PSTRESS 1, 3, 6, 10 GPa) relaxation series of CD, beta-Sn,
ST12, BC8 and R8 (PBE, VASP). Frames kept: hydrostatic (deviatoric stress < 0.5 GPa, shear < 0.05 GPa),
max force < 0.6 eV/A. Most single-point records are rattled training configurations (forces 1-3 eV/A, shear
stress) and are not equilibrium energies. Several runs never reached their PSTRESS target (CD at 10 GPa stops
at ~5.6 GPa), and printed stresses of intermediate steps scatter by up to ~4 GPa from -dE/dV, so pressures
come from E(V) fits (third-order Birch-Murnaghan), not from the printed stress or from PSTRESS.

ISIF = 6 relaxes the cell only. CD and beta-Sn have no free internal coordinates, but the compressed ST12, BC8
and R8 frames carry forces up to 0.5 eV/A: their energies are upper bounds. I relaxed the ions of every such
frame at fixed cell with the released potential and added that relaxation energy (up to ~15 meV/atom for ST12
at 10 GPa) to the DFT energy. Without it the beta-Sn -> ST12 / BC8 / R8 crossings would be
{raw_tp['ST12']:.2f} / {raw_tp['BC8']:.2f} / {raw_tp['R8']:.2f} GPa; with it {tp['bSn_to_ST12']:.2f} / {tp['bSn_to_BC8']:.2f} / {r8:.2f} GPa.
The potential alone reproduces the DFT pressures to a few tenths of a GPa but overestimates E(ST12) - E(BC8)
by ~2 meV/atom, so it was used only for the relaxation correction.

Results: CD -> beta-Sn equilibrium at {tp['CD_to_bSn']:.2f} GPa (enthalpy crossing). On decompression beta-Sn
becomes higher in enthalpy than ST12 below {tp['bSn_to_ST12']:.2f} GPa, than R8 below {r8:.2f} GPa and than BC8
below {tp['bSn_to_BC8']:.2f} GPa, but CD is lower than all of them below {tp['CD_to_bSn']:.2f} GPa. Thermodynamics alone
therefore predicts reversion to CD; the metastable products require the CD reconstruction to be kinetically
blocked. At 0 GPa the order is {', '.join(order)} (E - E_BC8 = {de['ST12_minus_BC8']:.1f} and {de['R8_minus_BC8']:.1f} meV/atom), with
ST12 becoming lower than BC8 above about 1 GPa, so which metastable phase forms is decided by kinetics and
decompression rate, which these static data cannot resolve (no barriers are graded or computed here). PBE
underbinds Ge and the CD -> beta-Sn pressure (~8.8 GPa) lies below the experimental ~10-11 GPa; the PBE
pressures should not be compared with experiment without that caveat.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
