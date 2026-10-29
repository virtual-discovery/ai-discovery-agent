"""Hidden references, oracle-route values, robustness variants and ablations (called by make_task.py).

Everything is computed through the lab's own request handler (task_lab._run), so the references use exactly
the code path the agent pays for; hidden ingredients (fourth-order rates at 40/80 GPa and 60 GPa/600 K, the
80 GPa phonon-electron table) enter as extra diagonal rates, which is what the lab does internally for
released ones.
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import make_task as mt     # noqa: E402
import tanbte as tb        # noqa: E402

sys.path.insert(0, str(mt.LAB))
import task_lab            # noqa: E402

TARGETS = ['kappa_xx_80GPa_300K_W_mK', 'kappa_zz_80GPa_300K_W_mK', 'kappa_xx_40GPa_300K_W_mK',
           'kappa_xx_60GPa_600K_W_mK', 'phel_resistivity_fraction_80GPa_300K']


def kap(p, T, channels=('3ph', 'isotope'), solver='iterative', extra=None, sb=0.1, mesh=16):
    r = {'kind': 'kappa', 'pressure_GPa': p, 'mesh': mesh, 'T_K': T, 'scalebroad': sb, 'channels': list(channels),
         'solver': solver}
    if extra is not None:
        r['extra_rates_ps'] = np.asarray(extra).tolist()
    task_lab.validate(r)
    out = task_lab._run(r)
    if solver == 'iterative' and not out['converged']:
        print('WARNING: not converged', p, T, channels, solver)
    return np.array(out['kappa_W_mK'])


def f4(p, T):
    d = mt.fourph_rates(p, T)
    return d['pp'] + d['pm'] + d['mm']


def ratio(a, b):
    return np.where(b > 0, a / np.where(b > 0, b, 1), 1.0)


def summarize(K80, K80n, K40, K60h, name, store):
    store[name] = {'kappa_xx_80GPa_300K_W_mK': K80[0, 0], 'kappa_zz_80GPa_300K_W_mK': K80[2, 2],
                   'kappa_xx_40GPa_300K_W_mK': K40[0, 0], 'kappa_xx_60GPa_600K_W_mK': K60h[0, 0],
                   'phel_resistivity_fraction_80GPa_300K': 1 - K80[0, 0] / K80n[0, 0]}
    print(name, {k: round(float(v), 4) for k, v in store[name].items()}, flush=True)


def main():
    g = task_lab.grid(0, 16)
    phel = {p: tb.read_phel(mt.RAW/str(p)/'BTE.sr_phel_metal', g)[1] for p in mt.PRESSURES}
    raw = {p: np.loadtxt(mt.RAW/str(p)/'BTE.sr_phel_metal')[:, 1] for p in mt.PRESSURES}
    F = {(p, T): f4(p, T) for (p, T) in mt.FOURPH_RUNS}
    S = {}
    # ---- reference (hidden ingredients)
    K80 = kap(80, 300, extra=F[(80, 300)] + phel[80])
    K80n = kap(80, 300, extra=F[(80, 300)])
    K40 = kap(40, 300, extra=F[(40, 300)] + phel[40])
    K60h = kap(60, 600, extra=F[(60, 600)] + phel[60])
    summarize(K80, K80n, K40, K60h, 'reference', S)
    # MC noise: single seeds at 80 GPa
    for s in (1, 2):
        d = np.load(mt.FOURPH/('r4_80_16_0.1_8192_%d_300.npz' % s))
        f = d['pp'] + d['pm'] + d['mm']
        k = kap(80, 300, extra=f + phel[80])
        S['mc_seed%d_kappa_xx_80' % s] = k[0, 0]
        print('seed', s, k[0, 0], k[2, 2])

    # ---- oracle route: mode-wise linear interpolation of rates in pressure; T-ratio of 4ph from 20/100 GPa
    def route(f80, p80, f40, f60h, extra_scale=1.0, solver='iterative', sb=0.1, name=None, ch=('3ph', 'isotope')):
        a = kap(80, 300, ch, solver, f80 + p80 * extra_scale, sb)
        b = kap(80, 300, ch, solver, f80, sb)
        c = kap(40, 300, ch, solver, f40 + phel[40] * extra_scale, sb)
        d = kap(60, 600, ch, solver, f60h + phel[60] * extra_scale, sb)
        summarize(a, b, c, d, name, S)

    lin80 = 0.5 * (F[(60, 300)] + F[(100, 300)])
    lin40 = 0.5 * (F[(20, 300)] + F[(60, 300)])
    pel80 = 0.5 * (phel[60] + phel[100])
    rT = 0.5 * (ratio(F[(20, 600)], F[(20, 300)]) + ratio(F[(100, 600)], F[(100, 300)]))
    f60h = F[(60, 300)] * rT
    route(lin80, pel80, lin40, f60h, name='oracle')
    # robustness
    geo = lambda a, b: np.sqrt(a * b)                                       # noqa: E731
    route(geo(F[(60, 300)], F[(100, 300)]), geo(phel[60], phel[100]), geo(F[(20, 300)], F[(60, 300)]),
          F[(60, 300)] * np.sqrt(ratio(F[(20, 600)], F[(20, 300)]) * ratio(F[(100, 600)], F[(100, 300)])),
          name='robust_geometric_interpolation')
    # quadratic in P through three released pressures, mode by mode
    def quad(ps, ys, x):
        c = np.polyfit(ps, np.stack([y.ravel() for y in ys]), 2)
        return np.clip((c[0] * x ** 2 + c[1] * x + c[2]).reshape(ys[0].shape), 0, None)
    route(quad([20, 60, 100], [F[(20, 300)], F[(60, 300)], F[(100, 300)]], 80),
          quad([40, 60, 100], [phel[40], phel[60], phel[100]], 80),
          quad([0, 20, 60], [F[(0, 300)], F[(20, 300)], F[(60, 300)]], 40),
          F[(60, 300)] * ratio(F[(100, 600)], F[(100, 300)]), name='robust_quadratic_and_T_ratio_from_100GPa')
    route(lin80, pel80, lin40, F[(60, 300)] * ratio(F[(20, 600)], F[(20, 300)]), name='robust_T_ratio_from_20GPa')
    # ---- ablations
    route(lin80, pel80, lin40, f60h, solver='RTA', name='abl_RTA')
    z = np.zeros_like(lin80)
    route(z, pel80, z, z, name='abl_no_4ph')
    route(lin80, z, lin40, f60h, name='abl_no_phel')
    route(lin80, pel80, lin40, f60h, extra_scale=2 * np.pi, name='abl_phel_times_2pi')
    # archived rows read in the lab's own (i-fastest) order instead of ShengBTE's (k-fastest) order
    def naive(p):
        return raw[p].reshape(6, -1).T
    a = kap(80, 300, extra=lin80 + 0.5 * (naive(60) + naive(100)))
    b = kap(80, 300, extra=lin80)
    c = kap(40, 300, extra=lin40 + naive(40))
    d = kap(60, 600, extra=f60h + naive(60))
    summarize(a, b, c, d, 'abl_phel_row_order', S)
    route(lin80, pel80, lin40, F[(60, 300)], name='abl_4ph_not_rescaled_to_600K')
    route(lin80, pel80, lin40, F[(60, 300)] * 4.0, name='abl_4ph_classical_T2_scaling')
    route(lin80, pel80, lin40, f60h, sb=1.0, name='abl_scalebroad_1')
    # kappa-level interpolation of released full-physics conductivities
    full = {p: kap(p, 300, ('3ph', 'isotope', '4ph', 'phel')) for p in (20, 60, 100)}
    noel = {p: kap(p, 300, ('3ph', 'isotope', '4ph')) for p in (60, 100)}
    k80 = 0.5 * (full[60] + full[100])
    k80n = 0.5 * (noel[60] + noel[100])
    k40 = 0.5 * (full[20] + full[60])
    k60h = full[60] * 0.5                       # kappa ~ 1/T
    summarize(k80, k80n, k40, k60h, 'abl_kappa_level_interpolation', S)
    # anisotropy: reporting the trace average for the xx/zz fields
    ref = S['reference']
    S['abl_trace_average'] = dict(ref, **{
        'kappa_xx_80GPa_300K_W_mK': np.trace(K80) / 3, 'kappa_zz_80GPa_300K_W_mK': np.trace(K80) / 3,
        'kappa_xx_40GPa_300K_W_mK': np.trace(K40) / 3, 'kappa_xx_60GPa_600K_W_mK': np.trace(K60h) / 3})
    # extra context: released full-physics values and 3ph-only at both broadenings (report, not graded)
    ctx = {}
    for p in (0, 20, 60, 100):
        ctx['full_%d' % p] = np.diag(kap(p, 300, ('3ph', 'isotope', '4ph', 'phel'))).tolist()
    for p in (60, 80):
        ctx['3ph_sb1_%d' % p] = np.diag(kap(p, 300, ('3ph', 'isotope'), sb=1.0)).tolist()
    S['context'] = ctx
    print(json.dumps(ctx))
    (HERE/'calibration_raw.json').write_text(json.dumps(S, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
