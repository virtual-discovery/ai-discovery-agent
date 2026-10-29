"""Fit the model families to the released data and record the wrong-path / defensible predictions of the withheld data.

    CEIN3_RAW=<raw dir> python3 -B calibrate.py        (run inside ai-theorist-benchmark-runtime:20261007)
Writes calibration_fits.json: per variant the parameters, velocity/bandwidth, the scores against the withheld data
(RM map, RM 1.4-meV cut, Q-integrated chi'') and the predictions themselves (used by make_task.py).
"""
import json
import time
import numpy as np
import cmodel as M
from cfit import L, he_ds, hr_ds, cut_ds, integrated_ds, fit, predict, chi2r, score_all, project

HERE = __import__('pathlib').Path(__file__).resolve().parent
QRM, ERM = L('HighResolution_Q_RM.dat'), L('HighResolution_E_RM.dat')


def withheld(hp=0.08, int_scale=1.0):
    return hr_ds('RM', half_perp=hp), cut_ds('RM', '1p4meV', hp), integrated_ds(int_scale)


def visible(hp=0.08, perp_scale=1.0, he=True, hr=True):
    out = [he_ds(p, perp_scale=perp_scale) for p in ['RG', 'GX', 'XM', 'MG']] if he else []
    if hr:
        out += [hr_ds(p, half_perp=hp) for p in ['RG', 'RX']]
        out += [cut_ds(p, t, hp) for p, t in [('RG', '0p6meV'), ('RG', '1p4meV'), ('RX', '0p6meV'), ('RX', '1p4meV')]]
    return out


def evaluate(p, hp=0.08, int_scale=1.0):
    rm, cut, integ = withheld(hp, int_scale)
    full = M.HRRows('RM', QRM, ERM, half_perp=hp).model(p)
    pc, pi = cut.obj.model(p), integ.obj.model(p)
    pr = M.props(p)
    return {'params': p, 'props': pr,
            'scores': {'rm_map': round(chi2r(full[rm.mask], rm.y, rm.e), 4), 'rm_cut_1p4meV': round(chi2r(pc, cut.y, cut.e), 4),
                       'integrated': round(chi2r(pi, integ.y, integ.e, offset=False), 4)},
            'pred': {'rm_map': np.round(full, 6).tolist(), 'rm_cut_1p4meV': np.round(pc, 6).tolist(), 'integrated': np.round(pi, 6).tolist()}}


LR0 = {'SJ1': 0.3, 'A': 1.0, 'kappa': 0.06, 'gamma_hr': 0.3, 'sigma_he': 0.4}
VARIANTS = [
    # name, data selection, N shells, start, fixed, geometry
    ('J1_HEonly', dict(hr=False), 1, {'SJ1': 0.45, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05, 'gamma_hr': 0.05, 'gamma_he': 0.0, 'gap': 0}, {}),
    ('J1_forced_all', {}, 1, {'SJ1': 0.45, 'gamma_hr': 0.2, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05, 'gap': 0}, {}),
    ('J123_short', {}, 3, {'SJ1': 0.4, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05, 'gap': 0}, {}),
    ('J6shell_short', {}, 6, {'SJ1': 0.4, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05, 'gap': 0}, {}),
    ('noperp_J1LR', {}, 1, LR0, {'gap': 0}, {'hp': 0.0, 'perp_scale': 1e-6}),
    ('J1LR_HRonly', dict(he=False), 1, LR0, {'gap': 0, 'sigma_he': 0.3, 'gamma_he': 0.0}, {}),
    ('J1LR', {}, 1, LR0, {'gap': 0}, {}),
    ('J123LR', {}, 3, {'SJ1': 0.2, 'A': 1.5, 'kappa': 0.06, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {'gap': 0}, {}),
    ('J1LR_h004', {}, 1, {'SJ1': 0.3, 'A': 0.8, 'kappa': 0.045, 'gamma_hr': 0.5, 'sigma_he': 0.5}, {'gap': 0}, {'hp': 0.04}),
    ('J10shell', {}, 10, {'SJ1': 0.4, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05, 'gap': 0}, {}),
]


def main():
    out = {}
    for name, sel, N, x0, fixed, geo in VARIANTS:
        t = time.time()
        hp, ps = geo.get('hp', 0.08), geo.get('perp_scale', 1.0)
        dsets = visible(hp=hp, perp_scale=ps, **sel)
        p, cost = fit(dsets, N, x0, fixed, nfev=300)
        ev = evaluate(p, hp=hp, int_scale=ps)
        ev['visible_scores'] = score_all(visible(hp=hp, perp_scale=ps), p)
        ev['fit_cost'] = round(cost, 4)
        out[name] = ev
        pr = ev['props']
        print(f"{name:16s} {time.time()-t:5.0f}s SJ={np.round(p['SJ'],3).tolist()} A={p['A']:.3f} kappa={p['kappa']:.3f} "
              f"g_hr={p['gamma_hr']:.3f} sig_he={p['sigma_he']:.3f} | v={pr['v']:.1f} W={pr['W']:.3f} eta={pr['eta']:.3f} | "
              f"withheld {ev['scores']} | visible {ev['visible_scores']}", flush=True)
    # reference: long-range model fitted to all released data, withheld sets included
    allsets = visible() + [hr_ds('RM'), cut_ds('RM', '0p4meV'), cut_ds('RM', '1p4meV')]
    p, cost = fit(allsets, 3, {'SJ1': 0.2, 'A': 1.5, 'kappa': 0.06, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {'gap': 0}, nfev=300)
    out['reference_all_data'] = evaluate(p)
    pr = out['reference_all_data']['props']
    print(f"reference_all_data v={pr['v']:.1f} W={pr['W']:.3f} eta={pr['eta']:.3f} {out['reference_all_data']['scores']}")
    # constant (non-dispersive continuum) predictions
    rm, cut, integ = withheld()
    out['continuum'] = {'scores': {'rm_map': round(chi2r(np.ones(len(rm.y)) + 1e-9 * np.arange(len(rm.y)), rm.y, rm.e), 4),
                                   'rm_cut_1p4meV': round(chi2r(np.ones(len(cut.y)) + 1e-9 * np.arange(len(cut.y)), cut.y, cut.e), 4),
                                   'integrated': round(chi2r(np.ones(len(integ.y)), integ.y, integ.e, offset=False), 4)}}
    print('continuum', out['continuum']['scores'])
    (HERE/'calibration_fits.json').write_text(json.dumps(out) + '\n')


if __name__ == '__main__':
    main()
