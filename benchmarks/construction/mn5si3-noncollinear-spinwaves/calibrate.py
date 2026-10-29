"""Calibration of mn5si3-noncollinear-spinwaves: the reference route, defensible variants and wrong paths.

    python3 -B calibrate.py            # after make_task.py; writes solution/ablation_values.json and calibration_table.json

All routes start from data the agent can buy (task_lab.run), never from the hidden targets or paper parameters.
"""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mn5si3-noncollinear-spinwaves'
sys.path.insert(0, str(TASK/'solution'))
sys.path.insert(0, str(TASK/'environment/lab'))
import task_lab          # noqa: E402
import solve             # noqa: E402
import mnsi_lswt as L    # noqa: E402

CATALOG = json.loads((TASK/'environment/agent/catalog.json').read_text())
HID = json.loads((TASK/'environment/lab/private/hidden.json').read_text())


def buy(req):
    return task_lab.run(req, None)


def all_data(cuts=solve.PLAN['cuts']):
    return {'exchange': {n: buy({'kind': 'exchange_set', 'name': n}) for n in task_lab.DATA['exchange']},
            'dispersion': {k: buy({'kind': 'dispersion_table', 'direction': k}) for k in 'hk'},
            'cuts': [buy({'kind': 'constant_energy_cut', 'direction': k, 'energy_meV': E}) for k, E in cuts],
            'gap_table': buy({'kind': 'gap_table'})}


def l_rms(pred):
    dev = []
    for E, t in HID['l_targets'].items():
        lo, hi = sorted(pred[E])
        dev += [lo - t['low']['q_rlu'], hi - t['high']['q_rlu']]
    return float(np.sqrt(np.mean(np.square(dev))))


def l_pred_from_model(m, scale, sq, weights=True, lowest=False):
    f = L.DIRS['l'][0]
    if lowest:
        spec = []
        for x in L.line_grid('l'):
            w, i = L.sqw(np.array(f(x), float), *solve.as_tuple(m))
            spec.append((w[:1], np.ones(min(1, len(w)))))
    else:
        spec = L.line_spectrum('l', solve.as_tuple(m), weights=weights)
    energies = CATALOG['to_predict']['l_dispersion']['energies_meV']
    x = np.linspace(0.5, 1.5, 401)
    guesses = []
    for E in energies:
        prof = L.cut_profile(spec, 'l', E, x, solve.SIGMA_E, sq, scale)
        guesses.append(max(0.01, abs(x[x < 1][np.argmax(prof[x < 1])] - 1)))
    pk = L.predicted_peaks(spec, 'l', energies, guesses, solve.SIGMA_E, sq, scale)
    return {('%g' % E): [round(float(pk[E][0]), 4), round(float(pk[E][1]), 4)] for E in energies}


def iso_cone(data):
    pts = []
    for key, c, ax in (('h', 2.0, 0), ('k', 3.0, 1)):
        pts += [(abs(p['q_rlu'] - c) * 2 * np.pi / L.LAT[ax], p['E_meV']) for p in data['dispersion'][key]['points']
                if p['q_err_rlu'] > 0.001]
    q, E = np.array(pts).T
    c = np.linalg.lstsq(np.c_[np.ones_like(q), q ** 2], E ** 2, rcond=None)[0]
    out = {}
    for Et in CATALOG['to_predict']['l_dispersion']['energies_meV']:
        dq = np.sqrt(max(Et ** 2 - c[0], 0) / c[1]) * L.LAT[2] / (2 * np.pi)
        out['%g' % Et] = [round(1 - dq, 4), round(1 + dq, 4)]
    return out


def verdict(row):
    th = HID['thresholds']
    ok = {'l': row['l_rms'] <= th['l_rms_rlu'],
          'optical': th['optical_window_meV'][0] <= row['optical'] <= th['optical_window_meV'][1],
          'mn1': row['mn1'] is True,
          'canting': th['canting_window_deg'][0] <= row['canting'] <= th['canting_window_deg'][1]}
    return ok


def main():
    data = all_data()
    rows = {}

    def record(name, r, **over):
        row = {'pred': r['pred'], 'optical': r['optical'], 'mn1': r['mn1'], 'canting': r['canting']}
        row.update(over)
        row['l_rms'] = l_rms(row['pred'])
        row['passes'] = verdict(row)
        row['all_pass'] = all(row['passes'].values())
        for k in ('scale', 'sigQ', 'K_y'):
            if k in r:
                row[k] = r[k]
        rows[name] = row
        print('%-40s l_rms %.4f optical %6.2f mn1 %-5s canting %6.2f  pass %s %s' % (
            name, row['l_rms'], row['optical'], row['mn1'], row['canting'], row['all_pass'],
            {k: v for k, v in row['passes'].items() if not v}), flush=True)

    ref = solve.analyse(CATALOG, data)
    record('reference', ref)
    print('   scale scan', {k: round(v, 3) for k, v in ref['scale_scan'].items()}, 'sigQ', ref['sigQ'], 'K_y', ref['K_y'])
    # ---------------- defensible variants (must pass)
    record('robust_collinear_reference_table', solve.analyse(CATALOG, data, {'exchange': 'collinear_AFMz_r4.5'}))
    record('robust_long_range_r17_table', solve.analyse(CATALOG, data, {'exchange': 'noncollinear_r17'}))
    record('robust_sigQ_0.04', solve.analyse(CATALOG, data, {'sigQ': 0.04}))
    record('robust_sigQ_0.08', solve.analyse(CATALOG, data, {'sigQ': 0.08}))
    record('robust_scale_0.9', solve.analyse(CATALOG, data, {'scale': 0.9}))
    record('robust_scale_1.1', solve.analyse(CATALOG, data, {'scale': 1.1}))
    record('robust_all_cuts_for_resolution', solve.analyse(CATALOG, all_data([(c['direction'], c['energy_meV'])
                                                                               for c in task_lab.DATA['cuts'].values()])))
    tables = {'mn2': data['exchange']['noncollinear_r4.5']['bonds'], 'mn1': data['exchange']['mn1_couplings']['bonds']}
    aniso = [(ref['K_y'], [0, 1, 0])] if ref['K_y'] > 0 else []
    both = solve.build(CATALOG, tables, True, aniso)
    small = dict(both, S=np.r_[np.full(8, solve.S_SPIN), np.full(8, 0.6)])
    record('robust_mn1_spin_0.6', {'pred': l_pred_from_model(small, ref['scale'], ref['sigQ']), 'mn1': True,
                                   'canting': solve.canting_deg(both), 'optical': solve.zone_centre_mode(small, ref['scale'])})
    # ---------------- informational
    record('info_isotropic_cone_from_h_k', dict(ref, pred=iso_cone(data)))
    record('info_no_anisotropy', solve.analyse(CATALOG, data, {'aniso': 'none'}))
    # ---------------- wrong paths (must fail)
    record('dp1_mn2_only_sites', solve.analyse(CATALOG, data, {'mn1': False}))
    record('dp2_bonds_double_counted', solve.analyse(CATALOG, data, {'scale': 2.0}))
    record('dp3_free_ion_spin_5_2', solve.analyse(CATALOG, data, {'scale': solve.S_SPIN / 2.5}))
    record('dp4_lowest_branch_only', dict(ref, pred=l_pred_from_model(both, ref['scale'], ref['sigQ'], lowest=True)))
    record('dp4_branch_weights_ignored', solve.analyse(CATALOG, data, {'weights': False}))
    record('dp5_collinear_ground_state_not_relaxed', ref, canting=0.0)
    record('dp6_angle_to_mn1_axis', ref, canting=90.0 - ref['canting'])
    (HERE/'calibration_table.json').write_text(json.dumps(rows, indent=1) + '\n')
    abl = {k: {'pred': v['pred'], 'optical': v['optical'], 'mn1': v['mn1'], 'canting': v['canting']} for k, v in rows.items()}
    (TASK/'solution/ablation_values.json').write_text(json.dumps(abl, indent=1) + '\n')


def rescore():
    rows = json.loads((HERE/'calibration_table.json').read_text())
    for name, row in rows.items():
        row['passes'] = verdict(row)
        row['all_pass'] = all(row['passes'].values())
        print('%-40s l_rms %.4f optical %6.2f canting %6.2f pass %s' % (name, row['l_rms'], row['optical'], row['canting'],
                                                                        row['all_pass']))
    (HERE/'calibration_table.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    rescore() if sys.argv[1:] == ['rescore'] else main()
