"""Generate the staged, lab, private and ablation data for her-alloy-workfunction-descriptor.

    cd benchmarks/construction/her-alloy-workfunction-descriptor && python3 -B make_task.py

Reads pool.json.gz (written by extract.py). Writes into benchmarks/tasks/<id>/:
environment/agent/{catalog.json,descriptors.json}, environment/lab/data.json,
environment/lab/private/{hidden.json,hidden_values.json}, solution/ablation_values.json,
and prints the calibration numbers (also saved to calibration.json here).
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/her-alloy-workfunction-descriptor'
P = M.POOL
BUDGET = 60
COSTS = {'seccm_row': 1, 'lib3_spot': 2, 'coverage_path': 1, 'dft_composition': 1}
MAX_LIB3 = 5
N_ELIGIBLE = 40
THRESHOLDS = {'lib3_rmse_log10': 0.135, 'volcano_centre_window_eV': [-0.52, -0.36],
              'coverage_shift_window_eV': [0.07, 0.45], 'alpha_wf_window': [-0.06, 0.13]}
# Oracle purchase plan (mirrored by solution/solve.py).
ORACLE_ROWS = [3, 7, 12, 16]
YS = sorted({round(s['y_mm'], 1) for s in P['libraries']['Lib-1']['spots']})


def row_of(s):
    return YS.index(round(s['y_mm'], 1))


def r3(x):
    return float(round(x, 3))


def eligible_and_graded():
    lib3 = P['libraries']['Lib-3']['spots']
    rng = np.random.default_rng(20261026)
    idx = sorted(rng.choice(len(lib3), N_ELIGIBLE, replace=False).tolist())
    elig = [lib3[i]['mdi'] for i in idx]
    graded = [s['mdi'] for s in lib3 if s['mdi'] not in elig]
    return elig, graded


def fps(comp, k, start=2):
    """Farthest-point sampling in composition space starting from the spot richest in element `start`."""
    idx = [int(np.argmax(comp[:, start]))]
    while len(idx) < k:
        d = np.min(np.linalg.norm(comp[:, None, :] - comp[None, idx, :], axis=2), axis=1)
        idx.append(int(np.argmax(d)))
    return idx


def coverage_record(p):
    sites = {}
    path = []
    for n, (s, t, gd, gr) in enumerate(zip(p['site_sequence'], p['theta_ml'], p['dG_diff_ev'], p['dG_rep_ev'])):
        sites[f'hollow{s:02d}'] = round(gd - gr, 6)       # single-H (1/9 ML) energy of that hollow
        path.append({'n_H': n + 1, 'theta_ML': t, 'site_added': f'hollow{s:02d}',
                     'differential_adsorption_energy_eV': gd})
    return {'single_site_adsorption_energies_eV': dict(sorted(sites.items())), 'filling_path': path}


def main():
    elig, graded = eligible_and_graded()
    libs = P['libraries']
    # ---------------- staged (agent-visible) ----------------
    catalog = {
        'system': 'Ag-Au-Pd-Pt thin-film materials libraries (combinatorial co-sputtering), acidic HER, scanning electrochemical cell microscopy (SECCM); DFT: RPBE/GPAW fcc(111) 3x3x4 slabs. Released by Liasi et al. (2026).',
        'budget': BUDGET, 'purchase_costs': COSTS, 'max_lib3_spot_purchases': MAX_LIB3,
        'potentials_mV': [-100, -200, -300, -400],
        'current_convention': 'Current density j in A/cm^2 at the stated potential vs RHE; cathodic (HER) current is negative.',
        'libraries': {}, 'dft_compositions': [], 'coverage_paths': [],
        'descriptors_file': '/app/descriptors.json',
    }
    for name, lib in libs.items():
        spots = []
        for s in lib['spots']:
            row = {'mdi': s['mdi'], 'row': row_of(s), 'x_mm': s['x_mm'], 'y_mm': s['y_mm'],
                   'composition_at_frac': {k: round(v, 3) for k, v in s['composition'].items()}}
            spots.append(row)
        entry = {'archive_id': lib['archive_dir'], 'n_spots': len(spots), 'spots': spots}
        if name == 'Lib-3':
            entry['calibration_spots_available'] = elig
            entry['spots_to_predict'] = graded
            entry['note'] = 'Lib-3 currents are not released except for purchases of calibration_spots_available (lib3_spot).'
        else:
            entry['rows_available'] = sorted({r['row'] for r in spots})
        catalog['libraries'][name] = entry
    adsorption = {}
    for rec in P['dft_adsorption']:
        adsorption.setdefault(rec['comp_id'], []).append(rec)
    wfs = P['dft_work_functions']
    comp_ids = sorted(wfs)
    # map catalog dft entries to comp ids by composition
    catalog['dft_compositions'] = [{'id': cid, 'composition_at_frac': {k: round(v, 4) for k, v in wfs[cid]['composition'].items()},
                                    'n_variants': len(wfs[cid]['variants'])} for cid in comp_ids]
    paths = [p for p in P['coverage_isotherms']['paths'] if p['direction'] == 'adsorption']
    for p in paths:
        catalog['coverage_paths'].append({'id': f"{p['comp_id']}/{p['variant']}", 'composition_id': p['comp_id'],
                                          'decoration': p['variant'], 'n_sites': 9})
    descriptors = {'note': 'Per measurement area (MDI). gnn: DimeNet++ surrogate applied to 1,000 random fcc(111) slabs at the EDX composition; '
                           'hads_histogram counts the 9,000 dilute (1/9 ML) hydrogen adsorption energies (eV, relative to 1/2 E(H2), RPBE, no ZPE or entropy) '
                           'in bins of width bin_width_eV, bin i covering [bin_start_eV + (first_bin + i)*w, +w); wf: mean and spread of the 1,000 slab work functions (eV). '
                           'linear_combination: scalar baseline estimates of the adsorption energy (meV) and work function (eV), the composition-weighted sum of per-element coefficients fitted to the DFT set.',
                   'bin_start_eV': P['bin_edges_start_eV'], 'bin_width_eV': P['bin_width_eV'], 'libraries': {}}
    for name, lib in libs.items():
        d = {}
        for s in lib['spots']:
            g, l = s['dimenet'], s['lincomb']
            d[str(s['mdi'])] = {'gnn': {'hads_histogram': {'first_bin': g['hads']['first_bin'], 'counts': g['hads']['counts']},
                                        'wf_mean_eV': r3(g['wf_mean_eV']), 'wf_std_eV': r3(g['wf_std_eV'])},
                                'linear_combination': {'hads_meV': int(round(1000 * l['hads']['mean'])), 'wf_eV': r3(l['wf_mean_eV'])}}
        descriptors['libraries'][name] = d
    # ---------------- lab data ----------------
    data = {'seccm': {}, 'lib3': {}, 'dft': {}, 'coverage': {}}
    for name in ('Lib-1', 'Lib-2'):
        for s in libs[name]['spots']:
            for mv in M.POTS:
                data['seccm'].setdefault(f'{name}|{row_of(s)}|{mv}', []).append({'mdi': s['mdi'], 'j_A_cm2': s['j_A_cm2'][str(mv)]})
    lib3 = {s['mdi']: s for s in libs['Lib-3']['spots']}
    for m in elig:
        data['lib3'][str(m)] = {str(-mv): lib3[m]['j_A_cm2'][str(mv)] for mv in M.POTS}
    for cid in comp_ids:
        variants = {}
        for rec in adsorption.get(cid, []):
            variants.setdefault(rec['variant'], {'sites': {}})['sites'][rec['site']] = rec['e_ads_ev']
        for v, w in wfs[cid]['variants'].items():
            variants.setdefault(v, {'sites': {}})['work_function_eV'] = w['wf_top_ev']
        data['dft'][cid] = {'composition_at_frac': wfs[cid]['composition'], 'variants': variants}
    for p in paths:
        data['coverage'][f"{p['comp_id']}/{p['variant']}"] = coverage_record(p)
    # ---------------- private ----------------
    hidden = {'lib3_targets_A_cm2': {str(m): [lib3[m]['j_A_cm2'][str(mv)] for mv in M.POTS] for m in graded},
              'thresholds': THRESHOLDS}
    # ---------------- reference analysis on the oracle's purchases ----------------
    rows = set(ORACLE_ROWS)
    L = {n: M.Lib(n, {s['mdi'] for s in libs[n]['spots'] if row_of(s) in rows}) for n in ('Lib-1', 'Lib-2')}
    Lf = {n: M.Lib(n) for n in ('Lib-1', 'Lib-2')}
    Ll = {n: M.Lib(n, {s['mdi'] for s in libs[n]['spots'] if row_of(s) in rows}, sur='lincomb') for n in ('Lib-1', 'Lib-2')}
    E3 = M.Lib('Lib-3', set(elig))
    cal_mdi = sorted(E3.mdi[fps(E3.comp, MAX_LIB3)].tolist())
    cal = M.Lib('Lib-3', set(cal_mdi))
    calL = M.Lib('Lib-3', set(cal_mdi), sur='lincomb')
    G = M.Lib('Lib-3', set(graded))
    GL = M.Lib('Lib-3', set(graded), sur='lincomb')
    CW = ('alpha_c', 'alpha_ads', 'e_opt', 'alpha_wf')

    def calibrate(q, c, kw, both=True):
        def resid(x):
            qq = dict(q, logk=x[0], **({'logjlim': x[1]} if both else {}))
            pr = M.predict(c, qq, **kw)
            obs = -c.j
            ok = obs > 1e-9
            return np.log10(pr[ok]) - np.log10(obs[ok])
        x0 = [q['logk'], q['logjlim']] if both else [q['logk']]
        r = least_squares(resid, x0, loss='huber', f_scale=0.3)
        return dict(q, logk=r.x[0], **({'logjlim': r.x[1]} if both else {}))

    def transfer(fkw, pkw, libs_=None, c=cal, g=G, mode='T2'):
        q, _ = M.fit(list((libs_ or L).values()), **fkw)
        q0 = dict(q[0], logk=float(np.mean([x['logk'] for x in q])), logjlim=float(np.mean([x['logjlim'] for x in q])))
        if mode == 'T2':
            qq = calibrate(q0, c, pkw, both=pkw.get('cap', True))
        elif mode == 'T1':
            qq = calibrate(q0, c, pkw, both=False)
        else:
            qq = q0
        pred = np.log10(M.predict(g, qq, **pkw))
        return q[0], {str(m): [float(v) for v in row] for m, row in zip(g.mdi, pred)}

    def rmse(pred):
        e = []
        for m, row in pred.items():
            for v, t in zip(row, hidden['lib3_targets_A_cm2'][m]):
                if t < -1e-9:
                    e.append(v - np.log10(-t))
        return float(np.sqrt(np.mean(np.square(e))))

    out = {}
    qA, pA = transfer({}, {})
    qC, pC = transfer({'free': CW, 'use_wf': True}, {'use_wf': True})
    q046, p046 = transfer({'free': ('alpha_c', 'alpha_ads'), 'fixed': {'e_opt': -0.24 - 0.2207}}, {})
    qC046, _ = transfer({'free': ('alpha_c', 'alpha_ads', 'alpha_wf'), 'fixed': {'e_opt': -0.4607}, 'use_wf': True}, {'use_wf': True})
    qAf, _ = transfer({}, {}, libs_=Lf)
    qCf, _ = transfer({'free': CW, 'use_wf': True}, {'use_wf': True}, libs_=Lf)
    q024, p024 = transfer({'free': ('alpha_c', 'alpha_ads'), 'fixed': {'e_opt': -0.24}}, {})
    qmean, pmean = transfer({'mode': 'mean'}, {'mode': 'mean'})
    qnc, pnc = transfer({'cap': False}, {'cap': False})
    _, pT0 = transfer({}, {}, mode='T0')
    _, pT1 = transfer({}, {}, mode='T1')
    qlin, plin = transfer({}, {}, libs_=Ll, c=calL, g=GL)
    qB, pB = transfer({'free': ('alpha_c', 'alpha_wf'), 'fixed': {'alpha_ads': 0.0}, 'use_wf': True}, {'use_wf': True})
    qC024, _ = transfer({'free': ('alpha_c', 'alpha_ads', 'alpha_wf'), 'fixed': {'e_opt': -0.24}, 'use_wf': True}, {'use_wf': True})
    qsh, _ = M.fit(list(L.values()), free=CW, use_wf=True, shared_k=True)
    qcov, pcov = transfer({'free': ('alpha_c', 'alpha_ads'), 'fixed': {'e_opt': -0.24 - 0.7018}}, {})
    # alternative calibration design: farthest-point sampling started from the most Pt-rich spot
    alt_mdi = sorted(E3.mdi[fps(E3.comp, MAX_LIB3, start=3)].tolist())
    alt = M.Lib('Lib-3', set(alt_mdi))
    _, pAalt = transfer({}, {}, c=alt)
    # coverage shift estimators (equimolar decorations)
    eq = [p for p in paths if p['comp_id'].startswith('ag025')]
    om_all = np.mean([np.polyfit(p['theta_ml'], p['dG_rep_ev'], 1)[0] for p in eq])
    om_win = np.mean([np.polyfit(np.array(p['theta_ml'])[0:8], np.array(p['dG_rep_ev'])[0:8], 1)[0] for p in eq])
    om_org = np.mean([((np.array(p['theta_ml']) - 1 / 9) * p['dG_rep_ev']).sum() / ((np.array(p['theta_ml']) - 1 / 9) ** 2).sum() for p in eq])
    al_all = np.mean([np.polyfit(p['theta_ml'], p['dG_diff_ev'], 1)[0] for p in eq])
    last_first = np.mean([p['dG_diff_ev'][-1] - p['dG_diff_ev'][0] for p in eq])
    f = 8 / 9
    shift = float(om_all * f)
    ref = {'volcano_centre_eV': qAf['e_opt'], 'coverage_shift_eV': shift, 'alpha_wf': qCf['alpha_wf']}
    hidden['reference'] = ref
    hidden['calibration_spots_oracle'] = cal_mdi
    vals = [v for row in hidden['lib3_targets_A_cm2'].values() for v in row] + list(ref.values())
    abl = {'ablations': {
        'coverage_total_slope': {'coverage_shift_eV': float(al_all * f), 'volcano_centre_eV': -0.24 - float(al_all * f), 'lib3': pcov},
        'coverage_first_last': {'coverage_shift_eV': float(last_first)},
        'literature_centre': {'volcano_centre_eV': -0.24, 'lib3': p024},
        'mean_energy_descriptor': {'lib3': pmean, 'volcano_centre_eV': qmean['e_opt']},
        'no_transport_cap': {'lib3': pnc},
        'no_library_calibration': {'lib3': pT0},
        'prefactor_only_calibration': {'lib3': pT1},
        'linear_combination_descriptor': {'lib3': plin, 'volcano_centre_eV': qlin['e_opt']},
        'wf_only_coefficient': {'alpha_wf': qB['alpha_wf'], 'lib3': pB},
        'wf_with_dilute_binding': {'alpha_wf': qC024['alpha_wf']},
        'pooled_without_library_offsets': {'alpha_wf': qsh[0]['alpha_wf']},
    }, 'robustness': {
        'combined_model_predictions': {'lib3': pC, 'volcano_centre_eV': qC['e_opt'], 'alpha_wf': qC['alpha_wf']},
        'centre_fixed_from_dft_coverage': {'lib3': p046, 'volcano_centre_eV': -0.4607, 'alpha_wf': qC046['alpha_wf']},
        'alternative_calibration_spots': {'lib3': pAalt},
        'coverage_shift_window_fit': {'coverage_shift_eV': float(om_win * f)},
        'coverage_shift_through_dilute_point': {'coverage_shift_eV': float(om_org * f)},
        'oracle_subset_reference': {'volcano_centre_eV': qA['e_opt'], 'alpha_wf': qC['alpha_wf']},
    }}
    calib = {'lib3_rmse_log10': {'oracle_A': rmse(pA), 'combined_C': rmse(pC), 'centre_fixed_-0.461': rmse(p046),
                                 'alt_calibration': rmse(pAalt), 'mean_energy': rmse(pmean), 'no_cap': rmse(pnc),
                                 'T0': rmse(pT0), 'T1_prefactor_only': rmse(pT1), 'literature_centre': rmse(p024),
                                 'linear_combination': rmse(plin), 'wf_only': rmse(pB), 'total_slope_centre': rmse(pcov),
                                 'tolerance': THRESHOLDS['lib3_rmse_log10']},
             'volcano_centre_eV': {'reference_all_rows': qAf['e_opt'], 'oracle_subset': qA['e_opt'], 'combined': qC['e_opt'],
                                   'mean_energy': qmean['e_opt'], 'linear_combination': qlin['e_opt'],
                                   'literature': -0.24, 'total_slope': -0.24 - al_all * f, 'tolerance': THRESHOLDS['volcano_centre_window_eV']},
             'coverage_shift_eV': {'reference_all9': shift, 'window_0.1_0.9': om_win * f, 'through_dilute_point': om_org * f,
                                   'total_slope_alpha': al_all * f, 'first_last_differential': last_first,
                                   'tolerance': THRESHOLDS['coverage_shift_window_eV']},
             'alpha_wf': {'reference_all_rows': qCf['alpha_wf'], 'oracle_subset': qC['alpha_wf'], 'centre_fixed': qC046['alpha_wf'],
                          'wf_only': qB['alpha_wf'], 'dilute_binding': qC024['alpha_wf'], 'pooled_no_offsets': qsh[0]['alpha_wf'],
                          'tolerance': THRESHOLDS['alpha_wf_window']},
             'oracle_calibration_spots': cal_mdi, 'alternative_calibration_spots': alt_mdi}
    # ---------------- write ----------------
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/agent/descriptors.json').write_text(json.dumps(descriptors, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(vals) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps(abl, indent=1) + '\n')
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    print(json.dumps(calib, indent=1))


if __name__ == '__main__':
    main()
