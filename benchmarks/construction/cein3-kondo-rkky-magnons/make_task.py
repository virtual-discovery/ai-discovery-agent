"""Generate the staged, lab, private and ablation data (and design.json) for cein3-kondo-rkky-magnons.

    python3 -B fetch.py && (run calibrate.py in the runtime image) && python3 -B make_task.py
Reads the raw archive files (CEIN3_RAW, default ~/.cache/ai-theorist/cein3-kondo-rkky-magnons/raw) and
calibration_fits.json; writes catalog.json, data.json, private/hidden.json, private/hidden_values.json,
solution/ablation_values.json, both public.json files and design.json, and prints the calibration table.
"""
import json
import os
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/cein3-kondo-rkky-magnons'
RAW = Path(os.environ.get('CEIN3_RAW', Path.home()/'.cache/ai-theorist/cein3-kondo-rkky-magnons/raw'))
L = lambda n: np.loadtxt(RAW/n)
CAL = json.loads((HERE/'calibration_fits.json').read_text())
BUDGET = 60
COSTS = {'he_spectrum': 1, 'hr_spectrum': 1, 'hr_cut': 2, 'band_structure': 1}
THRESHOLDS = {'rm_map_chi2r': 1.7, 'rm_cut_1p4meV_chi2r': 1.9, 'integrated_chi2r': 5.5,
              'bandwidth_meV': [2.35, 3.25], 'velocity_bandwidth_ratio': [0.75, 3.0], 'J1_nearest_neighbour_meV': [0.8, 1.2]}
S_EFF = 0.5
R = np.array([0.5, 0.5, 0.5])
R3 = np.sqrt(0.75)
HE_FILES = {'RG': 'HighEnergy_Int_RG.dat', 'GX': 'HighEnergy_Int_GX.dat', 'XM': 'HighEnergy_IntXM.dat', 'MG': 'HighEnergy_Int_MG.dat'}
HE_PERP = {'RG': 0.09, 'MG': 0.09, 'GX': 0.17, 'XM': 0.17}
HR_DIR = {'RG': [-1, -1, -1], 'RX': [0, -1, -1], 'RM': [0, 0, -1]}


def unit(v):
    v = np.array(v, float); return v / np.linalg.norm(v)


def he_center(path, s):
    if path == 'RG': return R * (1 - s / R3)
    if path == 'GX': return np.array([s - R3, 0, 0])
    if path == 'XM': return np.array([0.5, s - R3 - 0.5, 0])
    f = (s - R3 - 1.0) / np.sqrt(0.5); return np.array([0.5, 0.5, 0]) * (1 - f)


def clean(row, nd=6):
    return [None if not np.isfinite(v) else float('%.*g' % (nd, v)) for v in row]


def r5(x):
    if isinstance(x, list): return [r5(v) for v in x]
    if isinstance(x, dict): return {k: r5(v) for k, v in x.items()}
    return float('%.5g' % x) if isinstance(x, float) else x


def build_data():
    data = {'costs': COSTS, 'units': "chi''(Q,E) in the authors' absolute units (mu_B^2/meV); the 12-meV and the high-resolution data are normalised separately",
            'he': {}, 'hr': {}, 'cuts': {}}
    for path, f in HE_FILES.items():
        I, Q, E = L(f), L('HighEnergy_Q_%s.dat' % path), L('HighEnergy_E_%s.dat' % path)
        rows = {}
        for i, s in enumerate(Q):
            if np.isfinite(I[i]).any():
                rows[str(i)] = {'index': i, 's': round(float(s), 5), 'q_center_rlu': np.round(he_center(path, s), 5).tolist(), 'chi2': clean(I[i])}
        data['he'][path] = {'energies_meV': np.round(E, 4).tolist(), 'energy_bin_meV': round(float(E[2] - E[1]), 4),
                            'perpendicular_halfwidth_rlu': HE_PERP[path], 'rows': rows}
    for path in ['RG', 'RX']:
        I, Q, E = L('HighResolution_Int_%s.dat' % path), L('HighResolution_Q_%s.dat' % path), L('HighResolution_E_%s.dat' % path)
        data['hr'][path] = {'energies_meV': np.round(E, 4).tolist(), 'energy_bin_meV': round(float(E[1] - E[0]), 4),
                            'rows': {str(i): {'index': i, 'q_rlu': round(float(q), 4), 'q_center_rlu': np.round(R + q * unit(HR_DIR[path]), 5).tolist(),
                                                  'chi2': clean(I[i])} for i, q in enumerate(Q)}}
        data['cuts'][path] = {}
        for tag, e in [('0p6meV', 0.6), ('1p4meV', 1.4)]:
            q, y, err = L('HighResolution_%s_ConstEcut_%s.dat' % (path, tag)).T
            data['cuts'][path]['%.1f' % e] = {'energy_meV': e, 'energy_halfwindow_meV': 0.1, 'q_rlu': np.round(q, 5).tolist(),
                                              'chi2': clean(y), 'error': clean(err)}
    fb, ob = L('fbands.dat'), L('allotherbands.dat')
    data['bands'] = {'path_coordinate_f': np.round(fb[:, 0], 4).tolist(), 'f_bands_eV': np.round(fb[:, 1:], 5).tolist(),
                     'path_coordinate_other': np.round(ob[:, 0], 4).tolist(), 'other_bands_eV': np.round(ob[:, 1:], 4).tolist(),
                     'note': 'Archived band structure of the authors\' multi-orbital periodic Anderson model on the path R-Gamma-X-M-Gamma '
                             '(path coordinate in r.l.u., R = 0). f_bands_eV: the 14 Ce 4f bands; other_bands_eV: all remaining bands. '
                             'Energies in eV as archived (relative to the chemical potential).'}
    return data


def build_catalog(data):
    starts = {'RG': ['sqrt(3)/2', -0.025], 'GX': ['sqrt(3)/2', 0.025], 'XM': ['sqrt(3)/2 + 1/2', 0.025], 'MG': ['sqrt(3)/2 + 1', 0.025]}
    he_avail = {p: {'s_of_index': 's = %s + index * (%g)' % tuple(starts[p]), 'indices': [r['index'] for r in data['he'][p]['rows'].values()]}
                for p in data['he']}
    hr_avail = {p: {'q_of_index': 'q = -0.33 + 0.02 * index', 'indices': [r['index'] for r in data['hr'][p]['rows'].values()]} for p in data['hr']}
    qrm, erm = L('HighResolution_Q_RM.dat'), L('HighResolution_E_RM.dat')
    qcut = L('HighResolution_RM_ConstEcut_1p4meV.dat')[:, 0]
    sint = L('Experiment_IntegralOverchi.dat')[:, 0]
    seg = lambda s: 'RG' if s <= R3 + 1e-9 else 'GX' if s <= R3 + 0.5 + 1e-9 else 'XM' if s <= R3 + 1 + 1e-9 else 'MG'
    return {
        'material': 'CeIn3, simple cubic (Pm-3m), a = 4.689 A; antiferromagnetic below T_N = 10 K with propagation vector R = (1/2,1/2,1/2); '
                    'Ce3+ crystal-field ground state Gamma7 doublet (Gamma8 quartet about 12 meV higher)',
        'experiment': {'instrument': 'CNCS time-of-flight spectrometer (SNS), co-aligned single crystals, [1,-1,0] vertical, T = 1.8 K',
                       'chi2_definition': "Archived chi''(Q,E): intensities corrected for absorption, put on an absolute scale with incoherent scattering, "
                                          'background-subtracted, divided by the Bose factor and by the Ce3+ magnetic form factor, and corrected for kf/ki. '
                                          'Background subtraction leaves negative values in places.',
                       'settings': {
                           'high_energy': {'incident_energy_meV': 12.0, 'energy_resolution': 'not archived',
                                           'binning': 'paths R-Gamma, Gamma-X, X-M, M-Gamma; bins of +-0.0125 r.l.u. along the path; intensity integrated '
                                                      'within +-0.17 r.l.u. along the two directions perpendicular to Gamma-X and X-M and within +-0.09 r.l.u. '
                                                      'perpendicular to R-Gamma and M-Gamma (authors\' Methods)'},
                           'high_resolution': {'incident_energy_meV': 3.315, 'elastic_energy_resolution_fwhm_ueV': 106,
                                               'binning': 'lines through R along R-Gamma, R-X and R-M; bins of +-0.01 r.l.u. along the line; authors\' Methods: '
                                                          '"intensity was integrated within a distance of 0.08 r.l.u. along the two reciprocal space directions '
                                                          'perpendicular to the reduced momentum transfer q. For the cuts in the lower panels, intensity was '
                                                          'integrated over energies +-0.1 meV"'}}},
        'coordinates': {'rlu': 'momenta in reciprocal-lattice units, 1 r.l.u. = 2 pi / a; reduced wave vectors q = (h,k,l)',
                        'path_coordinate_s': 'distance along R-Gamma-X-M-Gamma in r.l.u.: R = (1/2,1/2,1/2) at s = 0, Gamma at sqrt(3)/2, '
                                             'X = (1/2,0,0) at sqrt(3)/2 + 1/2, M = (1/2,1/2,0) at sqrt(3)/2 + 1, Gamma at sqrt(3)/2 + 1 + sqrt(1/2). '
                                             'On the R-Gamma file, s < 0 lies beyond R on the same (1,1,1) line. Each he_spectrum result gives s and the bin centre q.',
                        'high_resolution_q': 'signed distance from R in r.l.u. along the line: RG towards (0,0,0) for q > 0 (direction (-1,-1,-1)/sqrt3), '
                                             'RX towards (1/2,0,0) (direction (0,-1,-1)/sqrt2), RM towards (1/2,1/2,0) (direction (0,0,-1))'},
        'requests': {
            'he_spectrum': {'request': {'kind': 'he_spectrum', 'path': '<RG|GX|XM|MG>', 'index': '<integer>'}, 'cost': COSTS['he_spectrum'],
                            'returns': "one constant-Q energy spectrum of the 12-meV data: its path coordinate s, bin centre q, energies and chi'' "
                                       '(null where not measured)',
                            'available': he_avail},
            'hr_spectrum': {'request': {'kind': 'hr_spectrum', 'path': '<RG|RX>', 'index': '<integer>'}, 'cost': COSTS['hr_spectrum'],
                            'returns': "one constant-q energy spectrum of the high-resolution data along R-Gamma or R-X",
                            'available': hr_avail},
            'hr_cut': {'request': {'kind': 'hr_cut', 'path': '<RG|RX>', 'energy_meV': '<0.6|1.4>'}, 'cost': COSTS['hr_cut'],
                       'returns': "a high-resolution constant-energy cut (E +- 0.1 meV) along the line, with error bars"},
            'band_structure': {'request': {'kind': 'band_structure'}, 'cost': COSTS['band_structure'],
                               'returns': 'the archived band structure of the authors\' periodic Anderson model along R-Gamma-X-M-Gamma'}},
        'budget': BUDGET,
        'not_served': 'The high-resolution R-M map, the R-M constant-energy cuts and the Q-integrated chi\'\' are withheld for grading. '
                      'The authors\' theory files are not served.',
        'to_predict': {
            'rm_map': {'description': "high-resolution chi''(q,E) along R-M: rows = q_rlu (bins of +-0.01 r.l.u.), columns = energies_meV "
                                      '(bins of +-0.025 meV); binning and resolution as for the other high-resolution data. Pixels at E <= 0.25 meV '
                                      'and pixels outside the measured range are not graded.',
                       'q_rlu_of_row': 'q = -0.33 + 0.02 * row, row = 0..%d' % (len(qrm) - 1),
                       'energy_meV_of_column': 'E = 0.03 + 0.05 * column, column = 0..%d' % (len(erm) - 1)},
            'rm_cut_1p4meV': {'description': "high-resolution constant-energy cut along R-M at E = 1.4 +- 0.1 meV (bins of +-0.005 r.l.u. along the line)",
                              'q_rlu_of_point': 'q = %.6f + 0.0099 * point, point = 0..%d' % (qcut[0], len(qcut) - 1)},
            'integrated_chi': {'description': "chi''(Q,E) integrated over the magnon energies, binned like the 12-meV data at these path coordinates",
                               'points': [{'segment': seg(s), 'q_center_rlu': np.round(he_center(seg(s), s), 1).tolist(),
                                           'perpendicular_halfwidth_rlu': HE_PERP[seg(s)]} for s in sint]}},
    }


def build_private():
    I, Q, E = L('HighResolution_Int_RM.dat'), L('HighResolution_Q_RM.dat'), L('HighResolution_E_RM.dat')
    m = E > 0.25
    pairs = [(i, int(np.argmin(abs(Q + Q[i])))) for i in range(len(Q)) if Q[i] > 0]
    d = np.array([I[i] - I[j] for i, j in pairs if abs(Q[j] + Q[i]) < 1e-6])[:, m]
    sigma = float(np.nanstd(d) / np.sqrt(2))
    cq, cy, ce = L('HighResolution_RM_ConstEcut_1p4meV.dat').T
    s, iy, ie = L('Experiment_IntegralOverchi.dat').T
    ref = CAL['reference_all_data']['props']
    j1 = CAL['J1_HEonly']['params']['SJ'][0] / S_EFF
    hidden = {'rm_map': {'q_rlu': Q.tolist(), 'energies_meV': E.tolist(), 'emin_meV': 0.25, 'sigma': sigma, 'values': [clean(r, 8) for r in I]},
              'rm_cut_1p4meV': {'q_rlu': cq.tolist(), 'y': cy.tolist(), 'e': ce.tolist()},
              'integrated': {'s': s.tolist(), 'y': iy.tolist(), 'e': ie.tolist()},
              'thresholds': THRESHOLDS, 'exchange_range': 'long_range',
              'reference': {'bandwidth_meV': ref['W'], 'velocity_bandwidth_ratio': ref['eta'], 'velocity_meV_rlu': ref['v'],
                            'J1_nearest_neighbour_meV': j1, 'rm_map_sigma': sigma}}
    values = [v for r in hidden['rm_map']['values'] for v, e in zip(r, E) if v is not None and e > 0.25]
    values += cy.tolist() + ce.tolist() + iy.tolist() + ie.tolist() + [ref['W'], ref['eta'], ref['v'], j1]
    return hidden, values


def box_eta():
    """Simplest cut analysis: each visible cut fitted by a box (half-width r, erf edges of 0.01 r.l.u.); v = E/r."""
    from scipy.optimize import least_squares
    from scipy.special import erf
    vs = []
    for path in ['RG', 'RX']:
        for tag, e0 in [('0p6meV', 0.6), ('1p4meV', 1.4)]:
            q, y, e = L('HighResolution_%s_ConstEcut_%s.dat' % (path, tag)).T
            m = np.abs(q) < 0.3
            f = lambda x: (x[0] * 0.5 * (erf((q[m] - x[1] + x[2]) / (np.sqrt(2) * 0.01)) - erf((q[m] - x[1] - x[2]) / (np.sqrt(2) * 0.01))) + x[3] - y[m]) / e[m]
            best = min((least_squares(f, [0.8, 0.0, r0, 0.0], bounds=([0, -0.05, 0.001, -1], [5, 0.05, 0.3, 1])) for r0 in (0.02, 0.04, 0.08)),
                       key=lambda r: r.cost)
            vs.append(e0 / best.x[2])
    return float(np.mean(vs)), [round(v, 2) for v in vs]


def build_ablations():
    C = CAL
    pred = lambda k: {'rm_map_prediction': r5(C[k]['pred']['rm_map']), 'rm_cut_1p4meV_prediction': r5(C[k]['pred']['rm_cut_1p4meV']),
                      'integrated_chi_prediction': r5(C[k]['pred']['integrated'])}
    sc = lambda k: {'bandwidth_meV': r5(C[k]['props']['W']), 'velocity_bandwidth_ratio': r5(C[k]['props']['eta'])}
    qrm, erm = L('HighResolution_Q_RM.dat'), L('HighResolution_E_RM.dat')
    flat = {'rm_map_prediction': [[1.0] * len(erm) for _ in qrm], 'rm_cut_1p4meV_prediction': [1.0] * 70, 'integrated_chi_prediction': [1.0] * 16}
    abl = {'J1_HEonly': {**pred('J1_HEonly'), **sc('J1_HEonly'), 'exchange_range': 'nearest_neighbour'},
           'J123_short': {**pred('J123_short'), **sc('J123_short'), 'exchange_range': 'short_range'},
           'J6shell_short': {**pred('J6shell_short'), **sc('J6shell_short'), 'exchange_range': 'short_range'},
           'noperp_J1LR': {**pred('noperp_J1LR'), **sc('noperp_J1LR')},
           'J1_forced_all': {**pred('J1_forced_all'), **sc('J1_forced_all'), 'J1_nearest_neighbour_meV': r5(C['J1_forced_all']['params']['SJ'][0] / S_EFF),
                             'exchange_range': 'nearest_neighbour'},
           'J1LR_HRonly': {**pred('J1LR_HRonly'), **sc('J1LR_HRonly')},
           'free_ion_spin': {'J1_nearest_neighbour_meV': r5(C['J1_HEonly']['params']['SJ'][0] / 2.5)},
           'incoherent_continuum': {**flat, 'exchange_range': 'no_spin_waves'},
           'J10shell_free': {**pred('J10shell'), **sc('J10shell')}}
    rob = {'J123LR': {**pred('J123LR'), **sc('J123LR')},
           'J1LR_h004': {**pred('J1LR_h004'), **sc('J1LR_h004')},
           'J1LR': {**pred('J1LR'), **sc('J1LR')},
           'J1_from_bandwidth_2p75': {'J1_nearest_neighbour_meV': 2.75 / (6 * S_EFF), 'bandwidth_meV': 2.75},
           'eta_from_cut_box_widths': {'velocity_bandwidth_ratio': r5(BOX[0] / (2 * np.pi * C['reference_all_data']['props']['W']))}}
    return abl, rob


BOX = None


def write_design(hidden, abl, rob):
    C, th = CAL, THRESHOLDS
    sc = lambda k, f: C[k]['scores'][f]
    pr = lambda k, f: round(C[k]['props'][f], 3)
    oracle = json.loads((HERE/'oracle_result.json').read_text()) if (HERE/'oracle_result.json').exists() else {}
    design = {
        'task_id': 'cein3-kondo-rkky-magnons', 'catalog_id': 'cein3-kondo-rkky-magnons', 'start': 'released',
        'source_group': 'simeth_2023_cein3',
        'title': 'CeIn3: nearest-neighbour versus long-range RKKY exchange in the magnon spectrum',
        'decision_points': [
            {'id': 'dp1', 'fork': 'Which measurements discriminate: the 12-meV zone-boundary maps or the high-resolution data near R',
             'wrong_choice': 'Spend the budget on 12-meV spectra and fit a nearest-neighbour model, which describes them well',
             'consequence': 'eta fixed at 0.577; the predicted R-M map (chi2r %.2f), R-M 1.4-meV cut (%.2f) and integrated chi\'\' (%.1f) all fail'
                            % (sc('J1_HEonly', 'rm_map'), sc('J1_HEonly', 'rm_cut_1p4meV'), sc('J1_HEonly', 'integrated')),
             'ablation': 'nearest_neighbour_model_from_12meV_data_only', 'difficulty_type': 'discriminating experiment under a binding budget'},
            {'id': 'dp2', 'fork': 'Exchange range: short-range shells (up to third or sixth neighbours) or a long-range staggered term',
             'wrong_choice': 'Add second and third (or up to sixth) neighbour couplings until the visible spectra fit',
             'consequence': 'The R-M predictions pass, but the Q-integrated chi\'\' is predicted strongly modulated (chi2r %.1f and %.1f, limit %.1f) and eta stays below 1 (%.2f, %.2f)'
                            % (sc('J123_short', 'integrated'), sc('J6shell_short', 'integrated'), th['integrated_chi2r'], pr('J123_short', 'eta'), pr('J6shell_short', 'eta')),
             'ablation': 'short_range_J1_J3_model', 'difficulty_type': 'model class; transfer to a withheld observable'},
            {'id': 'dp3', 'fork': 'Momentum binning: the archived spectra are averaged over wide perpendicular windows',
             'wrong_choice': 'Model each spectrum as taken exactly on the line (no perpendicular integration)',
             'consequence': 'The broad constant-q profiles near R force a cusp-like J(q): eta %.2f (limit %.1f) and integrated chi2r %.2f'
                            % (pr('noperp_J1LR', 'eta'), th['velocity_bandwidth_ratio'][1], sc('noperp_J1LR', 'integrated')),
             'ablation': 'momentum_integration_window_ignored', 'difficulty_type': 'evidence in the reduction metadata that must be discovered'},
            {'id': 'dp4', 'fork': 'Bandwidth: zone-boundary data versus a compromise fit',
             'wrong_choice': 'Force one nearest-neighbour model through the 12-meV and high-resolution data together',
             'consequence': 'The steep low-energy data inflate J1: W %.2f meV (window %s), J1 %.2f meV, eta 0.577'
                            % (pr('J1_forced_all', 'W'), th['bandwidth_meV'], C['J1_forced_all']['params']['SJ'][0] / S_EFF),
             'ablation': 'nearest_neighbour_model_forced_through_all_data', 'difficulty_type': 'cascade'},
            {'id': 'dp5', 'fork': 'Effective spin of the Ce3+ ground state',
             'wrong_choice': 'Use the free-ion J = 5/2 instead of the Gamma7 doublet (S = 1/2)',
             'consequence': 'J1 = %.3f meV instead of about 1 meV' % abl['free_ion_spin']['J1_nearest_neighbour_meV'],
             'ablation': 'free_ion_J_5_2_as_effective_spin', 'difficulty_type': 'convention'},
            {'id': 'dp6', 'fork': 'Spin waves or an incoherent Kondo continuum',
             'wrong_choice': 'Read the broad, plateau-like spectra near R as a Q-independent continuum',
             'consequence': 'Non-dispersive predictions fail the R-M map (chi2r %.2f) and cut (%.2f)'
                            % (C['continuum']['scores']['rm_map'], C['continuum']['scores']['rm_cut_1p4meV']),
             'ablation': 'broad_spectra_read_as_incoherent_continuum', 'difficulty_type': 'competing explanation'},
            {'id': 'dp7', 'fork': 'How many free exchange parameters',
             'wrong_choice': 'Fit ten free real-space shells and trust the result',
             'consequence': 'Underdetermined: the fit places a spurious branch off the measured paths, W %.2f meV and eta %.2f (integrated chi2r %.2f); '
                            'another start of the same fit lands near W 2.9 meV, so the shell signs and the bandwidth are not stable'
                            % (pr('J10shell', 'W'), pr('J10shell', 'eta'), sc('J10shell', 'integrated')),
             'ablation': 'ten_free_shells_overfit', 'difficulty_type': 'underdetermined fit'},
        ],
        'graded_fields': [
            'rm_map_prediction: reduced chi2 (best scale and offset) <= %.2f against the withheld high-resolution R-M map' % th['rm_map_chi2r'],
            'rm_cut_1p4meV_prediction: reduced chi2 (best scale and offset) <= %.2f against the withheld R-M 1.4-meV cut' % th['rm_cut_1p4meV_chi2r'],
            'integrated_chi_prediction: reduced chi2 (best scale) <= %.1f against the withheld Q-integrated chi\'\'' % th['integrated_chi2r'],
            'bandwidth_meV in %s' % th['bandwidth_meV'], 'velocity_bandwidth_ratio in %s' % th['velocity_bandwidth_ratio'],
            'J1_nearest_neighbour_meV in %s' % th['J1_nearest_neighbour_meV'], 'exchange_range == long_range'],
        'length_estimate': {'oracle_lab_evaluations': oracle.get('evaluations'), 'oracle_charge': oracle.get('charge'), 'budget': BUDGET,
                            'catalog_worth_units': CATALOG_WORTH, 'tool_calls': '50-80', 'gpt6astra_cost_usd': '3-5',
                            'oracle_container_seconds': oracle.get('seconds'), 'pilot': 'pending (coordinator)'},
        'runtime_dependencies': ['numpy and scipy (runtime image); the lab is a lookup of archived data, grading is numpy'],
        'data_provenance': {'archive_url': 'https://zenodo.org/records/10146787', 'doi': '10.5281/zenodo.10146787',
                            'paper_doi': '10.1038/s41467-023-43947-z', 'files_used': ['all 48 .dat files (members.json)'],
                            'md5': {f['name']: f['md5'] for f in json.loads((HERE/'members.json').read_text())['files']}, 'license': 'CC BY 4.0'},
        'exceptions': {},
        'notes': ('Adapted from the design: the lab serves archived data (constant-Q spectra, cuts, band structure) instead of an LSWT '
                  'calculator, because the agent can run LSWT itself and the budget would not bind. Graded targets are the withheld '
                  'archived R-M map, R-M 1.4-meV cut and Q-integrated chi\'\' (reduced chi2 after best scale), plus windows on W, eta and J1. '
                  'The gap at R is not graded: the archived bins below ~0.25 meV are dominated by the subtracted elastic line and fits '
                  'give 0-0.43 meV depending on the energy cut. The J1..J4 sign pattern is not graded: it is not uniquely recoverable. '
                  'The RX constant-energy cuts are served (the design withheld them) because they are the only visible cuts besides RG. '
                  'eta windows are wide because the velocity depends on the reading of the 0.08 r.l.u. window (22-33 meV/r.l.u.), and the '
                  'paper\'s v = 38.6 used cuts that are not archived.'),
        'calibration': {
            'rm_map_chi2r': {'oracle': oracle.get('rm_map_chi2r'), 'reference_all_data': sc('reference_all_data', 'rm_map'),
                             'robustness': {k: sc(k, 'rm_map') for k in ['J1LR', 'J123LR', 'J1LR_h004']},
                             'nearest_ablation': {'nearest_neighbour_12meV_only': sc('J1_HEonly', 'rm_map'), 'continuum': C['continuum']['scores']['rm_map'],
                                                  'no_integration_window (passes)': sc('noperp_J1LR', 'rm_map')},
                             'tolerance': '<= %.2f' % th['rm_map_chi2r']},
            'rm_cut_1p4meV_chi2r': {'oracle': oracle.get('rm_cut_1p4meV_chi2r'), 'reference_all_data': sc('reference_all_data', 'rm_cut_1p4meV'),
                                    'robustness': {k: sc(k, 'rm_cut_1p4meV') for k in ['J1LR', 'J123LR', 'J1LR_h004']},
                                    'nearest_ablation': {'nearest_neighbour_12meV_only': sc('J1_HEonly', 'rm_cut_1p4meV'), 'continuum': C['continuum']['scores']['rm_cut_1p4meV'],
                                                         'no_integration_window (passes)': sc('noperp_J1LR', 'rm_cut_1p4meV')},
                                    'tolerance': '<= %.2f' % th['rm_cut_1p4meV_chi2r']},
            'integrated_chi2r': {'oracle': oracle.get('integrated_chi2r'), 'reference_all_data': sc('reference_all_data', 'integrated'),
                                 'robustness': {k: sc(k, 'integrated') for k in ['J1LR', 'J123LR', 'J1LR_h004']},
                                 'nearest_ablation': {'no_integration_window': sc('noperp_J1LR', 'integrated'), 'J1_J3': sc('J123_short', 'integrated'),
                                                      'nearest_neighbour': sc('J1_HEonly', 'integrated'), 'six_shells': sc('J6shell_short', 'integrated'),
                                                      'flat_continuum (passes)': C['continuum']['scores']['integrated']},
                                 'tolerance': '<= %.1f' % th['integrated_chi2r']},
            'bandwidth_meV': {'reference_all_data': pr('reference_all_data', 'W'), 'oracle': oracle.get('W'),
                              'robustness': {**{k: pr(k, 'W') for k in ['J1LR', 'J123LR', 'J1LR_h004', 'J1_HEonly']}, 'paper': 2.75},
                              'nearest_ablation': {'J1_forced_all': pr('J1_forced_all', 'W'), 'six_shells': pr('J6shell_short', 'W'), 'ten_free_shells': pr('J10shell', 'W')},
                              'tolerance': str(th['bandwidth_meV'])},
            'velocity_bandwidth_ratio': {'reference_all_data': pr('reference_all_data', 'eta'), 'oracle': oracle.get('eta'),
                                         'robustness': {**{k: pr(k, 'eta') for k in ['J1LR', 'J123LR', 'J1LR_h004']},
                                                        'cut_box_widths': rob['eta_from_cut_box_widths']['velocity_bandwidth_ratio'], 'paper': 2.23},
                                         'nearest_ablation': {'nearest_neighbour': 0.577, 'no_integration_window': pr('noperp_J1LR', 'eta'), 'ten_free_shells': pr('J10shell', 'eta'),
                                                              'six_shells (passes)': pr('J6shell_short', 'eta'), 'J1_J3 (passes)': pr('J123_short', 'eta')},
                                         'tolerance': str(th['velocity_bandwidth_ratio'])},
            'J1_nearest_neighbour_meV': {'reference': round(hidden['reference']['J1_nearest_neighbour_meV'], 4), 'oracle': oracle.get('J1'),
                                         'robustness': {'from_W_2.75': round(2.75 / 3, 4)},
                                         'nearest_ablation': {'forced_all_data': round(C['J1_forced_all']['params']['SJ'][0] / S_EFF, 4),
                                                              'free_ion_J_5_2': abl['free_ion_spin']['J1_nearest_neighbour_meV']},
                                         'tolerance': str(th['J1_nearest_neighbour_meV'])},
        },
    }
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')


CATALOG_WORTH = None


def main():
    global BOX
    BOX = box_eta()
    print('box-width velocities', BOX)
    data = build_data()
    catalog = build_catalog(data)
    hidden, values = build_private()
    abl, rob = build_ablations()
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}) + '\n')
    for f, extra in [('environment/lab/public.json', {'reserved_validation_charge': 0}), ('tests/public.json', {'min_report_chars': 400})]:
        (TASK/f).write_text(json.dumps({'task_id': 'cein3-kondo-rkky-magnons', 'solver_budget': BUDGET, **extra}, indent=1) + '\n')
    total = sum(len(v['rows']) for v in data['he'].values()) * COSTS['he_spectrum'] + \
        sum(len(v['rows']) for v in data['hr'].values()) * COSTS['hr_spectrum'] + 4 * COSTS['hr_cut'] + COSTS['band_structure']
    print('catalog worth', total, 'units; budget', BUDGET)
    global CATALOG_WORTH
    CATALOG_WORTH = total
    write_design(hidden, abl, rob)
    print('%-22s %8s %8s %8s %7s %7s' % ('variant', 'rm_map', 'rm_cut', 'integ', 'W', 'eta'))
    for k, v in CAL.items():
        if isinstance(v, dict) and 'scores' in v:
            s = v['scores']; pr = v.get('props', {})
            print('%-22s %8.3f %8.3f %8.3f %7s %7s' % (k, s['rm_map'], s['rm_cut_1p4meV'], s['integrated'],
                                                    '%.3f' % pr['W'] if pr else '-', '%.3f' % pr['eta'] if pr else '-'))
    print('reference', hidden['reference'])


if __name__ == '__main__':
    main()
