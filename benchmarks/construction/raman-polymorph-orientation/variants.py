"""Wrong paths (ablations) and defensible alternatives (robustness), computed from the oracle's purchases.

Called by make_task.py with the replayed purchases (identical to what the lab serves the reference solution).
Writes solution/ablation_values.json and calibration.json (this folder) and prints the calibration table.
"""
import json
import math
import sys
import numpy as np

sys.dont_write_bytecode = True
import analysis as A        # noqa: E402
import model as M           # noqa: E402

MAIN, ALPHA, GAMMA = 'mp-7456', 'mp-2843', 'mp-7457'


def payload_from(out, minority=None):
    return {'main_phase': out['main_phase'], 'minority_phase': minority if minority is not None else out.get('minority_phase', 'none'),
            'frequency_scale': out['frequency_scale'], 'c_axis_lab': out['c_axis_lab'],
            'faceB_strongest_mode_ratio_tt_over_zz': out['ratio_tt_over_zz'],
            'faceB_highest_optical_phonon_cm1': out['highest_optical_cm1']}


def errors(p, ref):
    return {'s_rel': p['frequency_scale'] / ref['frequency_scale'] - 1,
            'axis_deg': M.axis_error_deg(p['c_axis_lab'], ref['c_axis_lab']),
            'ratio_log': math.log(p['faceB_strongest_mode_ratio_tt_over_zz'] / ref['faceB_ratio_tt_over_zz']),
            'lo_rel': p['faceB_highest_optical_phonon_cm1'] / ref['faceB_highest_optical_cm1'] - 1,
            'main': p['main_phase'], 'minority': p['minority_phase']}


def nearest_line_matching(catalog, results):
    """Wrong path: compare candidates by the distance of observed peaks to the nearest predicted Raman line."""
    spectra = A.prepare(results)
    peaks = set()
    for sp in spectra:
        y, nu = sp['y'] / sp['accum'], sp['nu']
        for i in range(2, len(y) - 2):
            if nu[i] > 30 and y[i] == y[i - 2:i + 3].max() and y[i] - np.percentile(y, 20) > 150:
                peaks.add(round(float(nu[i]), 1))
    obs = np.array(sorted(peaks))
    best = {}
    for mpid, rec in catalog['polymorph_records'].items():
        e = A.parse(rec)
        f, R = A.modes(e)
        f = f[(R ** 2).sum((1, 2)) > 0.05]
        grid = np.linspace(0.95, 1.05, 1001)
        score = [np.mean(np.min(np.abs(obs[:, None] - g * f[None, :]), axis=1)) for g in grid]
        k = int(np.argmin(score))
        best[mpid] = {'scale': float(grid[k]), 'mean_abs_dev': float(score[k]),
                      'raw_mean_abs_dev': float(np.mean(np.min(np.abs(obs[:, None] - f[None, :]), axis=1)))}
    return best


def run(catalog, results, ref, thresholds, lab_axes_only, task, here):
    log = print
    oracle = A.analyse(catalog, results, log=log)
    P0 = payload_from(oracle)
    parsed = {m: A.parse(r) for m, r in catalog['polymorph_records'].items()}
    e = parsed[MAIN]
    Q, s = np.array(oracle['Q']), oracle['frequency_scale']
    spectra = A.prepare(results)
    abl, rob, info = {}, {}, {}

    # dp1: raw database frequencies (no global scale factor)
    abl['raw_dft_frequencies'] = {'frequency_scale': 1.0,
                                  'faceB_highest_optical_phonon_cm1': oracle['highest_optical_cm1'] / s}
    # dp2: pure scale through the origin; the axis offset is absorbed into s
    pure = A.analyse(catalog, results, offset='zero', candidates=[MAIN], minority=False, log=log)
    abl['scale_without_axis_offset'] = {'frequency_scale': pure['frequency_scale'],
                                        'faceB_highest_optical_phonon_cm1': pure['highest_optical_cm1']}
    # dp3: nearest-line peak matching (positions only)
    match = nearest_line_matching(catalog, results)
    info['nearest_line_matching'] = match
    winner = min(match, key=lambda m: match[m]['mean_abs_dev'])
    raw_winner = min(match, key=lambda m: match[m]['raw_mean_abs_dev'])
    abl['nearest_line_peak_matching'] = {'main_phase': winner, 'frequency_scale': match[winner]['scale']}
    info['nearest_line_winner'], info['nearest_line_raw_winner'] = winner, raw_winner
    # dp4: the tetragonal near-twin (gamma) fitted to the same spectra
    gam = A.analyse(catalog, results, candidates=[GAMMA], minority=False, log=log)
    abl['tetragonal_near_twin'] = {k: v for k, v in payload_from(gam).items() if k != 'minority_phase'}
    # dp5: parallel/crossed settings along the lab axes only -> mirror branch
    lab_only = A.analyse(catalog, results, use=lab_axes_only, candidates=[MAIN], minority=False, log=log)
    Qa = np.array(lab_only['Q'])
    Qm = np.diag([-1.0, 1.0, -1.0]) @ Qa
    sub = [spectra[i] for i in lab_axes_only]
    fitter = A.CandidateFit(e, sub)
    from scipy.spatial.transform import Rotation
    branches = {}
    for name, QQ in (('fit', Qa), ('mirror', Qm)):
        o = fitter.fit(fitter.start(Rotation.from_matrix(QQ).as_rotvec(), lab_only['frequency_scale'], lab_only['axis_offset_cm1']), max_nfev=40)
        pred = A.face_b_predictions(e, o.Q, float(o.x[4]))
        branches[name] = {'chi2': o.chi2, 'pred': pred, 's': float(o.x[4]),
                          'axis_err': M.axis_error_deg(pred['c_axis_lab'], ref['c_axis_lab'])}
    info['lab_axes_only_branches'] = {k: {'chi2': v['chi2'], 'axis_err_deg': v['axis_err']} for k, v in branches.items()}
    wrong = max(branches, key=lambda k: branches[k]['axis_err'])
    right = min(branches, key=lambda k: branches[k]['axis_err'])
    bw = branches[wrong]['pred']
    abl['lab_axes_polarizations_only'] = {'c_axis_lab': bw['c_axis_lab'], 'faceB_strongest_mode_ratio_tt_over_zz': bw['ratio_tt_over_zz'],
                                          'faceB_highest_optical_phonon_cm1': bw['highest_optical_cm1']}
    br = branches[right]['pred']
    info['lab_axes_only_lucky_branch'] = {'axis_err_deg': branches[right]['axis_err'], 'ratio_tt_over_zz': br['ratio_tt_over_zz'],
                                          'highest_optical_cm1': br['highest_optical_cm1']}
    # dp6: no LO mixing in the face-A model (archived TO modes only)
    to = A.analyse(catalog, results, nac=False, candidates=[MAIN], minority=False, log=log)
    to_lo = A.face_b_predictions(e, np.array(to['Q']), to['frequency_scale'])
    abl['to_only_orientation_fit'] = {'c_axis_lab': to['c_axis_lab'], 'faceB_strongest_mode_ratio_tt_over_zz': to['ratio_tt_over_zz'],
                                      'faceB_highest_optical_phonon_cm1': to_lo['highest_optical_cm1']}
    abl['to_only_face_b_phonon'] = {'faceB_highest_optical_phonon_cm1': A.face_b_predictions(e, Q, s, nac=False)['highest_optical_cm1']}
    # dp7: archived Born rows used without the symmetry expansion / screening
    for born in ('unrotated', 'isotropic'):
        abl['born_rows_' + born] = {'faceB_highest_optical_phonon_cm1': A.face_b_predictions(e, Q, s, born=born)['highest_optical_cm1']}
    abl['nac_without_eps_inf'] = {'faceB_highest_optical_phonon_cm1': A.face_b_predictions(e, Q, s, eps='vacuum')['highest_optical_cm1']}
    # dp8: weak extra lines ignored or blamed on LO modes / noise
    abl['minority_phase_ignored'] = {'minority_phase': 'none'}
    # dp9: intensity prefactors omitted in the global fit
    nopref = A.analyse(catalog, results, prefactors=False, candidates=[MAIN], minority=False, log=log)
    info['no_prefactor_fit'] = errors(payload_from(nopref, minority=ALPHA), ref)
    info['no_prefactor_fit']['chi2_per_point'] = nopref['chi2_per_point'][MAIN]

    # robustness: orientation from non-polar line areas only, calibration from the Stokes/anti-Stokes centroid
    d0, nu_main = A.calibrate(spectra)
    s_c = nu_main / A.strongest_above(e)
    rv, _ = A.coarse_orientation(e, spectra, s_c, d0, keep=4)
    cands = [A.area_orientation(e, spectra, s_c, d0, r) for r in rv]
    Qa2 = cands[0]
    pa2 = A.face_b_predictions(e, Qa2, s_c)
    rob['nonpolar_areas_centroid_calibration'] = {'frequency_scale': s_c, 'c_axis_lab': pa2['c_axis_lab'],
                                                  'faceB_strongest_mode_ratio_tt_over_zz': pa2['ratio_tt_over_zz'],
                                                  'faceB_highest_optical_phonon_cm1': pa2['highest_optical_cm1']}
    # robustness: low window only (no high window, no long accumulation)
    low = [i for i, r in enumerate(results) if r['center_cm1'] in (200, -200)]
    lw = A.analyse(catalog, results, use=low, candidates=[MAIN], minority=False, log=log)
    rob['low_window_spectra_only'] = {k: v for k, v in payload_from(lw).items() if k not in ('minority_phase', 'main_phase')}
    # robustness: constant line width (g1 = 0) instead of a width growing with frequency
    fitc = A.CandidateFit(e, spectra)
    lo_b, hi_b = [-np.inf] * 4 + [0.93, -10.0, 0.2, 0.0], [np.inf] * 4 + [1.07, 10.0, 6.0, 1e-6]
    from scipy.optimize import least_squares
    pc = np.array(oracle['fit_params'])
    pc[7] = 5e-7
    oc = least_squares(fitc.resid, pc, bounds=(lo_b, hi_b), x_scale='jac', max_nfev=100)
    Qc = Rotation.from_rotvec(oc.x[:3]).as_matrix()
    pcb = A.face_b_predictions(e, Qc, float(oc.x[4]))
    rob['constant_linewidth_fit'] = {'frequency_scale': float(oc.x[4]), 'c_axis_lab': pcb['c_axis_lab'],
                                     'faceB_strongest_mode_ratio_tt_over_zz': pcb['ratio_tt_over_zz'],
                                     'faceB_highest_optical_phonon_cm1': pcb['highest_optical_cm1']}
    values = {'ablations': abl, 'robustness': rob}
    (task/'solution/ablation_values.json').write_text(json.dumps(values, indent=1) + '\n')

    table = {'oracle': errors(P0, ref)}
    for name, ch in abl.items():
        p = dict(P0)
        p.update(ch)
        table['ablation:' + name] = errors(p, ref)
    for name, ch in rob.items():
        p = dict(P0)
        p.update(ch)
        table['robustness:' + name] = errors(p, ref)
    calib = {'thresholds': thresholds, 'reference': ref, 'oracle': {k: oracle[k] for k in oracle if k not in ('minority_test',)},
             'minority_test': oracle.get('minority_test'), 'table': table, 'info': info}
    (here/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    print('%-46s %9s %8s %9s %9s %10s %10s' % ('path', 's_rel%', 'axis', 'ratio_log', 'lo_rel%', 'main', 'minority'))
    for name, r in table.items():
        print('%-46s %9.3f %8.2f %9.3f %9.3f %10s %10s' % (name, 100 * r['s_rel'], r['axis_deg'], r['ratio_log'],
                                                          100 * r['lo_rel'], r['main'], r['minority']))
    print(json.dumps(info, indent=1, default=float))
