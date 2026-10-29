"""Write design.json (with the calibration table) from calibration.json and check.json.

    python3 -B write_design.py
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/raman-polymorph-orientation'
C = json.loads((HERE/'calibration.json').read_text())
T = C['table']
check = json.loads((HERE/'check.json').read_text())['summary'] if (HERE/'check.json').exists() else {}


def r(x, n=4):
    return float('%.*g' % (n, x))


def col(field, scale=1.0):
    out = {'reference_oracle': r(scale * T['oracle'][field])}
    out['robustness'] = {k.split(':', 1)[1]: r(scale * v[field]) for k, v in T.items() if k.startswith('robustness:')}
    abl = {k.split(':', 1)[1]: r(scale * v[field]) for k, v in T.items() if k.startswith('ablation:') and abs(v[field] - T['oracle'][field]) > 1e-9}
    out['ablations'] = abl
    return out


th = C['thresholds']
calibration = {
    'frequency_scale (relative error, %)': {**col('s_rel', 100), 'reference_value': C['reference']['frequency_scale'],
                                            'tolerance': '|s/s_ref - 1| <= %.2f %%' % (100 * th['scale_rel'])},
    'c_axis (angle, deg; modulo sign and the 180-degree turn about the face normal)': {
        **col('axis_deg'), 'tolerance': '<= %.0f deg' % th['axis_deg']},
    'faceB_ratio_tt_over_zz (log error)': {**col('ratio_log'), 'reference_value': r(C['reference']['faceB_ratio_tt_over_zz']),
                                           'tolerance': '|ln(ratio/ref)| <= ln(1.10) = %.4f' % th['ratio_log']},
    'faceB_highest_optical_phonon (relative error, %)': {**col('lo_rel', 100), 'reference_value': r(C['reference']['faceB_highest_optical_cm1'], 6),
                                                         'tolerance': '|f/f_ref - 1| <= %.1f %%' % (100 * th['lo_rel'])},
    'main_phase': {'reference': C['reference']['main_phase'], 'oracle_chi2_per_point': {k: r(v) for k, v in C['oracle']['chi2_per_point'].items()},
                   'ablations': {'nearest_line_peak_matching': 'mp-2843', 'tetragonal_near_twin': 'mp-7457'}, 'tolerance': 'exact'},
    'minority_phase': {'reference': C['reference']['minority_phase'],
                       'oracle_chi2_drop_fraction': {k: r(v['chi2_drop_fraction']) for k, v in C['minority_test'].items()},
                       'decision_rule': 'largest fractional chi2 drop above 5 %', 'ablations': {'minority_phase_ignored': 'none'},
                       'tolerance': 'exact'},
}

design = {
    'task_id': 'raman-polymorph-orientation',
    'catalog_id': 'raman-polymorph-orientation',
    'start': 'released',
    'source_group': 'bagheri_2023_raman_db',
    'title': 'SnF2 crystal: polymorph, orientation and unmeasured geometries from polarized Raman',
    'decision_points': [
        {'id': 'dp1', 'fork': 'Match raw database (PBEsol) frequencies or fit a global frequency scale',
         'wrong_choice': 'Raw frequencies (s = 1)',
         'consequence': 's off by -2.3 %; face-B phonon off by -2.3 %; raw nearest-line matching also prefers the line-dense alpha record',
         'ablation': 'raw_dft_frequencies', 'difficulty_type': 'convention'},
        {'id': 'dp2', 'fork': 'Pure multiplicative scale through the origin, or a scale plus an axis offset revealed by the Stokes/anti-Stokes asymmetry (and the asymmetric notch edges)',
         'wrong_choice': 'Fit s through the origin to the Stokes lines',
         'consequence': 'The +2.12 cm-1 calibration offset of the spectrometer axis is absorbed into s: +0.75 % (tolerance 0.2 %), face-B phonon +0.94 %',
         'ablation': 'scale_without_axis_offset', 'difficulty_type': 'confounder discovered in the raw spectra'},
        {'id': 'dp3', 'fork': 'Compare candidates by peak positions alone or by a forward model of intensities and polarization selection rules',
         'wrong_choice': 'Nearest-line matching of observed peaks to each candidate (with or without a scale)',
         'consequence': 'The monoclinic alpha record (69 Raman lines) wins by line density (mean deviation 3.5 vs 6.1 cm-1); its forward-model chi2/point is 141 against 2.6',
         'ablation': 'nearest_line_peak_matching', 'difficulty_type': 'which analysis discriminates between competing explanations'},
        {'id': 'dp4', 'fork': 'Orthorhombic beta versus its tetragonal near-twin gamma, whose strong lines nearly coincide',
         'wrong_choice': 'Accept gamma (higher symmetry, same main lines) without resolving the split doublets, low-frequency modes and in-plane anisotropy',
         'consequence': 'Wrong polymorph; s -0.24 %, face-B phonon +1.9 %, ratio +13.7 % (gamma chi2/point 57 against 2.6)',
         'ablation': 'tetragonal_near_twin', 'difficulty_type': 'near-degenerate competing explanation; resolution and polarization evidence'},
        {'id': 'dp5', 'fork': 'Which polarizer/analyzer settings to buy: parallel and crossed along the lab axes only, or rotated settings as well',
         'wrong_choice': 'Lab-axis settings only (xx, yy, xy)',
         'consequence': 'A mirror-related orientation fits identically (chi2 equal to 1e-9); the wrong branch puts the c axis 65 deg off, the ratio off by -35 % and the face-B phonon off by +4.9 %. The consequence is visible only after buying rotated settings',
         'ablation': 'lab_axes_polarizations_only', 'difficulty_type': 'discriminating experiment under budget'},
        {'id': 'dp6', 'fork': 'Model the face-A spectra with the archived TO modes, or LO-mix the polar modes for q along the face normal (Born charges, eps_inf)',
         'wrong_choice': 'TO modes only, both in the orientation fit and for the face-B phonon',
         'consequence': 'Oblique LO lines (e.g. near 404 and 472 cm-1 on the axis) are unexplained; the orientation fit is biased by 12 deg (ratio +24 %); the face-B phonon taken from the TO list is -8.4 % off',
         'ablation': 'to_only_orientation_fit', 'difficulty_type': 'evidence in the data (lines absent from the TO list); cascade into the orientation'},
        {'id': 'dp7', 'fork': 'How to use the archived Born charges (3 rows for 12 atoms in beta: one per symmetry-independent atom)',
         'wrong_choice': 'Assign each row unrotated to all atoms of its orbit, use isotropic charges, or omit eps_inf screening',
         'consequence': 'Face-B phonon +5.7 % (unrotated; acoustic sum rule violated), -7.9 % (isotropic), +56 % (no eps_inf)',
         'ablation': 'born_rows_unrotated', 'difficulty_type': 'released data whose shape and convention must be discovered'},
        {'id': 'dp8', 'fork': 'Weak extra lines: a minority phase, LO-mixed lines of the main phase, or noise',
         'wrong_choice': 'Ignore the weak lines or attribute them all to LO mixing',
         'consequence': 'The randomly oriented alpha inclusions (3 % by volume; strongest line near 364 cm-1 on the axis, depolarized identically in every setting) are missed',
         'ablation': 'minority_phase_ignored', 'difficulty_type': 'confounder that mimics the target effect'},
    ],
    'graded_fields': [
        'main_phase (exact mpid)',
        'minority_phase (exact mpid or "none")',
        'frequency_scale within 0.2 %',
        'c_axis_lab within 6 deg (modulo sign and the face-normal 180-degree turn)',
        'faceB_strongest_mode_ratio_tt_over_zz within a factor 1.10',
        'faceB_highest_optical_phonon_cm1 within 0.8 %',
    ],
    'length_estimate': {
        'oracle_lab_evaluations': 15, 'oracle_charge': 17, 'budget': 20,
        'catalog': 'unbounded (any polarizer/analyzer angle, window and accumulation count)',
        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2-4',
        'why': ('Parsing three database records, writing a Raman forward model (tensor rotation, thermal factors, line shapes), '
                'a Stokes/anti-Stokes calibration, a symmetry expansion of the Born charges and an LO-mixing calculation, '
                'planning polarization settings and windows under a 20-unit budget, a global orientation search with multi-start '
                'fits per candidate, and a minority-phase test.'),
        'oracle_container_seconds': (check.get('oracle') or {}).get('container_seconds'),
    },
    'runtime_dependencies': ['numpy, scipy (runtime image); no extra packages'],
    'data_provenance': {
        'archive_url': 'https://archive.materialscloud.org/records/0m9r0-g3486',
        'doi': '10.24435/materialscloud:pg-h3',
        'files_used': ['Computational_Raman_Database_v1.0.1.zip -> crd.json records 3216, 2564, 2390'],
        'md5': {'Computational_Raman_Database_v1.0.1.zip': 'e2aaea03df70fda95dbd970c60456882'},
        'record_sha256': {'3216': 'da176f4d3a76e91f399bd8ca8a1404f0d92c238d31555169d5b961917c860d4c',
                          '2564': '849d4bfea4de74b2e6022c51bfbd4bf999ee1e1bccea4ccd59f3465c2999c376',
                          '2390': 'e5735d7dc160445293f7135b0afd19a52d10710ce85306e1a26d42a6b06f7193'},
        'license': 'CC BY 4.0',
        'synthetic_measurements': 'forward model documented in SOURCE.md; hidden sample parameters in environment/lab/private/hidden.json',
    },
    'exceptions': {},
    'notes': ('One composition (SnF2, three database polymorphs) instead of the design\'s two compositions. No IR reflectance '
              'endpoint. The design\'s degenerate-partner and intensity-prefactor forks are not graded decision points here: '
              'the hidden beta phase has no degenerate modes, and a fit without the (n+1)/omega factors is visibly poor '
              '(chi2/point 84 against 2.6) but still recovers the graded values (c axis 1.1 deg, s +0.02 %). Both '
              'remain rubric items. The design\'s quality-flag fork is a rubric item as well: the flags do not decide the phase.'),
    'calibration': calibration,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('wrote design.json')
