"""Fill design.json -> calibration from calibration.json (make_task.py) and the gate's oracle report.

    python3 -B write_design.py      (after make_task.py and check_task)
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TASK = ROOT/'benchmarks/tasks/co3sn2s2-dmi-diffuse'
c = json.loads((HERE/'calibration.json').read_text())
state = ROOT/'runs/check/co3sn2s2-dmi-diffuse/lab-results/state.json'
oracle = json.loads(state.read_text())['submission']['report'] if state.exists() else {}
r3 = lambda x: round(float(x), 3)
rob, abl, ref = c['robustness'], c['ablations'], c['reference_recipe']
ham = c['hamiltonian']
cal = {}
for T in ('170', '175', '200'):
    key = 'chi2r_' + T
    cal['narrow_%s_K_profile_chi2r' % T] = {
        'reference_recipe': r3(ref[key]),
        'oracle': r3(oracle.get('scores', {}).get('narrow_%s_chi2r' % T, float('nan'))),
        'robustness': {k: r3(v[key]) for k, v in rob.items() if key in v},
        'nearest_ablation': {(k if k != 'temperature_independent_background' else
                              k + ' (passes: absorbed by the fit, so not an ablation)'): r3(v[key])
                             for k, v in abl.items() if key in v},
        'tolerance': '<= 12'}
cal['wide_101_180K_profile_chi2r'] = {
    'reference_recipe': r3(ref['chi2r_101']),
    'oracle': r3(oracle.get('scores', {}).get('wide_101_chi2r', float('nan'))),
    'robustness': {k: r3(v['chi2r_101']) for k, v in rob.items() if 'chi2r_101' in v},
    'nearest_ablation': {k: r3(v['chi2r_101']) for k, v in abl.items() if 'chi2r_101' in v},
    'tolerance': '<= 15'}
v = ham['variants']
cal['canting_angle_deg'] = {
    'reference': r3(ham['reference']['tilt_deg']),
    'robustness': {'rounded_to_0.01': round(ham['reference']['tilt_deg'], 2),
                   'no_interlayer_coupling (passes, not an ablation)': r3(v['no_interlayer']['tilt_deg'])},
    'nearest_ablation': {'reversed_generator_bond': r3(v['reversed_generator_bond']['tilt_deg']),
                         'in_plane_dm_only (fails via the gap)': r3(v['in_plane_dm_only']['tilt_deg']),
                         'out_of_plane_dm_only': r3(v['out_of_plane_dm_only']['tilt_deg']),
                         'copy_matrix_to_all_bonds': r3(v['copy_matrix_to_all_bonds']['tilt_deg']),
                         'measured_canting_paper (not graded)': 17.0},
    'tolerance': '+-0.5 deg'}
cal['spin_wave_gap_meV'] = {
    'reference': r3(ham['reference']['gap_meV']),
    'independent_check': {'linearised_landau_lifshitz_k0': r3(ham['landau_lifshitz_k0_meV'][0])},
    'robustness': {'conventional_9_site_cell': r3(ham['reference']['gap_meV']),
                   'no_interlayer_coupling (passes, not an ablation)': r3(v['no_interlayer']['gap_meV'])},
    'nearest_ablation': {'in_plane_dm_only': r3(v['in_plane_dm_only']['gap_meV']),
                         'reversed_generator_bond': r3(v['reversed_generator_bond']['gap_meV']),
                         'spin_0.15': r3(ham['gap_spin_0p15_meV']), 'bonds_counted_twice': r3(2*ham['reference']['gap_meV']),
                         'out_of_plane_dm_only / copy': 0.0, 'measured_gap_literature (not graded)': '2-3'},
    'tolerance': '+-6%'}
design = json.loads((TASK/'design.json').read_text())
design['calibration'] = cal
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps(cal, indent=1))
