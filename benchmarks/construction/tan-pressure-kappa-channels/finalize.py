"""Write private/hidden.json, hidden_values.json, solution/ablation_values.json and calibration.json from
calibration_raw.json (written by reference.py). Pure JSON; run on the host with python3 -B finalize.py."""
import json
from pathlib import Path

HERE = Path(__file__).parent
TASK = HERE.parents[2]/'benchmarks/tasks/tan-pressure-kappa-channels'
S = json.loads((HERE/'calibration_raw.json').read_text())
FIELDS = ['kappa_xx_80GPa_300K_W_mK', 'kappa_zz_80GPa_300K_W_mK', 'kappa_xx_40GPa_300K_W_mK',
          'kappa_xx_60GPa_600K_W_mK', 'phel_resistivity_fraction_80GPa_300K']
TOL = {'kappa_xx_80GPa_300K_W_mK': {'type': 'rel', 'value': 0.10},
       'kappa_zz_80GPa_300K_W_mK': {'type': 'rel', 'value': 0.10},
       'kappa_xx_40GPa_300K_W_mK': {'type': 'rel', 'value': 0.05},
       'kappa_xx_60GPa_600K_W_mK': {'type': 'rel', 'value': 0.04},
       'phel_resistivity_fraction_80GPa_300K': {'type': 'abs', 'value': 0.06}}
DPS = [('dp1', 'abl_RTA'), ('dp2', 'abl_no_4ph'), ('dp3', 'abl_no_phel'), ('dp4', 'abl_phel_times_2pi'),
       ('dp5', 'abl_phel_row_order'), ('dp6', 'abl_scalebroad_1'), ('dp7', 'abl_kappa_level_interpolation'),
       ('dp8', 'abl_4ph_not_rescaled_to_600K'), ('dp8', 'abl_4ph_classical_T2_scaling'), ('dp9', 'abl_trace_average')]
ref = {f: round(S['reference'][f], 6) for f in FIELDS}
# The quadratic-in-pressure variant is defensible for 80 GPa (interpolation between 40/60/100) but at 40 GPa it is
# an extrapolation-like fit through 0/20/60 that overshoots by 8%; the robustness variant therefore keeps the
# quadratic 80 GPa / 600 K values with the linear 40 GPa value.
quad = dict(S['robust_quadratic_and_T_ratio_from_100GPa'])
quad['kappa_xx_40GPa_300K_W_mK'] = S['oracle']['kappa_xx_40GPa_300K_W_mK']
robust = {'geometric_interpolation': S['robust_geometric_interpolation'], 'quadratic_80GPa_T_ratio_100GPa': quad,
          'T_ratio_from_20GPa_only': S['robust_T_ratio_from_20GPa']}
abl = {name: {f: S[name][f] for f in FIELDS} for _, name in DPS}
robust = {k: {f: v[f] for f in FIELDS} for k, v in robust.items()}
(TASK/'environment/lab/private/hidden.json').write_text(json.dumps({'reference': ref, 'tolerance': TOL}, indent=1) + '\n')
(TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(list(ref.values())) + '\n')
(TASK/'solution/ablation_values.json').write_text(json.dumps(
    {'ablation_decision_points': DPS, 'ablations': abl, 'robustness': robust}, indent=1) + '\n')


def err(f, v):
    return (v / ref[f] - 1) if TOL[f]['type'] == 'rel' else (v - ref[f])


cal = {}
for f in FIELDS:
    cal[f] = {'reference': ref[f], 'tolerance': TOL[f],
              'oracle_route': round(S['oracle'][f], 4), 'oracle_error': round(err(f, S['oracle'][f]), 4),
              'robustness_errors': {k: round(err(f, v[f]), 4) for k, v in robust.items()},
              'ablation_errors': {n: round(err(f, abl[n][f]), 4) for _, n in DPS}}
cal['mc_seed_spread_kappa_xx_80GPa'] = [S['mc_seed1_kappa_xx_80'], S['mc_seed2_kappa_xx_80']]
cal['released_full_physics_diag_W_mK'] = S['context']
(HERE/'calibration.json').write_text(json.dumps(cal, indent=1) + '\n')
for f in FIELDS:
    c = cal[f]
    fails = {n: e for n, e in c['ablation_errors'].items() if abs(e) > TOL[f]['value']}
    print(f, ref[f], 'oracle', c['oracle_error'], 'robust', c['robustness_errors'])
    print('   passing ablations on this field:', {n: e for n, e in c['ablation_errors'].items() if n not in fails})
for _, n in DPS:
    print(n, 'fails' if any(abs(err(f, abl[n][f])) > TOL[f]['value'] for f in FIELDS) else 'PASSES (bad)')
