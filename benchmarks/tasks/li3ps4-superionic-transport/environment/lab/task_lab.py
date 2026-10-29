"""Lab logic for li3ps4-superionic-transport: retrieval of author-archived ML-MD analyses and grading.

Purchases return archived values (conductivities, PS4/P sublattice dynamics, heat capacity, Li MSD, current
breakdown, pair correlations, order parameters). The conductivity cells listed in catalog.json `to_predict` are
never served; they and quantities derived from the full archive grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
WITHHELD = set(PRIVATE['targets']['sigma_S_per_m'])
FUNCTIONALS = ('PBEsol', 'r2SCAN', 'PBE0')
SOURCE = 'retrieved author-archived analysis of the released ML-MD runs (no new simulation)'


def temperature(request, key='T_K'):
    t = request.get(key)
    if isinstance(t, bool) or not isinstance(t, (int, float)) or not math.isfinite(t) or abs(t - round(t)) > 1e-9:
        raise ValueError(f'{key} must be one of the simulated temperatures listed in /app/catalog.json')
    return int(round(t))


def check_keys(request, keys, usage):
    if set(request) != set(keys):
        raise ValueError('Use ' + usage)


def validate(request):
    kind = request.get('kind')
    if kind == 'conductivity':
        usage = ('{"kind": "conductivity", "functional": ..., "initial_phase": ..., "T_K": ..., '
                 '"estimator": "nernst_einstein" or "green_kubo"}')
        check_keys(request, ('kind', 'functional', 'initial_phase', 'T_K', 'estimator'), usage)
        cell = f"{request['functional']}/{request['initial_phase']}/{temperature(request)}"
        if cell in WITHHELD:
            raise ValueError(f'{cell} is withheld for grading (listed in to_predict); it is not served')
        if request['estimator'] not in ('nernst_einstein', 'green_kubo') or f"{cell}/{request['estimator']}" not in DATA['conductivity']:
            raise ValueError('Not available: see requests.conductivity.available in /app/catalog.json. ' + usage)
        return COSTS[request['estimator']]
    if kind == 'sublattice_dynamics':
        check_keys(request, ('kind', 'functional', 'T_K'), '{"kind": "sublattice_dynamics", "functional": ..., "T_K": ...}')
        if f"{request['functional']}/{temperature(request)}" not in DATA['sublattice_dynamics']:
            raise ValueError('Not available: see requests.sublattice_dynamics.available in /app/catalog.json')
        return COSTS[kind]
    if kind == 'heat_capacity':
        check_keys(request, ('kind', 'functional', 'T_low_K', 'T_high_K'),
                   '{"kind": "heat_capacity", "functional": ..., "T_low_K": ..., "T_high_K": ...}')
        key = f"{request['functional']}/{temperature(request, 'T_low_K')}/{temperature(request, 'T_high_K')}"
        if key not in DATA['heat_capacity']:
            raise ValueError('Not available: T_low_K and T_high_K must be consecutive simulated temperatures listed in '
                             'requests.heat_capacity.available')
        return COSTS[kind]
    if kind in ('li_msd', 'li_pair_correlation'):
        check_keys(request, ('kind',), '{"kind": "%s"}' % kind)
        return COSTS[kind]
    if kind == 'current_breakdown':
        check_keys(request, ('kind', 'frame'), '{"kind": "current_breakdown", "frame": "matrix" or "barycentric"}')
        if request['frame'] not in ('matrix', 'barycentric'):
            raise ValueError('frame must be "matrix" or "barycentric"')
        return COSTS[kind]
    if kind == 'order_parameters':
        check_keys(request, ('kind', 'T_K'), '{"kind": "order_parameters", "T_K": 700, 750 or 800}')
        if str(temperature(request)) not in DATA['order_parameters']['samples']:
            raise ValueError('T_K must be 700, 750 or 800')
        return COSTS[kind]
    raise ValueError('kind must be one of: conductivity, sublattice_dynamics, heat_capacity, li_msd, current_breakdown, '
                     'li_pair_correlation, order_parameters')


def run(request, workdir):
    kind = request['kind']
    if kind == 'conductivity':
        cell = f"{request['functional']}/{request['initial_phase']}/{temperature(request)}"
        v = DATA['conductivity'][f"{cell}/{request['estimator']}"]
        return {'functional': request['functional'], 'initial_phase': request['initial_phase'], 'T_K': temperature(request),
                'estimator': request['estimator'], **v, 'units': 'S/m', 'source': SOURCE}
    if kind == 'sublattice_dynamics':
        T = temperature(request)
        return {'functional': request['functional'], 'T_K': T, **DATA['sublattice_dynamics'][f"{request['functional']}/{T}"],
                'note': 'quantities absent here are not archived at this temperature', 'source': SOURCE}
    if kind == 'heat_capacity':
        lo, hi = temperature(request, 'T_low_K'), temperature(request, 'T_high_K')
        return {'functional': request['functional'], 'T_low_K': lo, 'T_high_K': hi,
                **DATA['heat_capacity'][f"{request['functional']}/{lo}/{hi}"], 'source': SOURCE}
    if kind == 'li_msd':
        return {'run': 'ML-PBE0, started from gamma, 750 K', **DATA['li_msd'], 'source': SOURCE}
    if kind == 'current_breakdown':
        cb = DATA['current_breakdown']
        return {'frame': request['frame'], 'tau_ps': cb['tau_ps'], **cb[request['frame']], 'units': 'S/cm', 'source': SOURCE}
    if kind == 'li_pair_correlation':
        return {'run': 'alpha phase, 650 K, 190 ps, matrix frame', **DATA['li_pair_correlation'], 'source': SOURCE}
    T = temperature(request)
    op = DATA['order_parameters']
    return {'run': f'ML-PBE0, started from gamma, {T} K', 's1_s2_samples': op['samples'][str(T)],
            'reference_clusters_below_transition': op['reference_clusters'], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    tg, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['transition_temperature_K', 'sigma_predictions_S_per_m', 'haven_ratio_PBE0_gamma_750K',
              'activation_energy_eV_PBE0_gamma', 'li_s_cross_fraction', 'log10_hop_to_flip_rate_ratio_PBE0_800K',
              'paddle_wheel_supported', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    tc, sig, ea = payload['transition_temperature_K'], payload['sigma_predictions_S_per_m'], payload['activation_energy_eV_PBE0_gamma']
    if not isinstance(tc, dict) or set(tc) != {'PBE0', 'r2SCAN'} or not all(finite(v) for v in tc.values()):
        return False, {'error': 'transition_temperature_K must give finite numbers for PBE0 and r2SCAN'}
    if not isinstance(sig, dict) or set(sig) != set(tg['sigma_S_per_m']) or not all(finite(v) and v > 0 for v in sig.values()):
        return False, {'error': 'sigma_predictions_S_per_m must give a positive number for every to_predict cell'}
    if not isinstance(ea, dict) or set(ea) != {'below_transition', 'superionic_to_1000K'} or not all(finite(v) for v in ea.values()):
        return False, {'error': 'activation_energy_eV_PBE0_gamma must give finite below_transition and superionic_to_1000K'}
    scalars = ['haven_ratio_PBE0_gamma_750K', 'li_s_cross_fraction', 'log10_hop_to_flip_rate_ratio_PBE0_800K']
    if not all(finite(payload[k]) for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    if not isinstance(payload['paddle_wheel_supported'], bool):
        return False, {'error': 'paddle_wheel_supported must be true or false'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    dex = {k: math.log10(sig[k] / v) for k, v in tg['sigma_S_per_m'].items()}
    errors = {'transition_K': {f: tc[f] - tg['transition_temperature_K'][f] for f in tc},
              'sigma_log10': dex,
              'haven_ratio': payload['haven_ratio_PBE0_gamma_750K'] - tg['haven_ratio_PBE0_gamma_750K'],
              'ea_below_eV': ea['below_transition'] - tg['activation_energy_eV']['below_transition'],
              'ea_above_eV': ea['superionic_to_1000K'] - tg['activation_energy_eV']['superionic_to_1000K'],
              'li_s_fraction': payload['li_s_cross_fraction'] - tg['li_s_cross_fraction'],
              'log10_rate_ratio': payload['log10_hop_to_flip_rate_ratio_PBE0_800K'] - tg['log10_hop_to_flip_rate_ratio_PBE0_800K']}
    checks = {'evidence_executed': evidence,
              'transition_PBE0': abs(errors['transition_K']['PBE0']) <= th['transition_K'],
              'transition_r2SCAN': abs(errors['transition_K']['r2SCAN']) <= th['transition_K'],
              **{f'sigma_{k}': abs(v) <= th['sigma_log10'] for k, v in dex.items()},
              'haven_ratio': abs(errors['haven_ratio']) <= th['haven_ratio'],
              'activation_energy_below_transition': abs(errors['ea_below_eV']) <= th['ea_below_eV'],
              'activation_energy_superionic': abs(errors['ea_above_eV']) <= th['ea_above_eV'],
              'li_s_cross_fraction': abs(errors['li_s_fraction']) <= th['li_s_fraction'],
              'hop_to_flip_rate_ratio': abs(errors['log10_rate_ratio']) <= th['log10_rate_ratio'],
              'paddle_wheel_verdict': payload['paddle_wheel_supported'] is tg['paddle_wheel_supported']}
    return all(checks.values()), {'checks': checks, 'errors': errors}
