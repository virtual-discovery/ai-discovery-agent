"""Lab logic for fluorite-superionic-finite-size: retrieval of author-released MD runs and grading.

Purchases return the released input deck, MSD table (subset of lags) and a thermodynamic summary of one
run. The runs listed under `predict` in the catalog (and the other 6x6x6/7x7x7 PbF2 and 6x6x6 AgI runs)
are never served; their diffusion coefficients, and reference values derived from the released data,
grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
MATERIALS = ('CaF2', 'PbF2', 'UO2', 'AgI')
USAGE = 'Use {"kind": "md_run", "material": "CaF2"|"PbF2"|"UO2"|"AgI", "supercell": <n>, "temperature_K": <T>}'


def key_of(request):
    if request.get('kind') != 'md_run' or set(request) != {'kind', 'material', 'supercell', 'temperature_K'}:
        raise ValueError(USAGE)
    m, n, T = request['material'], request['supercell'], request['temperature_K']
    if m not in MATERIALS or type(n) is not int or type(T) is not int:
        raise ValueError(USAGE)
    key = f'{m}-{n}-{T}'
    if f'{m}-{n}x{n}x{n}-{T}K' in PRIVATE['diffusion_targets_A2_per_ps']:
        raise ValueError('This run is withheld: it is one of the runs you must predict')
    if key not in DATA:
        raise ValueError('No released run for this material, supercell and temperature (see runs_available)')
    return key


def validate(request):
    return max(1, math.ceil(DATA[key_of(request)]['atoms'] / 300))


def run(request, workdir):
    return DATA[key_of(request)]


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t = PRIVATE['thresholds']
    fields = ['transition_temperature_K', 'activation_energy_eV', 'diffusion_predictions_A2_per_ps', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    tc_ref, ea_ref, d_ref = (PRIVATE['transition_temperature_K'], PRIVATE['activation_energy_eV'],
                             PRIVATE['diffusion_targets_A2_per_ps'])
    tc, ea, dp = payload['transition_temperature_K'], payload['activation_energy_eV'], payload['diffusion_predictions_A2_per_ps']
    if not isinstance(tc, dict) or set(tc) != set(tc_ref) or not all(finite(v) for v in tc.values()):
        return False, {'error': 'transition_temperature_K must give a finite T_c for CaF2, PbF2 and UO2'}
    if not isinstance(ea, dict) or set(ea) != set(ea_ref) or not all(
            isinstance(ea[m], dict) and set(ea[m]) == set(ea_ref[m]) and all(finite(v) for v in ea[m].values()) for m in ea_ref):
        return False, {'error': 'activation_energy_eV must give finite ordered and superionic values for CaF2 and PbF2'}
    if not isinstance(dp, dict) or set(dp) != set(d_ref) or not all(finite(v) for v in dp.values()):
        return False, {'error': 'diffusion_predictions_A2_per_ps must give a finite D for every id in predict'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks = {'evidence_executed': evidence}
    report = {'Tc_error_K': {}, 'Ea_error_eV': {}, 'log10_D_error': {}}
    for m, ref in tc_ref.items():
        report['Tc_error_K'][m] = tc[m] - ref
        checks[f'Tc_{m}'] = abs(tc[m] - ref) <= t['Tc_tol_K']
    for m in ea_ref:
        for regime, ref in ea_ref[m].items():
            report['Ea_error_eV'][f'{m}_{regime}'] = ea[m][regime] - ref
            checks[f'Ea_{m}_{regime}'] = abs(ea[m][regime] - ref) <= t['Ea_tol_eV'][m][regime]
    for k, ref in d_ref.items():
        err = math.log10(dp[k] / ref) if dp[k] > 0 else None
        report['log10_D_error'][k] = err
        checks[f'D_{k}'] = err is not None and abs(err) <= t['logD_tol_dex']
    report['checks'] = checks
    return all(checks.values()), report
