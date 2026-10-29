"""Lab logic for gwbse-1d-chain-excitons: retrieval of author-archived DFT, GW and BSE results and grading.

Purchases return archived results of Grillo et al. (Materials Cloud 42t0v-jam38); nothing is recomputed.
The scGW QP table and the BSE_scGW spectrum of Bi2Te3 are never served. The graded quantities are derived
values (gaps, optical gaps, screening shifts, Bi2Te3 binding energy) held in private/hidden.json.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
CHAINS = ['S3', 'Te3', 'As2S3', 'Bi2Te3']
SOURCE = 'retrieved author-archived result (no new calculation)'


def check_chain(request):
    if request.get('chain') not in CHAINS:
        raise ValueError('chain must be one of S3, Te3, As2S3, Bi2Te3')


def validate(request):
    kind = request.get('kind')
    if kind == 'qp':
        if set(request) != {'kind', 'chain', 'level', 'k_indices'}:
            raise ValueError('Use {"kind": "qp", "chain": <chain>, "level": "G0W0" or "scGW", "k_indices": [integers 1-25]}')
        check_chain(request)
        if request['level'] not in ('G0W0', 'scGW'):
            raise ValueError('level must be "G0W0" or "scGW"')
        ks = request['k_indices']
        if (not isinstance(ks, list) or not ks or len(set(ks)) != len(ks)
                or not all(type(k) is int and 1 <= k <= 25 for k in ks)):
            raise ValueError('k_indices must be a non-empty list of distinct integers between 1 and 25')
        if request['level'] not in DATA['qp'][request['chain']]:
            raise ValueError(f"the {request['level']} QP table of {request['chain']} is withheld for grading")
        return COSTS['qp_per_k_point'] * len(ks)
    if kind == 'bse_spectrum':
        if set(request) != {'kind', 'chain', 'run'}:
            raise ValueError('Use {"kind": "bse_spectrum", "chain": <chain>, "run": "BSE_G0W0" or "BSE_scGW"}')
        check_chain(request)
        if request['run'] not in ('BSE_G0W0', 'BSE_scGW'):
            raise ValueError('run must be "BSE_G0W0" or "BSE_scGW"')
        if request['run'] not in DATA['bse'][request['chain']]:
            raise ValueError(f"the {request['run']} spectrum of {request['chain']} is withheld for grading")
        return COSTS['bse_spectrum']
    if kind == 'dft_bands':
        if set(request) != {'kind', 'chain', 'soc'} or not isinstance(request.get('soc'), bool):
            raise ValueError('Use {"kind": "dft_bands", "chain": <chain>, "soc": true or false}')
        check_chain(request)
        return COSTS['dft_bands']
    if kind == 'wire_screening':
        if set(request) != {'kind', 'chain'}:
            raise ValueError('Use {"kind": "wire_screening", "chain": <chain>}')
        check_chain(request)
        return COSTS['wire_screening']
    raise ValueError('kind must be "qp", "bse_spectrum", "dft_bands" or "wire_screening"')


def run(request, workdir):
    chain, kind = request['chain'], request['kind']
    if kind == 'qp':
        table = DATA['qp'][chain][request['level']]
        rows = {str(k): table['rows_by_k_index'][str(k)] for k in sorted(request['k_indices'])}
        return {'chain': chain, 'level': request['level'], 'file': table['file'], 'row_format': table['row_format'],
                'rows_by_k_index': rows, 'header': table['header'], 'source': SOURCE}
    if kind == 'bse_spectrum':
        return {'chain': chain, 'run': request['run'], **DATA['bse'][chain][request['run']], 'source': SOURCE}
    if kind == 'dft_bands':
        return {'chain': chain, **DATA['bands'][chain]['soc' if request['soc'] else 'nosoc'], 'source': SOURCE}
    return {'chain': chain, **DATA['screening'][chain], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['scgw_direct_gap_eV', 'scgw_optical_gap_eV', 'screening_update_shift_eV', 'bi2te3_binding_energy_eV',
              'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    gaps, opt, dw = payload['scgw_direct_gap_eV'], payload['scgw_optical_gap_eV'], payload['screening_update_shift_eV']
    if (not isinstance(gaps, dict) or set(gaps) != set(t['scgw_direct_gap_eV'])
            or not all(isinstance(v, dict) and set(v) == {'with_soc', 'without_soc'} and all(finite(x) for x in v.values())
                       for v in gaps.values())):
        return False, {'error': 'scgw_direct_gap_eV must give finite with_soc and without_soc values for Te3 and Bi2Te3'}
    for name, field in (('scgw_optical_gap_eV', opt), ('screening_update_shift_eV', dw)):
        if not isinstance(field, dict) or set(field) != set(t[name]) or not all(finite(v) for v in field.values()):
            return False, {'error': f'{name} must give a finite value for S3, Te3 and As2S3'}
    if not finite(payload['bi2te3_binding_energy_eV']):
        return False, {'error': 'bi2te3_binding_energy_eV must be a finite number'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for chain, ref in t['scgw_direct_gap_eV'].items():
        for key in ('with_soc', 'without_soc'):
            errors[f'gap:{chain}:{key}'] = gaps[chain][key] - ref[key]
            checks[f'gap:{chain}:{key}'] = abs(errors[f'gap:{chain}:{key}']) <= th['gap_abs_eV']
    for chain, ref in t['scgw_optical_gap_eV'].items():
        errors[f'optical_gap:{chain}'] = opt[chain] - ref
        checks[f'optical_gap:{chain}'] = abs(errors[f'optical_gap:{chain}']) <= th['optical_gap_abs_eV']
    for chain, ref in t['screening_update_shift_eV'].items():
        errors[f'screening_shift:{chain}'] = dw[chain] - ref
        checks[f'screening_shift:{chain}'] = abs(errors[f'screening_shift:{chain}']) <= th['screening_shift_abs_eV']
    low, high = th['bi2te3_binding_window_eV']
    eb = payload['bi2te3_binding_energy_eV']
    errors['bi2te3_binding_energy'] = eb - t['bi2te3_binding_energy_eV']
    checks['bi2te3_binding_energy'] = low <= eb <= high
    return all(checks.values()), {'checks': checks, 'errors_eV': errors}
