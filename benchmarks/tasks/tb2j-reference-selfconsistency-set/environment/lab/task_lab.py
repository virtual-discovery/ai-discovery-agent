"""Lab logic for tb2j-reference-selfconsistency-set: retrieval of author-archived TB2J output and grading.

Purchases return the archived isotropic exchange pairs of one of ten calculations (NiO, FeCl2, FePS3 with DFT and
DFT+U, FeP; FM reference or the authors' updated AFM reference) within one window of the archived pair distance,
together with the calculation's metadata. Nothing graded is served: the targets are derived from the archived pair
lists by the construction code (benchmarks/construction/tb2j-reference-selfconsistency-set).
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
CALCS = tuple(DATA['sets'])
SOURCE = 'retrieved author-archived TB2J output (no new calculation)'
PAIR_COLUMNS = ['i', 'j', 'Rx', 'Ry', 'Rz', 'J_iso_meV', 'dx_A', 'dy_A', 'dz_A', 'distance_A']
MATERIALS = ('NiO', 'FeCl2', 'FePS3_DFT', 'FePS3_DFT+U', 'FeP')
SOURCES = ('fm_reference', 'converged_reference')
USAGE = '{"kind": "exchange", "calculation": "FeP/FM", "window": "0-5"}'


def validate(request):
    if request.get('kind') != 'exchange':
        raise ValueError('kind must be "exchange". Use ' + USAGE)
    if set(request) != {'kind', 'calculation', 'window'}:
        raise ValueError('Use ' + USAGE)
    if request['calculation'] not in CALCS:
        raise ValueError('calculation must be one of %s' % ', '.join(CALCS))
    if request['window'] not in DATA['costs']:
        raise ValueError('window must be one of %s' % ', '.join(DATA['costs']))
    return DATA['costs'][request['window']]


def run(request, workdir):
    calc, window = request['calculation'], request['window']
    pairs = DATA['sets'][calc]['windows'][window]
    lo, hi = (float(x) for x in window.split('-'))
    return {'calculation': calc, 'window': window, 'distance_range_A': [lo, hi], 'metadata': DATA['sets'][calc]['meta'],
            'pair_columns': PAIR_COLUMNS, 'n_pairs': len(pairs), 'pairs': pairs, 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def vector(v):
    return isinstance(v, list) and len(v) == 3 and all(finite(x) for x in v)


def qdist(q, ref):
    """Distance (r.l.u.) after reducing components modulo 1 and taking absolute values (sign and mirror equivalence)."""
    return math.sqrt(sum((abs(a - round(a)) - abs(b - round(b)))**2 for a, b in zip(q, ref)))


def grade(payload, trace, validation):
    tg, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['fm_reference_kmin_FeP_rlu', 'converged_exchange_kmin_FeP_rlu', 'stabilization_energy_meV_per_ion', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    if not vector(payload['fm_reference_kmin_FeP_rlu']) or not vector(payload['converged_exchange_kmin_FeP_rlu']):
        return False, {'error': 'the k_min fields must each be three finite numbers'}
    en = payload['stabilization_energy_meV_per_ion']
    if (not isinstance(en, dict) or set(en) != set(MATERIALS)
            or not all(isinstance(en[m], dict) and set(en[m]) == set(SOURCES) for m in MATERIALS)
            or not all(en[m][s] is None or finite(en[m][s]) for m in MATERIALS for s in SOURCES)):
        return False, {'error': 'stabilization_energy_meV_per_ion must give fm_reference and converged_reference '
                                '(a finite number or null) for each of ' + ', '.join(MATERIALS)}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    d_fm = qdist(payload['fm_reference_kmin_FeP_rlu'], tg['fm_reference_kmin_FeP_rlu'])
    d_cv = qdist(payload['converged_exchange_kmin_FeP_rlu'], tg['converged_exchange_kmin_FeP_rlu'])
    energy_ok, energy_err = {}, {}
    for m in MATERIALS:
        for s in SOURCES:
            ref, val = tg['stabilization_energy_meV_per_ion'][m][s], en[m][s]
            key = m + '.' + s
            if ref is None:
                energy_ok[key] = val is None
                energy_err[key] = 'expected null' if val is not None else 0.0
            else:
                energy_ok[key] = val is not None and abs(val/ref - 1) <= th['energy_rel']
                energy_err[key] = None if val is None else val/ref - 1
    checks = {'evidence_executed': evidence,
              'fm_reference_kmin_FeP': d_fm <= th['kmin_fm_FeP_rlu'],
              'converged_exchange_kmin_FeP': d_cv <= th['kmin_converged_FeP_rlu'],
              'stabilization_energies': all(energy_ok.values())}
    report = {'checks': checks, 'kmin_distance_rlu': {'fm_reference': d_fm, 'converged_exchange': d_cv},
              'stabilization_relative_errors': energy_err,
              'failed_stabilization_entries': sorted(k for k, v in energy_ok.items() if not v)}
    return all(checks.values()), report
