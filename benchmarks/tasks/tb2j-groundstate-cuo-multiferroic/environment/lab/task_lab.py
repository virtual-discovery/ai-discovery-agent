"""Lab logic for tb2j-groundstate-cuo-multiferroic: retrieval of author-archived TB2J output for CuO, and grading.

Purchases return archived isotropic exchange pairs of one of four calculations (FM or AFM reference, with or
without the spin-ligand correction) within a distance window, with the calculation's metadata, or the authors'
archived h(k) eigenvalue files. Nothing graded is served: the targets are derived from the archived pair lists
by the construction code (benchmarks/construction/tb2j-groundstate-cuo-multiferroic).
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


def check_calculation(request, usage):
    if request.get('calculation') not in CALCS:
        raise ValueError('calculation must be one of %s. Use %s' % (', '.join(CALCS), usage))


def validate(request):
    kind = request.get('kind')
    if kind == 'exchange':
        usage = '{"kind": "exchange", "calculation": "FM/ligand", "window": "0-4"}'
        if set(request) != {'kind', 'calculation', 'window'}:
            raise ValueError('Use ' + usage)
        check_calculation(request, usage)
        if request['window'] not in DATA['costs']:
            raise ValueError('window must be one of %s' % ', '.join(DATA['windows']))
        return DATA['costs'][request['window']]
    if kind == 'author_dispersion':
        usage = '{"kind": "author_dispersion", "calculation": "AFM/ligand"}'
        if set(request) != {'kind', 'calculation'}:
            raise ValueError('Use ' + usage)
        check_calculation(request, usage)
        return 1
    raise ValueError('kind must be "exchange" or "author_dispersion"')


def run(request, workdir):
    calc = request['calculation']
    if request['kind'] == 'exchange':
        pairs = DATA['sets'][calc]['windows'][request['window']]
        lo, hi = (float(x) for x in request['window'].split('-'))
        return {'calculation': calc, 'window': request['window'], 'distance_range_A': [lo, hi],
                'metadata': DATA['sets'][calc]['meta'], 'pair_columns': PAIR_COLUMNS, 'n_pairs': len(pairs),
                'pairs': pairs, 'source': SOURCE}
    disp = DATA['dispersions'][calc]
    return {'calculation': calc, 'archived_header': disp['archived_header'], 'n_kpoints': len(disp['values']),
            'n_columns': len(disp['values'][0]), 'values': disp['values'],
            'note': 'Eigenvalues of the bosonic matrix h(k) of the paper (Eq. 9) for this calculation\'s reference '
                    'configuration, one row per k-point along the path in the header, as archived. The k-point '
                    'coordinates and the units are not archived.', 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def vector(v):
    return isinstance(v, list) and len(v) == 3 and all(finite(x) for x in v)


def q_distance(q, ref):
    """Smallest distance (r.l.u.) between q and +-ref + G for integer G."""
    best = float('inf')
    for sgn in (1, -1):
        d = [a - sgn*b for a, b in zip(q, ref)]
        best = min(best, math.sqrt(sum((x - round(x))**2 for x in d)))
    return best


def grade(payload, trace, validation):
    tg, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['kmin_FM_reference_rlu', 'stable_at_own_reference', 'converged_propagation_vector_rlu', 'chain_exchange_meV',
              'mean_magnon_energy_meV', 'magnon_band_maximum_meV', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    km, st = payload['kmin_FM_reference_rlu'], payload['stable_at_own_reference']
    if not isinstance(km, dict) or set(km) != {'ligand_corrected', 'uncorrected'} or not all(vector(v) for v in km.values()):
        return False, {'error': 'kmin_FM_reference_rlu must give three finite numbers for ligand_corrected and uncorrected'}
    if not isinstance(st, dict) or set(st) != set(CALCS) or not all(isinstance(v, bool) for v in st.values()):
        return False, {'error': 'stable_at_own_reference must give true or false for each of ' + ', '.join(CALCS)}
    if not vector(payload['converged_propagation_vector_rlu']):
        return False, {'error': 'converged_propagation_vector_rlu must be three finite numbers'}
    scalars = ['chain_exchange_meV', 'mean_magnon_energy_meV', 'magnon_band_maximum_meV']
    if not all(finite(payload[k]) for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    dist = {k: q_distance(km[k], tg['kmin_FM_reference_rlu'][k]) for k in km}
    dq = q_distance(payload['converged_propagation_vector_rlu'], tg['converged_propagation_vector_rlu'])
    rel = {k: payload[k]/tg[k] - 1 for k in scalars}
    checks = {'evidence_executed': evidence,
              'kmin_ligand_corrected': dist['ligand_corrected'] <= th['kmin_ligand_rlu'],
              'kmin_uncorrected': dist['uncorrected'] <= th['kmin_uncorrected_rlu'],
              'stability_verdicts': st == tg['stable_at_own_reference'],
              'converged_propagation_vector': dq <= th['converged_q_rlu'],
              'chain_exchange': abs(rel['chain_exchange_meV']) <= th['chain_J_rel'],
              'mean_magnon_energy': abs(rel['mean_magnon_energy_meV']) <= th['magnon_mean_rel'],
              'magnon_band_maximum': abs(rel['magnon_band_maximum_meV']) <= th['magnon_max_rel']}
    report = {'checks': checks, 'kmin_distance_rlu': dist, 'converged_q_distance_rlu': dq, 'relative_errors': rel,
              'wrong_stability_verdicts': sorted(k for k in st if st[k] != tg['stable_at_own_reference'][k])}
    return all(checks.values()), report
