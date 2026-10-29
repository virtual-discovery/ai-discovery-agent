"""Lab logic for kagome135-correlation-origin: slave-spin mean-field (SSMF) evaluations of the archived
135-kagome Wannier Hamiltonians, and grading.

Each request runs kag_ssmf.lab_ssmf on one archived model (Zenodo 17951035) at a requested electron count,
U and J/U. Graded values are never served directly: they are derived from several evaluations and from the
models themselves (private/hidden.json holds them).
"""
import json
import math
from pathlib import Path

import kag_ssmf as K

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
MODELS = ['TiSb', 'VSb', 'CrSb', 'TiBi', 'CrBi', 'CrBi_9.2']
ORBITALS = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
FORMAT = ('Use {"kind": "ssmf", "model": <one of %s>, "electrons": <26-36>, "U_eV": <0-8>, "J_over_U": <0-0.3>}'
          % MODELS)


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def validate(request):
    if request.get('kind') != 'ssmf' or set(request) != {'kind', 'model', 'electrons', 'U_eV', 'J_over_U'}:
        raise ValueError(FORMAT)
    if request['model'] not in MODELS:
        raise ValueError(FORMAT)
    if not (number(request['electrons']) and 26 <= request['electrons'] <= 36):
        raise ValueError('electrons must be a number between 26 and 36')
    if not (number(request['U_eV']) and 0 <= request['U_eV'] <= 8):
        raise ValueError('U_eV must be a number between 0 and 8')
    if not (number(request['J_over_U']) and 0 <= request['J_over_U'] <= 0.3):
        raise ValueError('J_over_U must be a number between 0 and 0.3')
    return 1


def run(request, workdir):
    model = K.load_model(HERE/'data'/(request['model'] + '.npz'))
    r = K.lab_ssmf(model, float(request['electrons']), float(request['U_eV']), float(request['J_over_U']))
    out = {'model': request['model'], 'electrons': request['electrons'], 'U_eV': request['U_eV'],
           'J_eV': request['U_eV'] * request['J_over_U'], 'metallic_solution': r['metallic'],
           'source': 'new SSMF evaluation of the archived Wannier Hamiltonian (lab evaluator, 9x9x6 k mesh)'}
    if r['metallic']:
        out['Z_by_wannier_function'] = {str(5 * s + i + 1): round(r['Z'][i], 6) for s in range(3) for i in range(5)}
        out['occupation_per_spin_by_wannier_function'] = {str(5 * s + i + 1): round(r['n'][i], 6)
                                                         for s in range(3) for i in range(5)}
        out['chemical_potential_eV'] = round(r['mu'], 6)
    else:
        out['note'] = ('no converged metallic solution could be continued from weak coupling to this U '
                       '(Mott-localized or beyond the metallic branch)')
    return out


def close(value, target, abs_tol=None, rel_tol=None):
    if not number(value):
        return False
    tol = abs_tol if abs_tol is not None else rel_tol * abs(target)
    return abs(value - target) <= tol


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerances']
    fields = ['electron_counts', 'CsCr3Sb5_U6', 'most_renormalized_orbital', 'filling_share', 'hund_ratio',
              'CsCr3Bi5_U4p5', 'sigma_eV', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    counts = payload['electron_counts']
    counts_ok = (isinstance(counts, dict) and set(counts) == set(MODELS)
                 and all(number(counts[m]) and abs(counts[m] - ref['electron_counts'][m]) < 0.5 for m in MODELS))
    cr = payload['CsCr3Sb5_U6'] if isinstance(payload['CsCr3Sb5_U6'], dict) else {}
    bi = payload['CsCr3Bi5_U4p5'] if isinstance(payload['CsCr3Bi5_U4p5'], dict) else {}
    sg = payload['sigma_eV'] if isinstance(payload['sigma_eV'], dict) else {}
    checks = {
        'evidence_executed': evidence,
        'electron_counts': counts_ok,
        'CsCr3Sb5_Z_mean': close(cr.get('Z_mean'), ref['CrSb_Z_mean'], abs_tol=tol['Z_abs']),
        'CsCr3Sb5_Z_min': close(cr.get('Z_min'), ref['CrSb_Z_min'], abs_tol=tol['Z_abs']),
        'most_renormalized_orbital': payload['most_renormalized_orbital'] == ref['most_renormalized_orbital'],
        'filling_share': close(payload['filling_share'], ref['filling_share'], abs_tol=tol['share_abs']),
        'hund_ratio': close(payload['hund_ratio'], ref['hund_ratio'], rel_tol=tol['hund_rel']),
        'CsCr3Bi5_Sb_lattice_Z_mean': close(bi.get('Sb_lattice_Z_mean'), ref['CrBi_U4p5_Z_mean'],
                                            abs_tol=tol['Z_abs']),
        'CsCr3Bi5_TiBi_lattice_Z_mean': close(bi.get('TiBi_lattice_Z_mean'), ref['CrBi92_U4p5_Z_mean'],
                                              abs_tol=tol['Z_abs']),
        'sigma_CsCr3Bi5': close(sg.get('CsCr3Bi5_Sb_lattice'), ref['sigma_CrBi'], rel_tol=tol['sigma_rel']),
        'sigma_CsTi3Sb5': close(sg.get('CsTi3Sb5'), ref['sigma_TiSb'], rel_tol=tol['sigma_rel']),
    }
    return all(checks.values()), {'checks': checks}
