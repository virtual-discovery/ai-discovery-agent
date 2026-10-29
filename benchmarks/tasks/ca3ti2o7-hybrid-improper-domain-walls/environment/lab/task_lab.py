"""Lab logic for ca3ti2o7-hybrid-improper-domain-walls.

The lab evaluates agent-built Ca3Ti2O7 structures with the authors' released neuroevolution
potentials (Skogvoll et al. 2026, Zenodo 22801877), using a numpy implementation of NEP/qNEP
(nep_numpy.py): single points, relaxations (fixed cell or cell length along one direction) and
Born effective charges. The reference domain-wall energies are never computed or served here;
they live in private/hidden.json and grade the submission.
"""
import json
import math
from pathlib import Path
import numpy as np

from nep_numpy import NEP
from relax import Relaxer

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
POTENTIALS = {'qNEP-w-BECs': 'Potential-qNEP-w-BECs.txt', 'qNEP': 'Potential-qNEP.txt', 'NEP': 'Potential-NEP.txt'}
_MODELS = {}
MAX_ATOMS = 1200
KINDS = ('single_point', 'relax', 'born_charges')


def model(name):
    if name not in _MODELS:
        _MODELS[name] = NEP(str(HERE/'potentials'/POTENTIALS[name]))
    return _MODELS[name]


def parse_structure(s):
    if not isinstance(s, dict) or set(s) != {'symbols', 'positions', 'cell'}:
        raise ValueError('structure must be {"symbols": [...], "positions": [[x, y, z], ...], "cell": [[a], [b], [c]]} (Angstrom)')
    sym = s['symbols']
    try:
        pos = np.array(s['positions'], float)
        cell = np.array(s['cell'], float)
    except (TypeError, ValueError):
        raise ValueError('positions and cell must be numeric arrays') from None
    if not isinstance(sym, list) or not all(x in ('Ca', 'Ti', 'O') for x in sym):
        raise ValueError('symbols must be a list of "Ca", "Ti", "O"')
    if pos.shape != (len(sym), 3) or cell.shape != (3, 3) or not np.all(np.isfinite(pos)) or not np.all(np.isfinite(cell)):
        raise ValueError('positions must be N x 3 and cell 3 x 3, finite')
    n = len(sym)
    if n > MAX_ATOMS:
        raise ValueError('at most %d atoms per structure' % MAX_ATOMS)
    counts = {e: sym.count(e) for e in ('Ca', 'Ti', 'O')}
    nfu = counts['Ti'] // 2
    if counts['Ti'] % 2 or counts != {'Ca': 3 * nfu, 'Ti': 2 * nfu, 'O': 7 * nfu} or nfu == 0:
        raise ValueError('composition must be an integer number of Ca3Ti2O7 formula units')
    vol = abs(np.linalg.det(cell))
    heights = [vol / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3])) for k in range(3)]
    if min(heights) < 3.0:
        raise ValueError('every cell height must be at least 3 Angstrom')
    _, I, J, S = NEP._pairs(pos, cell, 1.2)
    if len(I):
        raise ValueError('atoms closer than 1.2 Angstrom; check the structure')
    return sym, pos, cell, nfu


def validate(request):
    kind = request.get('kind')
    if kind not in KINDS:
        raise ValueError('kind must be one of %s' % (KINDS,))
    allowed = {'single_point': {'kind', 'structure', 'potential'},
               'relax': {'kind', 'structure', 'potential', 'cell', 'normal_axis', 'fmax', 'max_steps'},
               'born_charges': {'kind', 'structure'}}[kind]
    if not set(request) <= allowed or 'structure' not in request:
        raise ValueError('allowed keys for %s: %s' % (kind, sorted(allowed)))
    pot = request.get('potential', 'qNEP-w-BECs')
    if pot not in POTENTIALS:
        raise ValueError('potential must be one of %s' % sorted(POTENTIALS))
    sym, pos, cell, nfu = parse_structure(request['structure'])
    if kind == 'relax':
        if request.get('cell', 'fixed') not in ('fixed', 'normal'):
            raise ValueError('cell must be "fixed" or "normal"')
        if request.get('cell', 'fixed') == 'normal' and request.get('normal_axis') not in (0, 1, 2):
            raise ValueError('cell "normal" needs normal_axis 0, 1 or 2 (index of the cell vector that changes length)')
        fmax = request.get('fmax', 0.02)
        steps = request.get('max_steps', 600)
        if not isinstance(fmax, (int, float)) or isinstance(fmax, bool) or not 0.005 <= fmax <= 0.5:
            raise ValueError('fmax must be between 0.005 and 0.5 eV/A')
        if type(steps) is not int or not 1 <= steps <= 1500:
            raise ValueError('max_steps must be an integer between 1 and 1500')
        per = 800 if pot == 'NEP' else 400
        return int(math.ceil(len(sym) / per))
    return 1


def summary(o, sym, pos, cell, nfu):
    vol = abs(np.linalg.det(cell))
    f = np.sqrt((o['forces'] ** 2).sum(1))
    return {'energy_eV': round(float(o['energy']), 6), 'energy_per_formula_unit_eV': round(float(o['energy']) / nfu, 6),
            'formula_units': nfu, 'max_atomic_force_eV_A': round(float(f.max()), 5),
            'stress_GPa': np.round(-o['virial'] / vol * 160.21766, 4).tolist(),
            'volume_A3': round(float(vol), 4)}


def run(request, workdir):
    sym, pos, cell, nfu = parse_structure(request['structure'])
    kind = request['kind']
    pot = request.get('potential', 'qNEP-w-BECs')
    m = model('qNEP-w-BECs' if kind == 'born_charges' else pot)
    types = [m.elements.index(x) for x in sym]
    if kind == 'single_point':
        o = m.compute(types, pos, cell)
        return {'potential': pot, **summary(o, sym, pos, cell, nfu)}
    if kind == 'born_charges':
        o = m.compute(types, pos, cell, want_bec=True)
        return {'potential': 'qNEP-w-BECs', **summary(o, sym, pos, cell, nfu),
                'charges_e': np.round(o['charges'], 4).tolist(),
                'born_effective_charges_e': np.round(o['bec'], 4).tolist(),
                'convention': 'born_effective_charges_e[i][a][b] = dP_a / du_b for atom i (row = polarization, column = displacement)'}
    axis = request.get('normal_axis') if request.get('cell', 'fixed') == 'normal' else None
    r = Relaxer(m, types, pos, cell, relax_axis=axis)
    res = r.run(fmax=float(request.get('fmax', 0.02)), max_steps=int(request.get('max_steps', 600)))
    o = m.compute(types, res['positions'], res['cell'])
    out = {'potential': pot, 'cell_mode': request.get('cell', 'fixed'), **summary(o, sym, res['positions'], res['cell'], nfu),
           'converged': bool(res['converged']), 'iterations': res['iterations'],
           'normal_strain': round(float(res['strain']), 6),
           'final_structure': {'symbols': sym, 'cell': np.round(res['cell'], 6).tolist(),
                               'positions': np.round(res['positions'], 5).tolist()}}
    return out


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


NUMERIC = ['E_neutral_psi_wall_mJ_m2', 'E_neutral_phi_wall_mJ_m2', 'E_layer_psi_wall_mJ_m2',
           'charged_psi_wall_slope_mJ_m2_per_nm']
CATEGORICAL = {'psi_wall_profile': ('ising', 'rotational'), 'phi_wall_profile': ('ising', 'rotational')}


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerance']
    fields = NUMERIC + list(CATEGORICAL) + ['evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    if not all(finite(payload[f]) for f in NUMERIC):
        return False, {'error': 'numeric fields must be finite numbers'}
    if not all(payload[f] in allowed for f, allowed in CATEGORICAL.items()):
        return False, {'error': 'profile fields must be "ising" or "rotational"'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks = {'evidence_executed': evidence}
    report = {}
    for f in NUMERIC:
        lo, hi = tol[f]
        checks[f] = lo <= payload[f] <= hi
        report[f] = payload[f]
    for f in CATEGORICAL:
        checks[f] = payload[f] == ref[f]
    return all(checks.values()), {'checks': checks, 'values': report}
