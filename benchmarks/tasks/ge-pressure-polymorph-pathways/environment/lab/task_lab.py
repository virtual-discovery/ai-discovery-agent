"""Lab logic for ge-pressure-polymorph-pathways.

Two request kinds:
  record  - one archived VASP (PBE) calculation of the Fantasia et al. (2024) Ge database (Materials Cloud
            238ad-1bg44): settings and every ionic step (energy, stress, cell, positions, forces), as parsed
            from the released OUTCAR. No new calculation.
  dp      - the authors' released DeePMD potential (Fantasy-frz.pb), evaluated by a numpy re-implementation
            (dp_numpy.py) on structures the agent submits: single points, fixed-cell ionic relaxations, or
            cell + ion relaxations at a target pressure.
The graded transition pressures and energy differences live in private/hidden.json and are never served.
"""
import gzip
import json
import math
from pathlib import Path
import numpy as np

from dp_numpy import DP
from dprelax import relax_ions, relax_full

HERE = Path(__file__).parent
RECORDS = json.load(gzip.open(HERE/'records.json.gz', 'rt'))
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
KB_PER_EV_A3 = 1602.1766208
MAX_STRUCTURES, MAX_ATOMS = 20, 64
MODES = ('single', 'relax_ions', 'relax_cell')
_DP = []


def model():
    if not _DP:
        _DP.append(DP(HERE/'potential'))
    return _DP[0]


def parse_structure(s):
    if not isinstance(s, dict) or not {'cell', 'positions'} <= set(s) or not set(s) <= {'cell', 'positions', 'mode', 'pressure_GPa'}:
        raise ValueError('each structure is {"cell": 3x3 (A, rows = lattice vectors), "positions": Nx3 Cartesian A, '
                         '"mode": "single" | "relax_ions" | "relax_cell", "pressure_GPa": P (relax_cell only)}')
    try:
        cell = np.array(s['cell'], float)
        pos = np.array(s['positions'], float)
    except (TypeError, ValueError):
        raise ValueError('cell and positions must be numeric arrays') from None
    if cell.shape != (3, 3) or pos.ndim != 2 or pos.shape[1] != 3 or not 1 <= len(pos) <= MAX_ATOMS:
        raise ValueError('cell must be 3x3 and positions N x 3 with 1 <= N <= %d' % MAX_ATOMS)
    if not np.all(np.isfinite(cell)) or not np.all(np.isfinite(pos)):
        raise ValueError('non-finite numbers in the structure')
    vol = abs(np.linalg.det(cell))
    if vol / len(pos) < 8.0 or vol / len(pos) > 60.0:
        raise ValueError('volume per atom must be between 8 and 60 A^3')
    mode = s.get('mode', 'single')
    if mode not in MODES:
        raise ValueError('mode must be one of %s' % (MODES,))
    p = s.get('pressure_GPa', 0.0)
    if mode == 'relax_cell':
        if not isinstance(p, (int, float)) or isinstance(p, bool) or not -10.0 <= p <= 60.0:
            raise ValueError('relax_cell needs pressure_GPa between -10 and 60')
    elif 'pressure_GPa' in s:
        raise ValueError('pressure_GPa is only used with mode "relax_cell"')
    J, D = DP.neighbours(cell, pos, 1.6, 64)
    if np.any(J >= 0):
        raise ValueError('two atoms are closer than 1.6 A; check the structure')
    return cell, pos, mode, float(p)


def validate(request):
    kind = request.get('kind')
    if kind == 'record':
        if set(request) != {'kind', 'id'} or request['id'] not in RECORDS:
            raise ValueError('Use {"kind": "record", "id": <id from /app/catalog.json>}')
        return 1
    if kind == 'dp':
        if set(request) != {'kind', 'structures'} or not isinstance(request['structures'], list) or \
                not 1 <= len(request['structures']) <= MAX_STRUCTURES:
            raise ValueError('Use {"kind": "dp", "structures": [<1 to %d structures>]}' % MAX_STRUCTURES)
        for s in request['structures']:
            parse_structure(s)
        return 1
    raise ValueError('kind must be "record" or "dp"')


def r6(x):
    return float(round(float(x), 6))


def stress_kb(stress_eV_A3):
    """VASP convention and order: XX YY ZZ XY YZ ZX, positive = compressive (pressure)."""
    s = -np.asarray(stress_eV_A3) * KB_PER_EV_A3
    return [r6(s[0, 0]), r6(s[1, 1]), r6(s[2, 2]), r6(s[0, 1]), r6(s[1, 2]), r6(s[0, 2])]


def run_structure(s):
    cell, pos, mode, p = parse_structure(s)
    dp = model()
    n = len(pos)
    if mode == 'single':
        r = dp.compute(cell, pos)
        return {'mode': mode, 'natoms': n, 'energy_eV': r6(r['energy']), 'volume_A3': r6(abs(np.linalg.det(cell))),
                'stress_kB': stress_kb(r['stress_eV_A3']), 'forces_eV_A': np.round(r['forces'], 5).tolist(),
                'max_force_eV_A': r6(np.linalg.norm(r['forces'], axis=1).max())}
    if mode == 'relax_ions':
        o = relax_ions(dp, cell, pos)
        r = dp.compute(cell, o['positions'])
        return {'mode': mode, 'natoms': n, 'energy_initial_eV': r6(dp.energy(cell, pos)), 'energy_eV': r6(r['energy']),
                'volume_A3': r6(abs(np.linalg.det(cell))), 'stress_kB': stress_kb(r['stress_eV_A3']),
                'max_force_eV_A': r6(np.linalg.norm(r['forces'], axis=1).max()), 'iterations': o['iterations'],
                'positions_A': np.round(o['positions'], 6).tolist(), 'cell_A': np.round(cell, 6).tolist()}
    o = relax_full(dp, cell, pos, p)
    r = dp.compute(o['cell'], o['positions'])
    vol = abs(np.linalg.det(o['cell']))
    return {'mode': mode, 'natoms': n, 'target_pressure_GPa': p, 'energy_eV': r6(r['energy']), 'volume_A3': r6(vol),
            'enthalpy_eV': r6(r['energy'] + p / 160.21766208 * vol), 'stress_kB': stress_kb(r['stress_eV_A3']),
            'max_force_eV_A': r6(np.linalg.norm(r['forces'], axis=1).max()), 'iterations': o['iterations'],
            'cell_A': np.round(o['cell'], 6).tolist(), 'positions_A': np.round(o['positions'], 6).tolist()}


def run(request, workdir):
    if request['kind'] == 'record':
        return {'id': request['id'], **RECORDS[request['id']],
                'source': 'retrieved author-released VASP result (no new calculation)'}
    return {'model': 'released DeePMD potential (Fantasy-frz.pb), numpy evaluation',
            'results': [run_structure(s) for s in request['structures']]}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerances']
    fields = ['transition_pressures_GPa', 'zero_pressure_energy_differences_meV_per_atom',
              'metastable_order_at_0GPa', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    tp, de = payload['transition_pressures_GPa'], payload['zero_pressure_energy_differences_meV_per_atom']
    if not isinstance(tp, dict) or set(tp) != set(ref['transition_pressures_GPa']) or not all(finite(v) for v in tp.values()):
        return False, {'error': 'transition_pressures_GPa needs finite numbers for %s' % sorted(ref['transition_pressures_GPa'])}
    if not isinstance(de, dict) or set(de) != set(ref['zero_pressure_energy_differences_meV_per_atom']) or \
            not all(finite(v) for v in de.values()):
        return False, {'error': 'zero_pressure_energy_differences_meV_per_atom needs finite numbers for %s'
                       % sorted(ref['zero_pressure_energy_differences_meV_per_atom'])}
    order = payload['metastable_order_at_0GPa']
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for k, v in ref['transition_pressures_GPa'].items():
        errors['P_' + k] = tp[k] - v
        checks['P_' + k] = abs(tp[k] - v) <= tol['transition_pressures_GPa'][k]
    for k, v in ref['zero_pressure_energy_differences_meV_per_atom'].items():
        errors['dE_' + k] = de[k] - v
        checks['dE_' + k] = abs(de[k] - v) <= tol['zero_pressure_energy_differences_meV_per_atom'][k]
    checks['metastable_order'] = order == ref['metastable_order_at_0GPa']
    return all(checks.values()), {'checks': checks, 'errors': errors}
