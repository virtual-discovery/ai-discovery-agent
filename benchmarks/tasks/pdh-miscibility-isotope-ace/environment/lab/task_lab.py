"""Lab logic for pdh-miscibility-isotope-ace.

The lab evaluates agent-built Pd-H structures with the authors' released atomic cluster expansion
potential PdH_PBE_v03a (Mrovec, Qamar et al. 2026, Zenodo 20591788) through a numpy port of the ACE
C-tilde evaluator (ace_numpy.py, validated against python-ace to 1e-11 eV/atom): batched single
points, relaxations (fixed cell, isotropic volume or full cell) and frozen-host energy scans of one
atom. It also serves a subset of the authors' DFT (PBE) training/test records for validation.
Reference values are never computed here; they live in private/hidden.json and grade the submission.
"""
import json
import math
from pathlib import Path
import numpy as np

from ace_numpy import ACE, scan_atom
from relax import relax

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
DFT = {r['id']: r for r in json.loads((HERE/'dft_subset.json').read_text())}
_ACE = {}
MAX_ATOMS = 300
MAX_BATCH = 8
MAX_POINTS = 10000
EV_A3_GPA = 160.21766208


def ace():
    if 'm' not in _ACE:
        _ACE['m'] = ACE(str(HERE/'PdH_ace_ctilde.json'))
    return _ACE['m']


def parse_structure(s):
    if not isinstance(s, dict) or set(s) != {'symbols', 'positions', 'cell'}:
        raise ValueError('structure must be {"symbols": [...], "positions": [[x, y, z], ...], "cell": [[a], [b], [c]]} '
                         '(Cartesian Angstrom, periodic in all three directions)')
    sym = s['symbols']
    try:
        pos = np.array(s['positions'], float)
        cell = np.array(s['cell'], float)
    except (TypeError, ValueError):
        raise ValueError('positions and cell must be numeric arrays') from None
    if not isinstance(sym, list) or not sym or not all(x in ('Pd', 'H') for x in sym):
        raise ValueError('symbols must be a non-empty list of "Pd" and "H"')
    if pos.shape != (len(sym), 3) or cell.shape != (3, 3) or not np.all(np.isfinite(pos)) or not np.all(np.isfinite(cell)):
        raise ValueError('positions must be N x 3 and cell 3 x 3, finite')
    if len(sym) > MAX_ATOMS:
        raise ValueError('at most %d atoms per structure' % MAX_ATOMS)
    vol = abs(np.linalg.det(cell))
    if vol < 1.0:
        raise ValueError('cell is singular')
    heights = [vol / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3])) for k in range(3)]
    if min(heights) < 2.0:
        raise ValueError('every cell height must be at least 2 Angstrom')
    if len(sym) / vol > 0.5:
        raise ValueError('density above 0.5 atoms per cubic Angstrom; check units (Angstrom)')
    pos_w, I, J, V = ace().neighbours(cell, pos)
    if len(V) and np.linalg.norm(V, axis=1).min() < 0.5:
        raise ValueError('atoms closer than 0.5 Angstrom; check the structure')
    return sym, pos, cell


def validate(request):
    kind = request.get('kind')
    allowed = {'energy': {'kind', 'structures'},
               'relax': {'kind', 'structure', 'mode', 'fixed', 'fmax'},
               'scan': {'kind', 'structure', 'index', 'points', 'grid'},
               'dft_reference': {'kind', 'ids'}}
    if kind not in allowed:
        raise ValueError('kind must be one of %s' % sorted(allowed))
    if not set(request) <= allowed[kind]:
        raise ValueError('allowed keys for %s: %s' % (kind, sorted(allowed[kind])))
    if kind == 'energy':
        st = request.get('structures')
        if not isinstance(st, list) or not 1 <= len(st) <= MAX_BATCH:
            raise ValueError('structures must be a list of 1-%d structures' % MAX_BATCH)
        for s in st:
            parse_structure(s)
        return 1
    if kind == 'relax':
        parse_structure(request.get('structure'))
        if request.get('mode', 'positions') not in ('positions', 'volume', 'cell'):
            raise ValueError('mode must be "positions", "volume" or "cell"')
        fixed = request.get('fixed', [])
        n = len(request['structure']['symbols'])
        if not isinstance(fixed, list) or not all(type(i) is int and 0 <= i < n for i in fixed):
            raise ValueError('fixed must be a list of atom indices')
        fmax = request.get('fmax', 0.001)
        if not isinstance(fmax, (int, float)) or not 1e-4 <= fmax <= 0.1:
            raise ValueError('fmax must be between 1e-4 and 0.1 eV/A')
        return 1
    if kind == 'scan':
        sym, pos, cell = parse_structure(request.get('structure'))
        idx = request.get('index')
        if type(idx) is not int or not 0 <= idx < len(sym):
            raise ValueError('index must be the integer index of the atom to move')
        if ('points' in request) == ('grid' in request):
            raise ValueError('give exactly one of "points" (list of [x, y, z]) or "grid" ({"center": [x, y, z], "half_width": w, "spacing": h})')
        scan_points(request)
        return 2
    ids = request.get('ids')
    if not isinstance(ids, list) or not 1 <= len(ids) <= 4 or not all(i in DFT for i in ids):
        raise ValueError('ids must be a list of 1-4 ids from dft_reference_records in the catalog')
    return 1


def scan_points(request):
    """Validated scan points: an explicit list, or a cubic grid centre + k*spacing for |k*spacing| <= half_width."""
    if 'grid' in request:
        g = request['grid']
        if not isinstance(g, dict) or set(g) != {'center', 'half_width', 'spacing'}:
            raise ValueError('grid must be {"center": [x, y, z], "half_width": w, "spacing": h} (Angstrom)')
        try:
            c = np.array(g['center'], float)
            w, h = float(g['half_width']), float(g['spacing'])
        except (TypeError, ValueError):
            raise ValueError('grid values must be numeric') from None
        if c.shape != (3,) or not np.all(np.isfinite(c)) or not (0 < h <= w <= 2.0) or h < 0.02:
            raise ValueError('grid needs a 3-vector center, 0.02 <= spacing <= half_width <= 2 A')
        n = int(math.floor(w / h + 1e-9))
        if (2 * n + 1) ** 3 > MAX_POINTS:
            raise ValueError('grid has %d points; at most %d' % ((2 * n + 1) ** 3, MAX_POINTS))
        k = np.arange(-n, n + 1) * h
        X, Y, Z = np.meshgrid(k, k, k, indexing='ij')
        return c + np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    try:
        pts = np.array(request['points'], float)
    except (TypeError, ValueError):
        raise ValueError('points must be a list of [x, y, z] Cartesian positions (Angstrom)') from None
    if pts.ndim != 2 or pts.shape[1] != 3 or not 1 <= len(pts) <= MAX_POINTS or not np.all(np.isfinite(pts)):
        raise ValueError('points must be a list of 1-%d [x, y, z] positions' % MAX_POINTS)
    if np.ptp(pts, axis=0).max() > 4.0:
        raise ValueError('scan points must lie within a 4 Angstrom box')
    return pts


def _props(out, cell):
    return {'energy_eV': float(out['energy']), 'forces_eV_A': np.round(out['forces'], 6).tolist(),
            'stress_eV_A3': np.round(out['stress'], 8).tolist(),
            'pressure_GPa': float(-np.trace(out['stress']) / 3 * EV_A3_GPA),
            'volume_A3': float(abs(np.linalg.det(cell)))}


def run(request, workdir):
    kind = request['kind']
    if kind == 'energy':
        res = []
        for s in request['structures']:
            sym, pos, cell = parse_structure(s)
            res.append(_props(ace().compute(cell, sym, pos), cell))
        return {'results': res, 'potential': 'ACE PdH_PBE_v03a'}
    if kind == 'relax':
        sym, pos, cell = parse_structure(request['structure'])
        r = relax(ace(), cell, sym, pos, mode=request.get('mode', 'positions'), fixed=request.get('fixed', []),
                  fmax=request.get('fmax', 0.001), max_steps=1500)
        out = {'energy_eV': float(r['energy']), 'cell': np.round(r['cell'], 8).tolist(),
               'positions': np.round(r['positions'], 8).tolist(), 'symbols': sym,
               'forces_eV_A': np.round(r['forces'], 6).tolist(), 'stress_eV_A3': np.round(r['stress'], 8).tolist(),
               'pressure_GPa': float(-np.trace(r['stress']) / 3 * EV_A3_GPA),
               'volume_A3': float(abs(np.linalg.det(r['cell']))), 'max_force_eV_A': r['max_force'],
               'converged': r['converged'], 'optimizer_steps': r['steps'], 'mode': request.get('mode', 'positions')}
        return out
    if kind == 'scan':
        sym, pos, cell = parse_structure(request['structure'])
        pts = scan_points(request)
        e = scan_atom(ace(), cell, sym, pos, request['index'], pts)
        out = {'energies_eV': np.round(e, 8).tolist(), 'index': request['index'],
               'note': 'total energies with the moved atom at each point and every other atom frozen'}
        if 'grid' in request:
            out['grid_order'] = 'C order over (x, y, z) offsets -n..n times spacing (x slowest), n = floor(half_width/spacing)'
            out['n_per_axis'] = int(round(len(pts) ** (1 / 3)))
        return out
    return {'records': [{k: DFT[i][k] for k in ('id', 'label', 'energy_eV', 'forces_eV_A', 'stress_eV_A3')}
                        for i in request['ids']],
            'source': 'retrieved author-released DFT (VASP PBE) record; no new calculation'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerance']
    fields = list(ref) + ['evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bad = [f for f in ref if not finite(payload[f])]
    if bad:
        return False, {'error': 'fields must be finite numbers: %s' % bad}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    report = {}
    for f, r in ref.items():
        v = payload[f]
        kind, t = tol[f]
        if kind == 'abs':
            ok = abs(v - r) <= t
            report[f] = {'value': v, 'error': v - r}
        elif kind == 'rel':
            ok = abs(v / r - 1) <= t
            report[f] = {'value': v, 'relative_error': v / r - 1}
        else:   # 'log': factor tolerance on a positive quantity
            ok = v > 0 and abs(math.log(v / r)) <= math.log(t)
            report[f] = {'value': v, 'ratio_to_reference': v / r if r else None}
        checks[f] = bool(ok)
    report['checks'] = checks
    return all(checks.values()), report
