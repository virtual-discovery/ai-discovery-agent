"""Lab logic for gan-irradiation-defect-levels: retrieval of author-released VASP outputs and grading.

A purchase returns one released calculation of the GaN defect archive (Hauck, Jin & Tuttle, Zenodo 21893631):
its run settings as printed in the OUTCAR, k-point, electron count, energies, Fermi energy, initial and final
structures, final Kohn-Sham eigenvalues and occupations, per-atom core potentials and spin moments. A separate
purchase returns the builders' DFPT dielectric tensors. Transition levels, the band gap, the divacancy binding
energy and the Frenkel-pair energies are never served; reference values derived from the archive grade the
submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-released VASP output (no new calculation)'
LEVEL_DEFECTS = ('V_Ga', 'V_N', 'Ga_i', 'N_i', 'C_N')


def validate(request):
    kind = request.get('kind')
    if kind == 'calculation':
        if set(request) != {'kind', 'id'} or request['id'] not in DATA['calculations']:
            raise ValueError('Use {"kind": "calculation", "id": <id from calculations in /app/catalog.json>}')
        return 1
    if kind == 'dielectric':
        if set(request) != {'kind'}:
            raise ValueError('Use {"kind": "dielectric"}')
        return 1
    raise ValueError('kind must be "calculation" or "dielectric"')


def run(request, workdir):
    if request['kind'] == 'dielectric':
        return {**DATA['dielectric'], 'source': 'computed by the task builders (not part of the archive)'}
    return {'id': request['id'], **DATA['calculations'][request['id']], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _merge_narrow(levels, width):
    """Treat a charge state that is stable over less than `width` eV as unstable: replace the two levels around it
    by the single crossing of its neighbours. Levels in decreasing order are left untouched (not a stable ladder)."""
    out = sorted(levels, key=lambda s: s['level_eV'])
    k = 0
    while k + 1 < len(out):
        a, b = out[k], out[k + 1]
        if a['to'] == b['from'] and 0 <= b['level_eV'] - a['level_eV'] < width and a['from'] != b['to']:
            da, db = a['from'] - a['to'], b['from'] - b['to']
            out[k:k + 2] = [{'from': a['from'], 'to': b['to'],
                             'level_eV': (da * a['level_eV'] + db * b['level_eV']) / (da + db)}]
            k = max(k - 1, 0)
        else:
            k += 1
    return out


def _levels_ok(submitted, reference, gap, t):
    """Core-window matching: every reference level inside [w, gap - w] must be reported with the same pair of
    charge states within tol; a reported level inside [w + tol, gap - w - tol] that matches no reference core
    level is spurious."""
    w, tol = t['level_window_eV'], t['level_tol_eV']
    order = [s['level_eV'] for s in submitted]
    inverted = any(y < x for x, y in zip(order, order[1:]))
    submitted = submitted if inverted else _merge_narrow(submitted, t['merge_width_eV'])
    reference = _merge_narrow(reference, t['merge_width_eV'])
    core = [r for r in reference if w <= r['level_eV'] <= gap - w]
    used, errors = set(), []
    for r in core:
        hits = [k for k, s in enumerate(submitted) if k not in used and s['from'] == r['from'] and s['to'] == r['to']
                and abs(s['level_eV'] - r['level_eV']) <= tol]
        if not hits:
            errors.append('missing (%+d/%+d) near %.2f eV' % (r['from'], r['to'], r['level_eV']))
            continue
        used.add(min(hits, key=lambda k: abs(submitted[k]['level_eV'] - r['level_eV'])))
    if inverted:
        errors.append('levels not in increasing order')
    for k, s in enumerate(submitted):
        if k not in used and w + tol <= s['level_eV'] <= gap - w - tol:
            errors.append('unexpected (%+d/%+d) at %.2f eV' % (s['from'], s['to'], s['level_eV']))
    deviations = {}
    for r in core:
        best = [s['level_eV'] - r['level_eV'] for s in submitted if s['from'] == r['from'] and s['to'] == r['to']]
        if best:
            deviations['%+d/%+d' % (r['from'], r['to'])] = min(best, key=abs)
    return not errors, errors, deviations


def grade(payload, trace, validation):
    t, ref = PRIVATE['thresholds'], PRIVATE['reference']
    fields = ['band_gap_eV', 'transition_levels', 'divacancy_binding_energy_eV', 'frenkel_pair_energy_eV',
              'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    levels = payload['transition_levels']
    level_entry = lambda s: (isinstance(s, dict) and set(s) == {'from', 'to', 'level_eV'} and type(s['from']) is int
                             and type(s['to']) is int and finite(s['level_eV']))
    if not isinstance(levels, dict) or set(levels) != set(LEVEL_DEFECTS) or not all(
            isinstance(v, list) and all(level_entry(s) for s in v) for v in levels.values()):
        return False, {'error': 'transition_levels must map each of %s to a list of {"from", "to", "level_eV"} '
                                '(integer charges, finite level)' % (LEVEL_DEFECTS,)}
    fp = payload['frenkel_pair_energy_eV']
    if not isinstance(fp, dict) or set(fp) != {'N', 'Ga'} or not all(finite(v) for v in fp.values()):
        return False, {'error': 'frenkel_pair_energy_eV must give finite "N" and "Ga" values'}
    if not finite(payload['band_gap_eV']) or not finite(payload['divacancy_binding_energy_eV']):
        return False, {'error': 'band_gap_eV and divacancy_binding_energy_eV must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    report = {}
    gap = ref['band_gap_eV']
    d_gap = payload['band_gap_eV'] - gap
    checks['band_gap'] = abs(d_gap) <= t['gap_tol_eV']
    report['band_gap_error_eV'] = d_gap
    for d in LEVEL_DEFECTS:
        ok, errors, dev = _levels_ok(levels[d], ref['levels'][d], gap, t)
        checks['levels_' + d] = ok
        report['levels_' + d] = {'errors': errors, 'deviations_eV': dev}
    d_b = payload['divacancy_binding_energy_eV'] - ref['divacancy_binding_energy_eV']
    checks['divacancy_binding_energy'] = abs(d_b) <= t['binding_tol_eV']
    report['divacancy_binding_error_eV'] = d_b
    for k in ('N', 'Ga'):
        e = fp[k] - ref['frenkel_pair_energy_eV'][k]
        checks['frenkel_pair_energy_' + k] = abs(e) <= t['frenkel_tol_eV'][k]
        report['frenkel_pair_error_eV_' + k] = e
    report['checks'] = checks
    return all(checks.values()), report
