"""Lab logic for zno-hydrogen-shallow-donor: retrieval of author-released VASP results and grading.

Purchases return released calculations of the ZnO defect archive (Aicher et al. 2025): run settings, energies,
moments, structures and core potentials ("run"), or the final Kohn-Sham eigenvalues and occupations
("eigenvalues"). Nothing derived (formation energies, levels, Fermi levels, concentrations) is served;
the values derived from the released energies by the stated conventions grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-released VASP result (no new calculation)'
LEVEL_DEFECTS = ('OH_O', 'H_O', 'V_O')
SITES = ('BC_parallel', 'BC_perpendicular', 'AB_parallel', 'AB_perpendicular')


def validate(request):
    kind = request.get('kind')
    if kind not in ('run', 'eigenvalues'):
        raise ValueError('kind must be "run" or "eigenvalues"')
    if set(request) != {'kind', 'id'} or not isinstance(request.get('id'), str):
        raise ValueError('Use {"kind": "run"|"eigenvalues", "id": <run id from catalog runs>}')
    if request['id'] not in DATA['runs']:
        raise ValueError('Unknown run id (see the "runs" list in /app/catalog.json)')
    if kind == 'eigenvalues' and request['id'] not in DATA['eigenvalues']:
        raise ValueError('No eigenvalues are served for this run (see eigenvalues_available)')
    return 1


def run(request, workdir):
    rid = request['id']
    body = DATA['runs'][rid] if request['kind'] == 'run' else DATA['eigenvalues'][rid]
    return {'id': rid, 'kind': request['kind'], **body, 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _levels(entries, lo, hi):
    if not isinstance(entries, list):
        return None
    out = []
    for e in entries:
        if not isinstance(e, dict) or set(e) != {'from', 'to', 'level_eV'} or type(e['from']) is not int \
                or type(e['to']) is not int or not finite(e['level_eV']):
            return None
        if lo <= e['level_eV'] <= hi:
            out.append(e)
    return out


def grade(payload, trace, validation):
    t, ref, cat = PRIVATE['thresholds'], PRIVATE['reference'], PRIVATE['categorical']
    fields = ['band_gap_eV', 'delta_mu_eV', 'transition_levels', 'V_O_plus_ever_stable', 'epr_centre', 'oh_o_site',
              'fermi_level_eV', 'log10_concentration_cm3', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    p = payload
    if not finite(p['band_gap_eV']) or not finite(p['fermi_level_eV']):
        return False, {'error': 'band_gap_eV and fermi_level_eV must be finite numbers'}
    dm = p['delta_mu_eV']
    if not isinstance(dm, dict) or set(dm) != {'O', 'H'} or not all(finite(v) for v in dm.values()):
        return False, {'error': 'delta_mu_eV must give finite O and H values'}
    lc = p['log10_concentration_cm3']
    if not isinstance(lc, dict) or set(lc) != set(ref['log10_concentration_cm3']) or not all(finite(v) for v in lc.values()):
        return False, {'error': 'log10_concentration_cm3 must give finite V_O, OH_O^+ and OH_O^0 values'}
    tl = p['transition_levels']
    lo, hi = t['ignore_levels_outside_eV']
    parsed = {d: _levels(tl.get(d), lo, hi) if isinstance(tl, dict) else None for d in LEVEL_DEFECTS}
    if not isinstance(tl, dict) or set(tl) != set(LEVEL_DEFECTS) or any(v is None for v in parsed.values()):
        return False, {'error': 'transition_levels must map OH_O, H_O and V_O to lists of {"from", "to", "level_eV"}'}
    if not isinstance(p['epr_centre'], dict) or set(p['epr_centre']) != {'defect', 'charge'}:
        return False, {'error': 'epr_centre must be {"defect": <name>, "charge": <int>}'}
    if p['oh_o_site'] not in SITES:
        return False, {'error': 'oh_o_site must be one of %s' % list(SITES)}
    if not isinstance(p['V_O_plus_ever_stable'], bool):
        return False, {'error': 'V_O_plus_ever_stable must be true or false'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = p['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    report = {}
    gap = ref['band_gap_eV']
    margin = t['interior_margin_eV']

    def ladder_ok(d, want_from, want_to, ok_level):
        entries = parsed[d]
        hit = [e for e in entries if e['from'] == want_from and e['to'] == want_to and ok_level(e['level_eV'])]
        others = [e for e in entries if not (e['from'] == want_from and e['to'] == want_to)
                  and margin <= e['level_eV'] <= gap - margin]
        return len(hit) == 1 and not others

    w0, w1 = t['shallow_level_window_eV']
    vo_ref = ref['transition_levels']['V_O'][0]['level_eV']
    checks['band_gap'] = abs(p['band_gap_eV'] - gap) <= t['band_gap_tol_eV']
    checks['delta_mu_O'] = abs(dm['O'] - ref['delta_mu_eV']['O']) <= t['delta_mu_tol_eV']
    checks['delta_mu_H'] = abs(dm['H'] - ref['delta_mu_eV']['H']) <= t['delta_mu_tol_eV']
    checks['level_OH_O'] = ladder_ok('OH_O', 1, 0, lambda x: w0 <= x <= w1)
    checks['level_H_O'] = ladder_ok('H_O', 1, 0, lambda x: w0 <= x <= w1)
    checks['level_V_O'] = ladder_ok('V_O', 2, 0, lambda x: abs(x - vo_ref) <= t['V_O_level_tol_eV'])
    checks['V_O_plus_stability'] = p['V_O_plus_ever_stable'] is cat['V_O_plus_ever_stable']
    checks['epr_centre'] = (p['epr_centre'].get('defect') == cat['epr_centre']['defect']
                            and p['epr_centre'].get('charge') == cat['epr_centre']['charge']
                            and type(p['epr_centre'].get('charge')) is int)
    checks['oh_o_site'] = p['oh_o_site'] == cat['oh_o_site']
    checks['fermi_level'] = abs(p['fermi_level_eV'] - ref['fermi_level_eV']) <= t['fermi_level_tol_eV']
    for k, tol in t['log10_tol'].items():
        checks['log10_' + k] = abs(lc[k] - ref['log10_concentration_cm3'][k]) <= tol
        report['log10_error_' + k] = lc[k] - ref['log10_concentration_cm3'][k]
    report.update(checks=checks, band_gap_error_eV=p['band_gap_eV'] - gap,
                  delta_mu_error_eV={k: dm[k] - ref['delta_mu_eV'][k] for k in dm},
                  fermi_level_error_eV=p['fermi_level_eV'] - ref['fermi_level_eV'])
    return all(checks.values()), report
