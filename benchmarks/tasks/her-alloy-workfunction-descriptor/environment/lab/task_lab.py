"""Lab logic for her-alloy-workfunction-descriptor: retrieval of author-released data and grading.

Purchases return released SECCM currents (Lib-1, Lib-2 by wafer row and potential; a small
number of Lib-3 calibration spots), released DFT adsorption energies and work functions, and
released DFT hydrogen-coverage paths. The Lib-3 currents of `spots_to_predict` are never
served; they grade the submission together with reference values derived from the released data.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
MAX_LIB3 = 5
_bought_lib3 = set()

USAGE = {
    'seccm_row': '{"kind": "seccm_row", "library": "Lib-1"|"Lib-2", "row": <row index>, "potential_mV": -100|-200|-300|-400}',
    'lib3_spot': '{"kind": "lib3_spot", "mdi": <id from calibration_spots_available>}',
    'coverage_path': '{"kind": "coverage_path", "id": <id from coverage_paths>}',
    'dft_composition': '{"kind": "dft_composition", "id": <id from dft_compositions>}',
}


def validate(request):
    kind = request.get('kind')
    if kind not in USAGE:
        raise ValueError('kind must be one of %s' % sorted(USAGE))
    try:
        if kind == 'seccm_row':
            ok = (set(request) == {'kind', 'library', 'row', 'potential_mV'} and request['library'] in ('Lib-1', 'Lib-2')
                  and type(request['row']) is int and type(request['potential_mV']) is int
                  and '%s|%d|%d' % (request['library'], request['row'], -request['potential_mV']) in DATA['seccm'])
        elif kind == 'lib3_spot':
            ok = set(request) == {'kind', 'mdi'} and type(request['mdi']) is int and str(request['mdi']) in DATA['lib3']
            if ok and request['mdi'] not in _bought_lib3 and len(_bought_lib3) >= MAX_LIB3:
                raise ValueError('At most %d Lib-3 calibration spots can be purchased' % MAX_LIB3)
        elif kind == 'coverage_path':
            ok = set(request) == {'kind', 'id'} and request['id'] in DATA['coverage']
        else:
            ok = set(request) == {'kind', 'id'} and request['id'] in DATA['dft']
    except TypeError:
        ok = False
    if not ok:
        raise ValueError('Use ' + USAGE[kind])
    return PUBLIC['requests'][kind]['cost']


def run(request, workdir):
    kind = request['kind']
    src = 'retrieved author-released data (no new measurement or calculation)'
    if kind == 'seccm_row':
        key = '%s|%d|%d' % (request['library'], request['row'], -request['potential_mV'])
        return {'library': request['library'], 'row': request['row'], 'potential_mV_vs_RHE': request['potential_mV'],
                'spots': DATA['seccm'][key], 'units': 'A/cm^2 (cathodic negative)', 'source': src}
    if kind == 'lib3_spot':
        _bought_lib3.add(request['mdi'])
        return {'library': 'Lib-3', 'mdi': request['mdi'],
                'j_A_cm2_by_potential_mV': DATA['lib3'][str(request['mdi'])], 'source': src}
    if kind == 'coverage_path':
        comp, deco = request['id'].split('/')
        return {'id': request['id'], 'composition_id': comp, 'decoration': deco, **DATA['coverage'][request['id']],
                'method': 'RPBE/GPAW fcc(111) 3x3x4 slab, 9 fcc hollows; lowest-energy nested filling path from 1 to 9 H. '
                          'Energies relative to 1/2 E(H2); no ZPE or entropy. single_site_adsorption_energies_eV: one H alone '
                          'in each hollow (1/9 ML). differential_adsorption_energy_eV: E(nH) - E((n-1)H) - 1/2 E(H2).',
                'source': src}
    rec = DATA['dft'][request['id']]
    return {'id': request['id'], **rec,
            'method': 'RPBE/GPAW fcc(111) 3x3x4 slab, random decorations (variants); dilute (1/9 ML) H adsorption energy per '
                      'fcc hollow relative to 1/2 E(H2), no ZPE or entropy; work_function_eV of the bare slab (adsorption side).',
            'source': src}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t = PRIVATE['thresholds']
    targets = PRIVATE['lib3_targets_A_cm2']
    fields = ['lib3_log10_abs_current', 'volcano_centre_eV', 'coverage_shift_eV', 'alpha_wf', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    pred = payload['lib3_log10_abs_current']
    if (not isinstance(pred, dict) or set(pred) != set(targets)
            or not all(isinstance(v, list) and len(v) == 4 and all(finite(x) for x in v) for v in pred.values())):
        return False, {'error': 'lib3_log10_abs_current must map every spots_to_predict MDI (as a string) to 4 finite numbers '
                                '(-100, -200, -300, -400 mV)'}
    for f in ('volcano_centre_eV', 'coverage_shift_eV', 'alpha_wf'):
        if not finite(payload[f]):
            return False, {'error': '%s must be a finite number' % f}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    err = [p - math.log10(-j) for m, js in targets.items() for p, j in zip(pred[m], js) if j < -1e-9]
    rmse = math.sqrt(sum(e * e for e in err) / len(err))
    lo_c, hi_c = t['volcano_centre_window_eV']
    lo_s, hi_s = t['coverage_shift_window_eV']
    lo_w, hi_w = t['alpha_wf_window']
    checks = {'evidence_executed': evidence,
              'lib3_predictions': rmse <= t['lib3_rmse_log10'],
              'volcano_centre': lo_c <= payload['volcano_centre_eV'] <= hi_c,
              'coverage_shift': lo_s <= payload['coverage_shift_eV'] <= hi_s,
              'alpha_wf': lo_w <= payload['alpha_wf'] <= hi_w}
    report = {'checks': checks, 'lib3_rmse_log10': rmse, 'n_graded_cells': len(err),
              'volcano_centre_eV': payload['volcano_centre_eV'], 'coverage_shift_eV': payload['coverage_shift_eV'],
              'alpha_wf': payload['alpha_wf']}
    return all(checks.values()), report
