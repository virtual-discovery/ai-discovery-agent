"""Lab logic for hzo-vacancy-local-environment: retrieval of author-released DFT results and grading.

Purchases return released total energies (VASP-PBE vacancy and pristine supercells), CASTEP vacancy enthalpies
for seven functionals with the authors' site labels, and CASTEP-PBE 3x3x3 vacancy enthalpies. The sites to
predict are never served; their archived formation energies, the random-alloy vacancy populations computed
from all archived sites, and the all-site LDA description grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-released result (no new calculation)'
COMPOSITIONS = ('x0.25', 'x0.75')


def _site(request):
    site = request.get('site')
    if type(site) is not int:
        raise ValueError('"site" must be an integer atom index')
    return str(site)


def validate(request):
    kind = request.get('kind')
    if kind == 'vacancy':
        if set(request) != {'kind', 'composition', 'site'} or request['composition'] not in DATA['vacancy']:
            raise ValueError('Use {"kind": "vacancy", "composition": "x0.25"|"x0.50"|"x0.75", "site": <O atom index>}')
        if _site(request) not in DATA['vacancy'][request['composition']]:
            raise ValueError('That site is not purchasable for this composition (see vacancy_sites_available)')
        return 1
    if kind == 'pristine':
        if set(request) != {'kind', 'composition'} or request['composition'] not in DATA['pristine']:
            raise ValueError('Use {"kind": "pristine", "composition": "x0.25"|"x0.50"|"x0.75"}')
        return 1
    if kind == 'castep_vacancy':
        if set(request) != {'kind', 'functional', 'site'} or request['functional'] not in DATA['castep_vacancy']:
            raise ValueError('Use {"kind": "castep_vacancy", "functional": one of %s, "site": <O index>}' % sorted(DATA['castep_vacancy']))
        if _site(request) not in DATA['castep_vacancy'][request['functional']]:
            raise ValueError('No released CASTEP calculation for that site')
        return 1
    if kind == 'castep_333_vacancy':
        if set(request) != {'kind', 'site'} or _site(request) not in DATA['castep_333_vacancy']:
            raise ValueError('Use {"kind": "castep_333_vacancy", "site": <index from castep_333.sites>}')
        return 2
    raise ValueError('kind must be "vacancy", "pristine", "castep_vacancy" or "castep_333_vacancy"')


def run(request, workdir):
    kind = request['kind']
    if kind == 'vacancy':
        comp, site = request['composition'], request['site']
        return {'composition': comp, 'site': site, 'total_energy_eV': DATA['vacancy'][comp][str(site)],
                'calculation': 'VASP-PBE, 96-atom supercell with one neutral O vacancy at this site, ions relaxed at the CIF lattice',
                'source': SOURCE}
    if kind == 'pristine':
        comp = request['composition']
        return {'composition': comp, 'total_energy_eV': DATA['pristine'][comp],
                'calculation': 'VASP-PBE, pristine 96-atom supercell at the CIF lattice', 'source': SOURCE}
    if kind == 'castep_vacancy':
        row = DATA['castep_vacancy'][request['functional']][str(request['site'])]
        return {'functional': request['functional'], 'site': request['site'], 'final_enthalpy_eV': row['final_enthalpy_eV'],
                'vacancy_type': row['vacancy_type'], 'local_environment': row['local_environment'],
                'n_Hf': row['n_Hf'], 'n_Zr': row['n_Zr'], 'frac_coord': row['frac_coord'],
                'calculation': 'CASTEP, 96-atom 2x2x2 supercell with one neutral O vacancy (labels as released by the authors)',
                'source': SOURCE}
    return {'site': request['site'], 'final_enthalpy_eV': DATA['castep_333_vacancy'][str(request['site'])],
            'calculation': 'CASTEP-PBE, 324-atom 3x3x3 supercell with one neutral O vacancy', 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t, ref, targets = PRIVATE['thresholds'], PRIVATE['reference'], PRIVATE['site_targets']
    fields = ['predicted_formation_energies_eV', 'vacancy_Zr_neighbour_distribution', 'vacancy_site_fraction',
              'lda_common_slope_model', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    pred, dist, frac, lda = (payload[f] for f in fields[:4])
    if not isinstance(pred, dict) or set(pred) != set(targets) or not all(finite(v) for v in pred.values()):
        return False, {'error': 'predicted_formation_energies_eV must give a finite energy for every sites_to_predict id'}
    if not isinstance(dist, dict) or set(dist) != set(COMPOSITIONS) or not all(
            isinstance(v, list) and len(v) == 5 and all(finite(p) and p >= 0 for p in v) and sum(v) > 0 for v in dist.values()):
        return False, {'error': 'vacancy_Zr_neighbour_distribution must map x0.25 and x0.75 to five non-negative fractions'}
    if not isinstance(frac, dict) or set(frac) != set(COMPOSITIONS) or not all(finite(v) and v > 0 for v in frac.values()):
        return False, {'error': 'vacancy_site_fraction must map x0.25 and x0.75 to positive numbers'}
    if not isinstance(lda, dict) or set(lda) != {'offset_O3_eV', 'slope_per_Zr_eV'} or not all(finite(v) for v in lda.values()):
        return False, {'error': 'lda_common_slope_model must give finite offset_O3_eV and slope_per_Zr_eV'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed}
    report = {}
    for x in COMPOSITIONS:
        err = [pred[k] - v for k, v in targets.items() if k.startswith(x + '-')]
        rmse = math.sqrt(sum(e * e for e in err) / len(err))
        bias = sum(err) / len(err)
        d = dist[x]
        total = sum(d)
        mean_nzr = sum(k * p for k, p in enumerate(d)) / total
        share4 = d[4] / total
        dlog = math.log10(frac[x] / ref[x]['site_fraction'])
        checks[f'predictions_{x}'] = rmse <= t['pred_rmse_max_eV'] and abs(bias) <= t['pred_bias_max_eV']
        checks[f'vacancy_distribution_{x}'] = (abs(total - 1) <= t['distribution_sum_tol']
                                               and abs(mean_nzr - ref[x]['mean_nZr']) <= t['mean_nZr_tol'][x]
                                               and abs(share4 - ref[x]['share_4Zr']) <= t['share_4Zr_tol'][x])
        checks[f'vacancy_site_fraction_{x}'] = abs(dlog) <= t['log10_site_fraction_tol'][x]
        report[x] = {'prediction_rmse_eV': rmse, 'prediction_mean_error_eV': bias, 'distribution_sum': total,
                     'mean_nZr_error': mean_nzr - ref[x]['mean_nZr'], 'share_4Zr_error': share4 - ref[x]['share_4Zr'],
                     'log10_site_fraction_error': dlog}
    d_off = lda['offset_O3_eV'] - ref['lda']['offset_O3_eV']
    d_slope = lda['slope_per_Zr_eV'] - ref['lda']['slope_per_Zr_eV']
    checks['lda_offset'] = abs(d_off) <= t['lda_offset_tol_eV']
    checks['lda_slope'] = abs(d_slope) <= t['lda_slope_tol_eV']
    report.update(checks=checks, lda_offset_error_eV=d_off, lda_slope_error_eV=d_slope)
    return all(checks.values()), report
