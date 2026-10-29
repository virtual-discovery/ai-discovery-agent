"""Lab logic for co3sn2s2-dmi-diffuse: retrieval of author-archived DMC diffuse-scattering arrays and grading.

Purchases return the archived DMCpy Cut2D channels (counts, monitor, normalization, normcount) with their bin
edges: whole narrow scans around (0 0 3), or square regions of the combined wide scans around a listed (h 0 l) position.
The narrow scans at 170, 175 and 200 K and the (1 0 1) region of the 180 K wide scan are never served; annular
profiles reduced from them grade the submission, together with two properties of the archived spin Hamiltonian.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
SOURCE = 'retrieved author-archived DMC data (no new measurement)'
WIDE_T = {1.5: '1.5', 165: '165', 180: '180', 300: '300'}
NOTE = ('Arrays are indexed [ix][iy]; qx_edges (len nx+1) run along a*, qy_edges (len ny+1) along c*, in 1/A. '
        'normalized intensity = counts*normcount/(monitor*normalization).')


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def resolve(request):
    kind = request.get('kind')
    if kind == 'narrow_scan':
        T = request.get('temperature_K')
        if set(request) != {'kind', 'temperature_K'} or not number(T):
            raise ValueError('Use {"kind": "narrow_scan", "temperature_K": <from narrow_scans_available>}')
        key = '%g' % T
        if key in ('170', '175', '200'):
            raise ValueError('The %s K narrow scan is withheld; it grades the submission' % key)
        if key not in DATA['narrow']:
            raise ValueError('temperature_K must be one of ' + ', '.join(sorted(DATA['narrow'], key=float)))
        return kind, key
    if kind == 'wide_region':
        T, hkl = request.get('temperature_K'), request.get('hkl')
        if set(request) != {'kind', 'temperature_K', 'hkl'} or not number(T) or not (
                isinstance(hkl, list) and len(hkl) == 3 and all(number(v) for v in hkl) and hkl[1] == 0):
            raise ValueError('Use {"kind": "wide_region", "temperature_K": 1.5|165|180|300, "hkl": [h, 0, l]}')
        if '%g' % T not in WIDE_T.values():
            raise ValueError('temperature_K must be 1.5, 165, 180 or 300')
        key = '%g|%g,%g' % (T, hkl[0], hkl[2])
        if key == '180|1,1':
            raise ValueError('The (1 0 1) region of the 180 K wide scan is withheld; it grades the submission')
        if key not in DATA['wide']:
            raise ValueError('hkl must be one of the wide_regions_available positions in /app/catalog.json')
        return kind, key
    raise ValueError('kind must be "narrow_scan" or "wide_region"')


def validate(request):
    kind, _ = resolve(request)
    return COSTS[kind]


def run(request, workdir):
    kind, key = resolve(request)
    if kind == 'narrow_scan':
        return {'temperature_K': float(key), **DATA['narrow'][key], 'note': NOTE, 'source': SOURCE}
    T, hl = key.split('|')
    return {'temperature_K': float(T), **DATA['wide'][key], 'note': NOTE, 'source': SOURCE}


def finite_list(x, n):
    return isinstance(x, list) and len(x) == n and all(number(v) for v in x)


def chi2r(pred, truth):
    return sum(((p - t)/e)**2 for p, t, e in zip(pred, truth['values'], truth['errors']))/len(pred)


def grade(payload, trace, validation):
    th = PRIVATE['thresholds']
    fields = ['narrow_003_profiles', 'wide_101_profile_180K', 'canting_angle_deg', 'spin_wave_gap_meV', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    narrow = payload['narrow_003_profiles']
    targets = PRIVATE['narrow_profiles']
    if not (isinstance(narrow, dict) and set(narrow) == set(targets) and all(finite_list(narrow[k], 8) for k in targets)):
        return False, {'error': 'narrow_003_profiles must map "170", "175" and "200" to lists of 8 finite numbers'}
    if not finite_list(payload['wide_101_profile_180K'], 8):
        return False, {'error': 'wide_101_profile_180K must be a list of 8 finite numbers'}
    if not (number(payload['canting_angle_deg']) and number(payload['spin_wave_gap_meV'])):
        return False, {'error': 'canting_angle_deg and spin_wave_gap_meV must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    scores = {'narrow_%s_chi2r' % k: chi2r(narrow[k], targets[k]) for k in sorted(targets)}
    scores['wide_101_chi2r'] = chi2r(payload['wide_101_profile_180K'], PRIVATE['wide_101_profile_180K'])
    tilt_err = payload['canting_angle_deg'] - PRIVATE['canting_angle_deg']
    gap_rel = payload['spin_wave_gap_meV']/PRIVATE['gap_meV'] - 1
    checks = {'evidence_executed': evidence}
    checks.update({'narrow_profile_%s' % k: scores['narrow_%s_chi2r' % k] <= th['narrow_profile_chi2r'] for k in sorted(targets)})
    checks['wide_101_profile'] = scores['wide_101_chi2r'] <= th['wide_101_profile_chi2r']
    checks['canting_angle'] = abs(tilt_err) <= th['canting_angle_tol_deg']
    checks['spin_wave_gap'] = abs(gap_rel) <= th['gap_rel_tol']
    return all(checks.values()), {'checks': checks, 'scores': scores, 'canting_angle_error_deg': tilt_err,
                                  'gap_relative_error': gap_rel}
