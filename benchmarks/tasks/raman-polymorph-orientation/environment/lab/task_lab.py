"""Lab logic for raman-polymorph-orientation: simulated polarized micro-Raman spectra of a hidden SnF2 crystal.

Each request returns one spectrometer window of a backscattering spectrum from the polished face (face A).
The spectra are synthesized from a forward model of the hidden sample built at construction time from the
Computational Raman Database records (private/lab_model.json): phonon frequencies and lab-frame Raman tensors
of the main phase (including LO mixing for the phonon wavevector along the face normal), orientation-averaged
tensors of a minority phase, thermal and scattered-photon factors, instrument resolution, a calibration offset
of the Raman-shift axis, background, laser-power drift and counting noise. Hidden targets are never served.
"""
import json
import math
import zlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
MODEL = json.loads((HERE/'private/lab_model.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
HC_K = 1.4387769
KEYS = {'kind', 'polarizer_deg', 'analyzer_deg', 'center_cm1', 'accumulations'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def validate(request):
    usage = ('Use {"kind": "raman", "polarizer_deg": <angle>, "analyzer_deg": <angle or null>, '
             '"center_cm1": <window centre>, "accumulations": <1-4>}')
    if request.get('kind') != 'raman' or not set(request) <= KEYS or not {'kind', 'polarizer_deg', 'analyzer_deg', 'center_cm1'} <= set(request):
        raise ValueError(usage)
    if not finite(request['polarizer_deg']) or not (request['analyzer_deg'] is None or finite(request['analyzer_deg'])):
        raise ValueError('polarizer_deg must be a number and analyzer_deg a number or null')
    lo, hi = MODEL['center_range_cm1']
    if not finite(request['center_cm1']) or not lo <= request['center_cm1'] <= hi:
        raise ValueError('center_cm1 must lie in [%g, %g]' % (lo, hi))
    acc = request.get('accumulations', 1)
    if type(acc) is not int or not 1 <= acc <= 4:
        raise ValueError('accumulations must be an integer from 1 to 4')
    return acc


def pol(angle_deg):
    a = math.radians(angle_deg)
    return np.array([math.cos(a), math.sin(a), 0.0])


def pseudo_voigt(x, fwhm_l, fwhm_g):
    f = (fwhm_g ** 5 + 2.69269 * fwhm_g ** 4 * fwhm_l + 2.42843 * fwhm_g ** 3 * fwhm_l ** 2
         + 4.47163 * fwhm_g ** 2 * fwhm_l ** 3 + 0.07842 * fwhm_g * fwhm_l ** 4 + fwhm_l ** 5) ** 0.2
    r = fwhm_l / f
    eta = 1.36603 * r - 0.47719 * r ** 2 + 0.11116 * r ** 3
    g, s = f / 2, f / 2.354820045
    return eta * (g / math.pi) / (x * x + g * g) + (1 - eta) * np.exp(-x * x / (2 * s * s)) / (s * math.sqrt(2 * math.pi))


def line_strengths(model, polarizer, analyzer):
    ei = pol(polarizer)
    if analyzer is None:
        outs = [pol(0.0), pol(90.0)]
    else:
        outs = [pol(analyzer)]
    main = [sum(float(es @ np.array(line['R']) @ ei) ** 2 for es in outs) for line in model['main_lines']]
    minority = []
    for line in model['minority_lines']:
        if analyzer is None:
            minority.append(line['ipar'] + line['icross'])
        else:
            d = math.radians(analyzer - polarizer)
            minority.append(line['ipar'] * math.cos(d) ** 2 + line['icross'] * math.sin(d) ** 2)
    return main, minority


def expected_counts(model, request):
    """Noise-free spectrum (counts above the detector offset) on the window grid."""
    acc = request.get('accumulations', 1)
    half, step = model['half_window_cm1'], model['pixel_cm1']
    n = int(round(2 * half / step)) + 1
    nu = request['center_cm1'] - half + step * np.arange(n)
    nu_true = nu - model['axis_offset_cm1']          # spectrometer axis -> true Raman shift
    laser = model['laser_true_cm1']
    T = model['T_K']
    main, minority = line_strengths(model, request['polarizer_deg'], request['analyzer_deg'])
    lines = [(l['f'], a, l['w']) for l, a in zip(model['main_lines'], main)]
    lines += [(l['f'], a, l['w']) for l, a in zip(model['minority_lines'], minority)]
    signal = np.zeros(n)
    for f, a, w in lines:
        if a <= 0:
            continue
        nb = 1.0 / math.expm1(HC_K * f / T)
        stokes = ((laser - f) / laser) ** 3 * (nb + 1) / f      # photon counting: omega_s^3 (n+1)/omega
        anti = ((laser + f) / laser) ** 3 * nb / f
        signal += a * (stokes * pseudo_voigt(nu_true - f, w, model['resolution_fwhm_cm1'])
                       + anti * pseudo_voigt(nu_true + f, w, model['resolution_fwhm_cm1']))
    b0, b1 = model['background']
    mean = acc * (model['intensity_scale'] * signal + b0 + b1 * nu_true)
    keep = np.abs(nu_true) >= model['notch_half_width_cm1']
    return nu, mean, keep


def synthesize(model, request, seed):
    rng = np.random.default_rng(seed)
    acc = request.get('accumulations', 1)
    nu, mean, keep = expected_counts(model, request)
    drift = 1.0 + model['power_drift_sd'] * rng.standard_normal()
    counts = rng.poisson(np.maximum(mean * drift, 0)) + acc * model['detector_offset_counts']
    counts = counts + rng.normal(0.0, model['read_noise_counts'] * math.sqrt(acc), len(nu))
    return nu[keep], np.round(counts[keep]).astype(int)


def request_seed(request, workdir):
    try:
        job = int(Path(workdir).name.split('-')[-1])
    except ValueError:
        job = 0
    text = json.dumps(request, sort_keys=True)
    return [MODEL['seed'], job, zlib.crc32(text.encode())]


def run(request, workdir):
    nu, counts = synthesize(MODEL, request, request_seed(request, workdir))
    return {'polarizer_deg': request['polarizer_deg'], 'analyzer_deg': request['analyzer_deg'],
            'center_cm1': request['center_cm1'], 'accumulations': request.get('accumulations', 1),
            'raman_shift_cm1': [round(float(x), 2) for x in nu], 'counts': [int(c) for c in counts],
            'note': 'face A backscattering; pixels inside the notch filter are omitted; counts include the detector offset'}


def axis_error_deg(c, c_true):
    c = np.asarray(c, float) / np.linalg.norm(c)
    ct = np.asarray(c_true, float)
    return min(math.degrees(math.acos(min(1.0, abs(float(c @ x))))) for x in (ct, ct * np.array([1.0, 1.0, -1.0])))


def norm_id(x):
    return x.strip().lower() if isinstance(x, str) else None


def grade(payload, trace, validation):
    t, ref = PRIVATE['thresholds'], PRIVATE['reference']
    fields = ['main_phase', 'minority_phase', 'frequency_scale', 'c_axis_lab',
              'faceB_strongest_mode_ratio_tt_over_zz', 'faceB_highest_optical_phonon_cm1', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    axis = payload['c_axis_lab']
    if not (isinstance(axis, list) and len(axis) == 3 and all(finite(v) for v in axis) and sum(v * v for v in axis) > 1e-12):
        return False, {'error': 'c_axis_lab must be a list of three finite numbers (not all zero)'}
    numbers = ['frequency_scale', 'faceB_strongest_mode_ratio_tt_over_zz', 'faceB_highest_optical_phonon_cm1']
    if not all(finite(payload[k]) for k in numbers) or payload['faceB_strongest_mode_ratio_tt_over_zz'] <= 0:
        return False, {'error': 'frequency_scale, faceB_* must be finite numbers (ratio > 0)'}
    if norm_id(payload['main_phase']) is None or norm_id(payload['minority_phase']) is None:
        return False, {'error': 'main_phase and minority_phase must be strings'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    s_err = payload['frequency_scale'] / ref['frequency_scale'] - 1
    ax_err = axis_error_deg(axis, ref['c_axis_lab'])
    r_err = math.log(payload['faceB_strongest_mode_ratio_tt_over_zz'] / ref['faceB_ratio_tt_over_zz'])
    lo_err = payload['faceB_highest_optical_phonon_cm1'] / ref['faceB_highest_optical_cm1'] - 1
    checks = {'evidence_executed': evidence,
              'main_phase': norm_id(payload['main_phase']) == ref['main_phase'],
              'minority_phase': norm_id(payload['minority_phase']) == ref['minority_phase'],
              'frequency_scale': abs(s_err) <= t['scale_rel'],
              'c_axis': ax_err <= t['axis_deg'],
              'faceB_ratio': abs(r_err) <= t['ratio_log'],
              'faceB_highest_optical_phonon': abs(lo_err) <= t['lo_rel']}
    report = {'checks': checks, 'frequency_scale_relative_error': s_err, 'c_axis_error_deg': ax_err,
              'faceB_ratio_log_error': r_err, 'faceB_highest_optical_relative_error': lo_err}
    return all(checks.values()), report
