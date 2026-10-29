"""Generate the staged, lab, private, solution and ablation files for raman-polymorph-orientation.

Run inside the runtime image from this folder (after fetch.py has written pool.json):
    docker run --rm --cpuset-cpus 6,7,22,23 --user $(id -u):$(id -g) -e HOME=/tmp -v <repo>:/repo \
        -w /repo/benchmarks/construction/raman-polymorph-orientation ai-theorist-benchmark-runtime:20261007 \
        python3 -B make_task.py [--skip-variants]
"""
import importlib.util
import json
import math
import sys
import zlib
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
import model as M          # noqa: E402

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/raman-polymorph-orientation'
LAB = TASK/'environment/lab'

# ----------------------------------------------------------------------------- hidden sample
HIDDEN = {
    'main_phase': 'mp-7456',            # beta-SnF2, P2_12_12_1
    'minority_phase': 'mp-2843',        # alpha-SnF2, C2/c, randomly oriented inclusions
    'minority_volume_fraction': 0.03,
    'euler_zyz_deg': [48.0, 53.0, 62.0],
    'frequency_scale': 1.0237,
    'laser_true_nm': 532.06,            # spectrometer axis assumes 532.00 nm
    'T_K': 295.0,
    'width_law': [1.0, 0.0065],         # phonon FWHM (cm-1) = (a + b omega) * U(0.85, 1.2)
    'minority_width_factor': 1.3,
    'intensity_scale': 8.0e8,
    'background': [150.0, 0.1],
    'seed': 20261019,
}
THRESHOLDS = {'scale_rel': 0.002, 'axis_deg': 6.0, 'ratio_log': math.log(1.10), 'lo_rel': 0.008}
BUDGET = 20
INSTRUMENT = {
    'laser_nominal_wavelength_nm': 532.0, 'detector_offset_counts': 500, 'read_noise_counts': 5.0,
    'resolution_fwhm_cm1': 1.2, 'pixel_cm1': 0.4, 'half_window_cm1': 160.0, 'center_range_cm1': [-450.0, 650.0],
    'notch_half_width_cm1': 18.0, 'power_drift_sd': 0.015,
}
# Purchase plan of the reference solution (also replayed here for the ablation values).
PLAN = [
    {'kind': 'raman', 'polarizer_deg': 0, 'analyzer_deg': 0, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 0, 'analyzer_deg': 0, 'center_cm1': -200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 0, 'analyzer_deg': 0, 'center_cm1': 470, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 90, 'analyzer_deg': 90, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 0, 'analyzer_deg': 90, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 45, 'analyzer_deg': 45, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 135, 'analyzer_deg': 135, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 45, 'analyzer_deg': 135, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 30, 'analyzer_deg': 30, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 60, 'analyzer_deg': 150, 'center_cm1': 200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 90, 'analyzer_deg': 90, 'center_cm1': 470, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 0, 'analyzer_deg': 90, 'center_cm1': 470, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 45, 'analyzer_deg': 45, 'center_cm1': 470, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 90, 'analyzer_deg': 90, 'center_cm1': -200, 'accumulations': 1},
    {'kind': 'raman', 'polarizer_deg': 0, 'analyzer_deg': None, 'center_cm1': 360, 'accumulations': 3},
]
LAB_AXES_ONLY = [i for i, r in enumerate(PLAN) if r['analyzer_deg'] is not None
                 and r['polarizer_deg'] % 90 == 0 and r['analyzer_deg'] % 90 == 0]


def sig(x, n=6):
    return float('%.*g' % (n, x))


def nested(x, n=6):
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if isinstance(x, list):
        return [nested(v, n) for v in x]
    if isinstance(x, float):
        return sig(x, n)
    return x


def catalog_record(e, rec):
    kv = rec['key_value_pairs']
    return {
        'mpid': e['mpid'], 'database_id': None, 'spacegroup': e['spacegroup'], 'pointgroup': e['pointgroup'],
        'cell_A': nested(e['cell']), 'atomic_numbers': e['numbers'].tolist(), 'masses_amu': nested(e['masses'], 9),
        'positions_A': nested(e['positions']), 'frequencies_cm1': nested(e['freq']),
        'irrep_labels': e['labels'], 'degenerate_sets': e['bands'],
        'raman_active': e['raman_active'].tolist(), 'ir_active': e['ir_active'].tolist(),
        'raman_tensors': [None if t is None else nested(t) for t in e['raman_tensors']],
        'eigenvectors': nested(e['eigvec']),
        'born_charges': None if e['born_rows'] is None else nested(e['born_rows']),
        'epsilon_inf': None if e['eps_inf'] is None else nested(e['eps_inf']),
        'e_above_hull_eV_per_atom': kv.get('mp_e_above_hull'), 'negative_freq_Gamma': kv.get('negative_freq_Gamma'),
        'negative_freq_path': kv.get('negative_freq_path'), 'bandgap_eV': kv.get('bandgap'),
        'dimensionality': kv.get('dimensionality'),
    }


def build_lab_model(main, minority, Q, h):
    rng = np.random.default_rng(h['seed'])
    s = h['frequency_scale']
    f, R = M.modes(main, Q.T @ np.array([0.0, 0.0, 1.0]))
    keep = (f > 5) & ((R ** 2).sum((1, 2)) > 1e-10)
    f, R = f[keep], R[keep]
    a, b = h['width_law']
    norm = math.sqrt((1 - h['minority_volume_fraction']) / main['volume'])
    main_lines = [{'f': float(s * fj), 'R': (Q @ Rj @ Q.T * norm).tolist(),
                   'w': float((a + b * s * fj) * rng.uniform(0.85, 1.2))} for fj, Rj in zip(f, R)]
    fm, Rm = M.modes(minority)
    keep = (fm > 5) & ((Rm ** 2).sum((1, 2)) > 1e-10)
    minority_lines = []
    for fj, Rj in zip(fm[keep], Rm[keep]):
        ip, ic = M.powder_invariants(Rj)
        scale = h['minority_volume_fraction'] / minority['volume']
        minority_lines.append({'f': float(s * fj), 'ipar': float(ip * scale), 'icross': float(ic * scale),
                               'w': float(h['minority_width_factor'] * (a + b * s * fj) * rng.uniform(0.85, 1.2))})
    laser_nom = 1e7 / INSTRUMENT['laser_nominal_wavelength_nm']
    laser_true = 1e7 / h['laser_true_nm']
    return {'seed': h['seed'], 'T_K': h['T_K'], 'laser_true_cm1': laser_true, 'axis_offset_cm1': laser_nom - laser_true,
            'intensity_scale': h['intensity_scale'], 'background': h['background'],
            'detector_offset_counts': INSTRUMENT['detector_offset_counts'], 'read_noise_counts': INSTRUMENT['read_noise_counts'],
            'resolution_fwhm_cm1': INSTRUMENT['resolution_fwhm_cm1'], 'pixel_cm1': INSTRUMENT['pixel_cm1'],
            'half_window_cm1': INSTRUMENT['half_window_cm1'], 'center_range_cm1': INSTRUMENT['center_range_cm1'],
            'notch_half_width_cm1': INSTRUMENT['notch_half_width_cm1'], 'power_drift_sd': INSTRUMENT['power_drift_sd'],
            'main_lines': main_lines, 'minority_lines': minority_lines}


def load_task_lab():
    spec = importlib.util.spec_from_file_location('task_lab_build', LAB/'task_lab.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def replay(lab, model, plan):
    """The oracle's purchases exactly as the lab serves them (job ids 1..n in plan order)."""
    out = []
    for job, req in enumerate(plan, start=1):
        seed = [model['seed'], job, zlib.crc32(json.dumps(req, sort_keys=True).encode())]
        nu, counts = lab.synthesize(model, req, seed)
        out.append({'polarizer_deg': req['polarizer_deg'], 'analyzer_deg': req['analyzer_deg'],
                    'center_cm1': req['center_cm1'], 'accumulations': req['accumulations'],
                    'raman_shift_cm1': [round(float(x), 2) for x in nu], 'counts': [int(c) for c in counts]})
    return out


def leaks(values, texts):
    needles = set()
    for v in values:
        if isinstance(v, str) and len(v) >= 4:
            needles.add(v)
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            for t in {f'{v:.4g}', f'{v:.5g}', f'{v:.6g}', repr(v)}:
                if len(t.replace('-', '').replace('.', '')) >= 4:
                    needles.add(t)
    return sorted({n for n in needles for t in texts if n in t})


DRIVER = r'''

# ============================================================================= driver (reference solution)
import json as _json
import os as _os
import subprocess as _subprocess
import time as _time
from pathlib import Path as _Path

APP = _Path(_os.environ.get('APP_DIR', '/app'))
PLAN = __PLAN__


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([_json.dumps(payload)] if payload is not None else [])
    result = _json.loads(_subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        _time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError('purchase failed: %s' % r)
    return r['id'], r['result']


def main():
    catalog = _json.loads((APP/'catalog.json').read_text())
    ids, results = [], []
    for request in PLAN:
        i, res = buy(request)
        ids.append(i)
        results.append(res)
    out = analyse(catalog, results)
    submission = {'main_phase': out['main_phase'], 'minority_phase': out['minority_phase'],
                  'frequency_scale': out['frequency_scale'], 'c_axis_lab': out['c_axis_lab'],
                  'faceB_strongest_mode_ratio_tt_over_zz': out['ratio_tt_over_zz'],
                  'faceB_highest_optical_phonon_cm1': out['highest_optical_cm1'], 'evidence_ids': ids}
    (APP/'submission.json').write_text(_json.dumps(submission, indent=1) + '\n')
    chi = ', '.join('%s %.2f' % (m, v) for m, v in out['chi2_per_point'].items())
    mt = ', '.join('%s %.3f' % (m, v['chi2_drop_fraction']) for m, v in out['minority_test'].items())
    (APP/'report.md').write_text("""# SnF2 crystal: reference analysis

**Calibration.** The strongest line and its anti-Stokes partner are not symmetric about zero on the spectrometer
axis: the axis is offset by %.2f cm-1 (laser line away from its nominal wavelength). The offset is fitted together
with the global frequency scale s = omega_sample/omega_PBEsol = %.5f; a pure scale through the origin would absorb
the offset and bias s by about +0.6 %%.

**Model.** Every purchased window is fitted at once with a forward model per candidate polymorph: lab-frame Raman
tensors Q R Q^T of the database modes, Stokes/anti-Stokes factors omega_s^3 (n+1)/omega and omega_as^3 n/omega at
295 K, pseudo-Voigt lines (instrumental 1.2 cm-1), linear background and the 1.5 %% power drift. Polar modes are
LO-mixed for the phonon wavevector along the face normal using the Born charges and eps_inf; the archived Born
charges list only the symmetry-independent atoms and were expanded with the space-group operations (acoustic sum
rule checked). Degenerate partners are summed. Orientations were searched globally (areas of isolated non-polar
lines) and refined; rotated polarizer/analyzer settings (45, 135, 30 degrees) remove the mirror ambiguity left by
parallel/crossed settings along the lab axes. Face-A backscattering cannot distinguish Q from a 180-degree turn
about the face normal; the face-B predictions are invariant under that turn.

**Polymorph.** chi2 per point: %s. The orthorhombic beta phase (mp-7456) fits; the tetragonal gamma phase has
nearly the same strong lines but not the split doublets, the low-frequency mode positions or the polarization
pattern; alpha only matches by line density.

**Minority phase.** Fractional chi2 drop when an orientation-averaged component is added: %s -> %s.

**Orientation and predictions.** c axis (lab) = %s. Face B: I(t,t)/I(z,z) of the strongest mode = %.3f; highest
optical phonon for q along n_B = %.1f cm-1 (LO, Born charges and eps_inf, scaled by s). Uncertainties: c axis about
1-2 degrees, s about 0.05 %%, LO frequency about 0.3 %%; the electro-optic contribution to LO Raman tensors is not
available from the database and is neglected, and the absolute scattering efficiencies are not tested.
""" % (out['axis_offset_cm1'], out['frequency_scale'], chi, mt, out['minority_phase'], [round(v, 4) for v in out['c_axis_lab']],
       out['ratio_tt_over_zz'], out['highest_optical_cm1']))
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(_json.dumps(status))
            return
        _time.sleep(2)


if __name__ == '__main__':
    main()
'''


def main():
    skip_variants = '--skip-variants' in sys.argv
    pool = M.load_pool()
    recs = {mpid: pool[key] for mpid, key in M.POLYMORPHS.items()}
    parsed = {mpid: M.record(rec) for mpid, rec in recs.items()}
    main_e, minor_e = parsed[HIDDEN['main_phase']], parsed[HIDDEN['minority_phase']]
    Q = M.euler(*HIDDEN['euler_zyz_deg'])
    lab_model = build_lab_model(main_e, minor_e, Q, HIDDEN)
    tg = M.targets(main_e, Q, HIDDEN['frequency_scale'])
    reference = {'main_phase': HIDDEN['main_phase'], 'minority_phase': HIDDEN['minority_phase'],
                 'frequency_scale': HIDDEN['frequency_scale'], 'c_axis_lab': tg['c_axis_lab'],
                 'faceB_ratio_tt_over_zz': tg['ratio_tt_over_zz'], 'faceB_highest_optical_cm1': tg['highest_optical_cm1']}
    hidden = {'thresholds': THRESHOLDS, 'reference': reference, 'sample': HIDDEN,
              'axis_offset_cm1': lab_model['axis_offset_cm1'], 'orientation_matrix': Q.tolist()}
    hidden_values = [HIDDEN['frequency_scale'], *tg['c_axis_lab'], tg['ratio_tt_over_zz'], tg['highest_optical_cm1'],
                     lab_model['axis_offset_cm1'], HIDDEN['laser_true_nm']]
    catalog = {
        'sample': ('Single crystal of SnF2 (composition known) with one polished face, face A. The spectra served by the lab '
                   'are simulated from a physical forward model of this hidden sample, built from the database records '
                   'below; treat them as measurements.'),
        'lab_frame': ('z is the outward normal of face A; the laser is incident along -z and the backscattered light is '
                      'collected along +z, so the phonon wavevector is along z. Polarizer and analyzer angles are measured '
                      'in the xy plane from +x.'),
        'face_B': {'outward_normal_lab': [round(math.cos(math.radians(M.FACE_B_DEG)), 6), round(math.sin(math.radians(M.FACE_B_DEG)), 6), 0.0],
                   'horizontal_polarization_t_lab': [round(-math.sin(math.radians(M.FACE_B_DEG)), 6), round(math.cos(math.radians(M.FACE_B_DEG)), 6), 0.0],
                   'note': 'face B is not polished and cannot be measured'},
        'instrument': {
            'laser_nominal_wavelength_nm': 532.0,
            'raman_shift_axis': 'computed by the spectrometer software from the nominal laser wavelength',
            'geometry': 'backscattering at normal incidence on face A, microscope objective, confocal',
            'polarization': ('polarizer_deg sets the incident polarization, analyzer_deg the analysed scattered polarization; '
                             'analyzer_deg null removes the analyzer (total scattered intensity). The crystal is treated as a '
                             'thin plate: the fields inside it follow the polarizer and analyzer directions (birefringence and '
                             'refraction neglected). The spectrometer throughput does not depend on polarization.'),
            'window': 'each request records one window of 320 cm-1 (centre +/- 160 cm-1) with 0.4 cm-1 pixels; centre in [-450, 650] cm-1',
            'instrumental_profile': 'Gaussian, FWHM 1.2 cm-1',
            'notch_filter': 'pixels within 18 cm-1 of the laser line are blocked and omitted',
            'detector': ('CCD: offset 500 counts and read noise 5 counts rms per accumulation, photon (Poisson) statistics, '
                         'flat spectral response; 60 s per accumulation'),
            'laser': 'constant power with 1.5 % rms power fluctuation between requests; laser heating is negligible',
            'stage_temperature_K': 295},
        'request': {'kind': 'raman', 'polarizer_deg': '<number>', 'analyzer_deg': '<number or null>',
                    'center_cm1': '<number>', 'accumulations': '<integer 1-4, default 1>'},
        'costs': 'one unit per accumulation', 'budget': BUDGET,
        'database': {
            'citation': ('M. Bagheri and H.-P. Komsa, High-throughput computation of Raman spectra from first principles, '
                         'Sci. Data 10, 80 (2023); Computational Raman Database v1.0.1, Materials Cloud Archive record '
                         '0m9r0-g3486, DOI 10.24435/materialscloud:pg-h3 (CC BY 4.0). Phonons from phonondb (PBEsol, VASP).'),
            'conventions': ('Records are the archived values (floats rounded to 6 significant digits). Frequencies in cm-1 at '
                            'the Gamma point without LO-TO splitting. eigenvectors[mode][atom][xyz] are normalized eigenvectors '
                            'of the mass-weighted dynamical matrix. raman_tensors[mode] (null if inactive) is '
                            'R = (V_c/4 pi) d chi/d xi, xi the normal coordinate along eigenvector/sqrt(mass); the Stokes '
                            'cross-section of Bagheri & Komsa is proportional to omega_s^4 (n+1)/omega |e_s.R.e_i|^2. '
                            'born_charges is archived in phonopy BORN form (each 3x3 row: first index electric field, second '
                            'atomic displacement); epsilon_inf is the electronic dielectric tensor. Cartesian frame and cell '
                            'as in cell_A (rows are lattice vectors, Angstrom).')},
        'polymorph_records': {mpid: catalog_record(parsed[mpid], recs[mpid]) for mpid in ('mp-2843', 'mp-7456', 'mp-7457')},
    }
    for mpid in catalog['polymorph_records']:
        catalog['polymorph_records'][mpid]['database_id'] = M.POLYMORPHS[mpid]
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (LAB/'private').mkdir(parents=True, exist_ok=True)
    (LAB/'private/lab_model.json').write_text(json.dumps(lab_model) + '\n')
    (LAB/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (LAB/'private/hidden_values.json').write_text(json.dumps(hidden_values) + '\n')
    public = {'task_id': 'raman-polymorph-orientation', 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
              'requests': {'raman': {'kind': 'raman', 'polarizer_deg': '<deg>', 'analyzer_deg': '<deg or null>',
                                     'center_cm1': '<cm-1>', 'accumulations': '<1-4>', 'cost': 'accumulations'}},
              'catalog': '/app/catalog.json'}
    (LAB/'public.json').write_text(json.dumps(public, indent=1) + '\n')
    (TASK/'tests/public.json').write_text(json.dumps({'task_id': 'raman-polymorph-orientation', 'solver_budget': BUDGET,
                                                      'min_report_chars': 400}, indent=1) + '\n')
    analysis_src = (HERE/'analysis.py').read_text()
    solve = analysis_src.replace('"""Analysis library for raman-polymorph-orientation (reference solution and construction variants).',
                                 '"""Reference solution for raman-polymorph-orientation (generated by make_task.py from analysis.py).', 1)
    (TASK/'solution/solve.py').write_text(solve + DRIVER.replace('__PLAN__', json.dumps(PLAN).replace('null', 'None')))
    texts = [(TASK/'instruction.md').read_text(), (LAB/'public.json').read_text(), (TASK/'tests/public.json').read_text(),
             (TASK/'environment/agent/catalog.json').read_text()]
    found = leaks(hidden_values, texts)
    print('targets', json.dumps(reference), 'axis offset', lab_model['axis_offset_cm1'])
    print('leaks', found)
    print('plan units', sum(r['accumulations'] for r in PLAN), 'evaluations', len(PLAN), 'lab-axes-only spectra', LAB_AXES_ONLY)
    if skip_variants:
        return
    import variants
    variants.run(catalog, replay(load_task_lab(), lab_model, PLAN), reference, THRESHOLDS, LAB_AXES_ONLY, TASK, HERE)


if __name__ == '__main__':
    main()
