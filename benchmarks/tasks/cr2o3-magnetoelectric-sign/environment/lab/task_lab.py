"""Lab logic for cr2o3-magnetoelectric-sign: retrieval of author-archived VASP/Elk/phonopy/QE results and grading.

Purchases return archived files (parsed into JSON) of Materials Cloud record xd06z-0pp29. No new calculation is run.
The magnetoelectric coefficient, Born charges, moment derivatives and domain assignments are never served; they
grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-archived result (no new calculation)'
MODES = ['Eu_1_x', 'Eu_2_x', 'Eu_3_x', 'Eu_4_x']
ATOMS = ['Cr1', 'Cr2', 'Cr3', 'Cr4', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6']
KINDS = {'born_reference': {'kind'}, 'born_displacement': {'kind', 'amplitude', 'atom'},
         'elk_moment': {'kind', 'mode', 'amplitude'}, 'vasp_moment': {'kind', 'mode', 'amplitude'},
         'phonon_gamma': {'kind', 'masses'}, 'force_constant_eigenmodes': {'kind'}, 'structure': {'kind', 'key'},
         'qe_inputs': {'kind', 'folder'}}


def validate(request):
    kind = request.get('kind')
    if kind not in KINDS:
        raise ValueError('kind must be one of %s' % sorted(KINDS))
    if set(request) != KINDS[kind]:
        raise ValueError('%s takes exactly the keys %s' % (kind, sorted(KINDS[kind])))
    if kind == 'born_displacement':
        if request['amplitude'] not in DATA['born']:
            raise ValueError('amplitude must be one of %s (use born_reference for 0_0)' % sorted(DATA['born']))
        if request['atom'] not in ATOMS:
            raise ValueError('atom must be one of %s' % ATOMS)
    elif kind in ('elk_moment', 'vasp_moment'):
        table = DATA['elk' if kind == 'elk_moment' else 'vasp']
        if request['mode'] not in table:
            raise ValueError('mode must be one of %s' % MODES)
        if request['amplitude'] not in table[request['mode']]:
            raise ValueError('amplitude must be one of %s for this mode' % sorted(table[request['mode']]))
    elif kind == 'phonon_gamma':
        if request['masses'] not in DATA['phonon_gamma']:
            raise ValueError('masses must be one of %s' % sorted(DATA['phonon_gamma']))
    elif kind == 'structure':
        if request['key'] not in DATA['structures']:
            raise ValueError('key must be one of %s' % sorted(DATA['structures']))
    elif kind == 'qe_inputs':
        if request['folder'] not in DATA['qe']:
            raise ValueError('folder must be one of %s' % sorted(DATA['qe']))
    return 1


def run(request, workdir):
    kind = request['kind']
    if kind == 'born_reference':
        return {'folder': 'VASP_calculations/Born_effective_charges/Disp_0_0', **DATA['born_reference'], 'source': SOURCE}
    if kind == 'born_displacement':
        a, atom = request['amplitude'], request['atom']
        return {'folders': {c: 'VASP_calculations/Born_effective_charges/Disp_%s/Pos_%s_%s' % (a, atom, c) for c in 'xyz'},
                'runs': DATA['born'][a][atom], 'source': SOURCE}
    if kind == 'elk_moment':
        m, a = request['mode'], request['amplitude']
        return {'folder': 'Elk_calculations/Induced_magnetic_moments/%s/Amplitude_%s/Calc_force' % (m, a),
                **DATA['elk'][m][a], 'source': SOURCE}
    if kind == 'vasp_moment':
        m, a = request['mode'], request['amplitude']
        return {'folder': 'VASP_calculations/Induced_magnetic_moments/%s/Amplitude_%s' % (m, a),
                **DATA['vasp'][m][a], **DATA['vasp_moment_inputs'], 'source': SOURCE}
    if kind == 'phonon_gamma':
        return {'file': 'Phonon_calculation/Getting_eigenvectors_force_constant_matrix/band_%s.yaml' % request['masses'],
                'q_position': [0.0, 0.0, 0.0], **DATA['phonon_gamma'][request['masses']],
                'eigenvector_format': 'phonopy band.yaml eigenvectors (real parts; imaginary parts are zero), '
                                      'indexed [band][atom][cartesian]', 'source': SOURCE}
    if kind == 'force_constant_eigenmodes':
        return {'folder': 'Phonon_calculation/Getting_eigenvectors_force_constant_matrix', **DATA['fc_eigen'],
                'note': 'Dataset names and values exactly as stored in the two hdf5 files; the files carry no units.',
                'source': SOURCE}
    if kind == 'structure':
        return {'key': request['key'], **DATA['structures'][request['key']], 'source': SOURCE}
    return {'folder': 'Quantum_espresso/' + request['folder'], 'files': DATA['qe'][request['folder']], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    ref, t = PRIVATE['reference'], PRIVATE['thresholds']
    fields = ['alpha_perp_ps_per_m', 'mode_fractions', 'dominant_mode', 'dM_du_out_pointing_muB_per_A',
              'born_charge_Cr1', 'dataset_domain', 'domain_selected_by_parallel_field_annealing', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    al, fr, dm, zb, dom = (payload['alpha_perp_ps_per_m'], payload['mode_fractions'],
                           payload['dM_du_out_pointing_muB_per_A'], payload['born_charge_Cr1'], payload['dataset_domain'])
    def numdict(d, keys):
        return isinstance(d, dict) and set(d) == set(keys) and all(finite(v) for v in d.values())
    if not numdict(al, ['out_pointing', 'in_pointing']):
        return False, {'error': 'alpha_perp_ps_per_m needs finite "out_pointing" and "in_pointing" values'}
    if not numdict(fr, MODES):
        return False, {'error': 'mode_fractions needs a finite value for each of %s' % MODES}
    if not numdict(dm, ['Eu_3_x', 'Eu_4_x']):
        return False, {'error': 'dM_du_out_pointing_muB_per_A needs finite "Eu_3_x" and "Eu_4_x" values'}
    if not numdict(zb, ['xx', 'xy', 'zz']):
        return False, {'error': 'born_charge_Cr1 needs finite "xx", "xy" and "zz" values'}
    if not (isinstance(dom, dict) and set(dom) == {'VASP', 'Elk'} and all(isinstance(v, str) for v in dom.values())):
        return False, {'error': 'dataset_domain needs "VASP" and "Elk" entries ("out_pointing" or "in_pointing")'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    a0 = ref['alpha_out_pointing_ps_m']
    norm = lambda s: s.strip().lower().replace('-', '_').replace(' ', '_') if isinstance(s, str) else s
    errs = {'alpha_out_rel': al['out_pointing'] / a0 - 1, 'alpha_in_rel': al['in_pointing'] / (-a0) - 1,
            'fractions_abs': {m: fr[m] - ref['mode_fractions'][m] for m in MODES},
            'dmdu_rel': {m: dm[m] / ref['dmdu_out_pointing_muB_per_A'][m] - 1 for m in dm},
            'born_abs': {k: zb[k] - ref['born_Cr1'][k] for k in zb}}
    checks = {'evidence_executed': evidence,
              'alpha_out_pointing': abs(errs['alpha_out_rel']) <= t['alpha_rel_tol'],
              'alpha_in_pointing': abs(errs['alpha_in_rel']) <= t['alpha_rel_tol'],
              'mode_fractions': all(abs(v) <= t['fraction_abs_tol'] for v in errs['fractions_abs'].values()),
              'dominant_mode': (isinstance(payload['dominant_mode'], str)
                                and norm(payload['dominant_mode']).removesuffix('_x') == norm(ref['dominant_mode']).removesuffix('_x')),
              'dM_du': all(abs(v) <= t['dmdu_rel_tol'] for v in errs['dmdu_rel'].values()),
              'born_charge_Cr1': (abs(errs['born_abs']['xx']) <= t['born_diag_abs_tol']
                                  and abs(errs['born_abs']['zz']) <= t['born_diag_abs_tol']
                                  and abs(errs['born_abs']['xy']) <= t['born_offdiag_abs_tol']),
              'dataset_domain': all(norm(dom[k]) == ref['dataset_domain'][k] for k in ('VASP', 'Elk')),
              'annealing_domain': norm(payload['domain_selected_by_parallel_field_annealing'])
                                  == ref['annealed_domain_parallel_EH']}
    return all(checks.values()), {'checks': checks, 'errors': errs}
