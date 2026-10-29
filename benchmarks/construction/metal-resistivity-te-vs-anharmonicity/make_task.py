"""Generate the staged, lab, private and ablation data for metal-resistivity-te-vs-anharmonicity from pool.json.

    bash run_in_image.sh make_task.py            (about 1 min)

Writes environment/agent/catalog.json, environment/lab/{data.json, phonons.py}, environment/lab/private/{hidden.json,
hidden_values.json}, solution/ablation_values.json, the budget in both public.json files, and calibration.json here.
"""
import json
import shutil
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/metal-resistivity-te-vs-anharmonicity'
COSTS = {'transport': 2, 'phonons': 1}
BUDGET = 46
WITHHELD = ['Pb/4-epw/2-with_TE/610K', 'Pb/4-epw/3-anharm/210K', 'Pb/4-epw/3-anharm/410K', 'Pb/4-epw/3-anharm/610K',
            'Pb/4-epw/4-TE_and_anharm/410K', 'Pb/4-epw/4-TE_and_anharm/510K', 'Pb/4-epw/4-TE_and_anharm/610K',
            'Nb/4-epw/3-TE_and_anharm/2000K']
NO_RESULTS = ['Nb/4-epw/3-TE_and_anharm/1500K']
THRESHOLDS = {'al_te_change_abs_pp': 1.0, 'pb_te_610_rel': 0.06, 'pb_anharm_abs_pp': 3.0, 'pb_te_anharm_510_rel': 0.05,
              'nb_te_anharm_2000_rel': 0.055}


def condition(run_id):
    m, _, folder = run_id.split('/')[:3]
    rest = run_id.split('/')[3:]
    if folder == '1-DFPT':
        return ('harmonic: DFPT force constants at the 0 K lattice parameter (archive: "EPW calculation of %s without TE '
                'or TDEP")' % m)
    T = rest[0]
    if folder == '2-with_TE':
        return 'thermal expansion: DFPT force constants at the lattice parameter of %s ("with lattice parameter at %s")' % (T, T)
    if folder == '3-anharm':
        return 'anharmonic at the 0 K lattice: TDEP force constants sampled at %s, 0 K lattice parameter' % T
    return 'thermal expansion + anharmonicity: TDEP force constants sampled at %s, lattice parameter of %s' % (T, T)


def ifc_view(rec):
    v = {'ifc': rec['id'], 'metal': rec['metal'], 'kind': rec['kind'], 'archive_file': rec['archive_file'],
         'lattice_temperature_K': rec['lattice_temperature_K'], 'sampling_temperature_K': rec['sampling_temperature_K']}
    if rec['kind'] == 'TDEP':
        v['tdep_iteration'] = rec['tdep_iteration']
    return v


def main():
    runs, ifcs = M.RUNS, M.IFCS
    # ---------------- catalog (agent-visible) ----------------
    lattice = {}
    for r in runs.values():
        m = r['metal']
        if r['id'].endswith('1-DFPT/epw') or '1-DFPT' in r['id']:
            key = '0K'
        else:
            key = r['id'].split('/')[3]
        lattice.setdefault(m, {})[key] = round(r['lattice_parameter_bohr'], 6)
    transport = []
    for rid, r in runs.items():
        temps = [float(t) for t in r['ibte_in']['temps'].replace(',', ' ').split()]
        entry = {'run': rid, 'metal': r['metal'], 'condition': condition(rid), 'temperatures_K': temps,
                 'archive_files': r['archive_files']}
        if rid in WITHHELD:
            entry['served'] = False
            entry['note'] = 'withheld by this task (its result is to be estimated)'
        elif rid in NO_RESULTS:
            entry['served'] = False
            entry['note'] = 'not purchasable: the archived ibte.out of this run repeats the input and contains no results'
        else:
            entry['served'] = True
        transport.append(entry)
    readme = (Path('/data/README.md').read_text() if Path('/data/README.md').exists() else '')
    catalog = {
        'title': 'Phonon-limited resistivity of Pb, Nb and Al with thermal expansion and TDEP anharmonicity',
        'source': 'A. Wang, J. Yin, F. A. Goudreault, M. Cote, O. Hellman, S. Ponce (2026), Materials Cloud Archive '
                  'record zxf0b-rs406 (DOI 10.24435/materialscloud:4h-r1), CC BY 4.0. Retrieved author calculations; '
                  'no new EPW calculation is run.',
        'archive_readme': readme,
        'methods': {
            'electrons': 'Quantum ESPRESSO 7.3 (Pb: fully relativistic LDA with spin-orbit coupling; Al, Nb: see the '
                         'archived inputs). One atom per primitive cell: Pb and Al fcc, Nb bcc.',
            'phonons': 'DFPT on an 8x8x8 (Pb) or 6x6x6 (Al, Nb) q-grid, force constants from q2r; TDEP (stochastic '
                       'self-consistent sampling with quantum statistics at the stated temperature, see the archive README) '
                       'converted to the same q2r format. The electron-phonon deformation potentials always come from DFPT at the run\'s '
                       'lattice parameter; the force constants (DFPT or TDEP) set the phonon frequencies and eigenvectors.',
            'transport': 'EPW (v5.8) linearised Boltzmann equation for the phonon-limited conductivity on fine k and q '
                         'meshes, with adaptive smearing; each output prints the self-energy relaxation-time (SERTA) '
                         'solution followed by the iterations of the full iterative solution (IBTE).',
            'lattice_parameters_bohr': lattice,
            'lattice_parameter_note': 'Temperature-dependent lattice parameters from the authors (thermal expansion); '
                                      'the 0 K value is the static relaxed lattice.'},
        'requests': {
            'transport': {'format': {'kind': 'transport', 'run': '<run with served = true>'}, 'cost': COSTS['transport'],
                          'returns': 'the run\'s input parameters, lattice parameter, fine meshes and, for each '
                                     'temperature, the Fermi level, printed DOS, SERTA conductivity tensor, IBTE '
                                     'iteration history and final IBTE conductivity tensor, in (Ohm cm)^-1'},
            'phonons': {'format': {'kind': 'phonons', 'ifc': '<ifc id>', 'grid': '<n, 4 to 20>',
                                   'eigenvectors': '<optional, default false>'}, 'cost': COSTS['phonons'],
                        'returns': 'harmonic frequencies (meV) of that force-constant set on a Gamma-centred n^3 grid '
                                   '(q folded into the first Brillouin zone, Cartesian, units of 2 pi / alat), '
                                   'optionally the polarisation vectors, a DOS histogram, standard DOS moments and the '
                                   'dispersion along high-symmetry lines; computed by the lab (numpy Fourier '
                                   'interpolation with the simple acoustic sum rule, identical to matdyn.x)'}},
        'budget': BUDGET,
        'transport_runs': transport,
        'force_constant_sets': [ifc_view(r) for r in ifcs.values()],
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    # ---------------- lab data ----------------
    served = {rid: {**r, 'condition': condition(rid)} for rid, r in runs.items() if rid not in WITHHELD + NO_RESULTS}
    lab = {'costs': COSTS, 'runs': served, 'ifcs': ifcs, 'withheld': WITHHELD, 'no_results': NO_RESULTS}
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab) + '\n')
    shutil.copy(HERE/'phonons.py', TASK/'environment/lab/phonons.py')
    # ---------------- hidden truth ----------------
    t = M.truths()
    withheld_rho = {rid: [M.rho(rid, row['T_K']) for row in runs[rid]['results']] for rid in WITHHELD}
    hidden = {'targets': t, 'thresholds': THRESHOLDS, 'withheld_resistivity_uohm_cm': withheld_rho}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [t['al_te_change_percent_481K'], t['pb_te_resistivity_610K'],
              *t['pb_anharmonic_change_percent_0K_lattice'].values(), t['pb_te_anharm_resistivity_510K'],
              t['nb_te_anharm_resistivity_2000K'], *[v for vs in withheld_rho.values() for v in vs]]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    # ---------------- ablations and robustness ----------------
    oracle = {'al': M.al_te_change_481(), 'te610': M.pb_te_610(), 'anh': {str(T): M.pb_anh0K(T) for T in (410, 610)},
              'ta510': M.pb_te_anh_510(), 'nb': M.nb_te_anh_2000()}
    serta = {'al': M.al_te_change_481(tensor='serta'), 'te610': M.pb_te_610(tensor='serta'),
             'anh': {str(T): M.pb_anh0K(T, tensor='serta') for T in (410, 610)},
             'ta510': M.pb_te_anh_510(tensor='serta'), 'nb': M.nb_te_anh_2000(tensor='serta')}
    abl = {
        'al_reference_mesh70': M.al_te_change_481(ref='Al/4-epw/1-DFPT/ibte.out-mesh70'),
        'al_reference_mesh60': M.al_te_change_481(ref='Al/4-epw/1-DFPT/ibte.out-mesh60'),
        'te_anh_510_separable_multiplicative': M.pb_te_anh_510('separable_multiplicative'),
        'te_anh_510_separable_additive': M.pb_te_anh_510('separable_additive'),
        'anharm_from_expanded_lattice': {str(T): M.pb_anh0K(T, mode='T_lattice_transfer') for T in (410, 610)},
        'anharm_inverse_square_moment': {str(T): M.pb_anh0K(T, weight='inv2_moment', calibrate=False) for T in (410, 610)},
        'te_anh_510_inverse_square_moment': M.pb_te_anh_510(weight='inv2_moment'),
        'nb_cancellation': M.nb_te_anh_2000('cancellation'),
        'nb_pb_calibrated_surrogate': M.nb_te_anh_2000('surrogate_pb_calibrated'),
        'nb_uncalibrated_surrogate': M.nb_te_anh_2000('surrogate_uncalibrated'),
        'nb_constant_anharmonic_fraction': M.nb_te_anh_2000('constant_anharmonic_fraction'),
        'te_610_linear_extrapolation': M.pb_te_610('rho_linear'),
        'te_610_phonon_only': M.pb_te_610('phonon_only'),
        'serta': serta,
        'te_anh_510_linear_fit_with_10K': M.pb_te_anh_510('linear_fit_with_10K'),
        'calibration_with_10K': {'anh': {str(T): M.pb_anh0K(T, calib_T=(10, 110, 210, 310)) for T in (410, 610)},
                                 'ta510': M.pb_te_anh_510(calib_T=(10, 110, 210, 310))},
    }
    rob = {
        'longitudinal_q2': {'anh': {str(T): M.pb_anh0K(T, weight='long_q2') for T in (410, 610)},
                            'ta510': M.pb_te_anh_510(weight='long_q2'), 'nb': M.nb_te_anh_2000(weight='long_q2'),
                            'te610': M.pb_te_610(weight='long_q2')},
        'calibration_110_310': {'anh': {str(T): M.pb_anh0K(T, calib_T=(110, 210, 310)) for T in (410, 610)},
                                'ta510': M.pb_te_anh_510(calib_T=(110, 210, 310))},
        'uncalibrated_q2_pb': {'anh': {str(T): M.pb_anh0K(T, calibrate=False) for T in (410, 610)}},
        'te_610_ratio_quadratic': M.pb_te_610('R_quadratic'),
        'te_610_rho_cubic': M.pb_te_610('rho_cubic'),
        'te_610_residual_constant': M.pb_te_610('phonon_e_const'),
        'te_610_rho_quadratic': M.pb_te_610('rho_quadratic'),
        'te_610_uniform_weight_with_residual': M.pb_te_610(weight='uniform'),
        'nb_ratio_to_harmonic': M.nb_te_anh_2000('ratio_to_harmonic'),
        'nb_proportional_T': M.nb_te_anh_2000('proportional_T'),
        'te_anh_510_ratio_to_harmonic_linear': M.pb_te_anh_510('R_harm_linear'),
        'te_anh_510_linear_fit_110_310': M.pb_te_anh_510('linear_fit'),
        'al_linear_interpolation': M.al_te_change_481(order=1),
    }
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'oracle_model': oracle, 'ablations': abl,
                                                                  'robustness': rob}, indent=1) + '\n')
    for rel in ('environment/lab/public.json', 'tests/public.json'):
        p = json.loads((TASK/rel).read_text())
        p['solver_budget'] = BUDGET
        (TASK/rel).write_text(json.dumps(p, indent=1) + '\n')
    (HERE/'calibration.json').write_text(json.dumps({'truth': t, 'oracle_model': oracle, 'ablations': abl,
                                                     'robustness': rob, 'thresholds': THRESHOLDS}, indent=1) + '\n')
    rel = lambda v, ref: round(100 * (v / ref - 1), 2)
    print('truth', json.dumps(t))
    print('oracle', {'al': round(oracle['al'], 3), 'te610 %': rel(oracle['te610'], t['pb_te_resistivity_610K']),
                     'anh': {k: round(v, 3) for k, v in oracle['anh'].items()},
                     'ta510 %': rel(oracle['ta510'], t['pb_te_anharm_resistivity_510K']),
                     'nb %': rel(oracle['nb'], t['nb_te_anharm_resistivity_2000K'])})
    print('serta', {'al': round(serta['al'], 3), 'te610 %': rel(serta['te610'], t['pb_te_resistivity_610K']),
                    'anh': {k: round(v, 3) for k, v in serta['anh'].items()},
                    'ta510 %': rel(serta['ta510'], t['pb_te_anharm_resistivity_510K']),
                    'nb %': rel(serta['nb'], t['nb_te_anharm_resistivity_2000K'])})
    print('long_q2 te610 %', rel(rob['longitudinal_q2']['te610'], t['pb_te_resistivity_610K']),
          'ta510 cal110-310 %', rel(rob['calibration_110_310']['ta510'], t['pb_te_anharm_resistivity_510K']))
    print('catalog bytes', (TASK/'environment/agent/catalog.json').stat().st_size,
          'lab data bytes', (TASK/'environment/lab/data.json').stat().st_size)
    print('catalog units', sum(COSTS['transport'] for e in transport if e['served']) + COSTS['phonons'] * len(ifcs))


if __name__ == '__main__':
    main()
