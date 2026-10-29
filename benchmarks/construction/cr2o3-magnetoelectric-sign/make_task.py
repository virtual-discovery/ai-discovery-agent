"""Write the catalog, lab data, private targets and ablation values for cr2o3-magnetoelectric-sign from pool.json.

    python3 -B make_task.py
"""
import json
from pathlib import Path

import numpy as np

import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/cr2o3-magnetoelectric-sign'
P = M.POOL
MODES = M.MODES
ELK_TAGS = ['0_005', '-0_005', '0_01', '-0_01', '0_02', '-0_02']
BUDGET = 44
COSTS = {'born_displacement': 1, 'born_reference': 1, 'elk_moment': 1, 'vasp_moment': 1, 'phonon_gamma': 1,
         'force_constant_eigenmodes': 1, 'structure': 1, 'qe_inputs': 1}
THRESH = {'alpha_rel_tol': 0.05, 'fraction_abs_tol': 0.025, 'dmdu_rel_tol': 0.05,
          'born_diag_abs_tol': 0.08, 'born_offdiag_abs_tol': 0.05}


def r(x, n=6):
    return float(f'{x:.{n}g}')


def lab_data():
    born = {}
    for tag, recs in P['born'].items():
        if tag == '0_0':
            continue
        born[tag] = {}
        for key, rec in recs.items():
            atom, c = key.split('_')
            born[tag].setdefault(atom, {})[c] = {
                'POSCAR': rec['poscar'], 'OUTCAR_ion_positions_cartesian_A': rec['outcar_positions_cartesian_A'],
                'OUTCAR_dipole_lines': rec['dipole_lines'], 'OUTCAR_energy_sigma0_eV': rec['energy_sigma0_eV'],
                'OUTCAR_finished': rec['finished']}
    ref = P['born']['0_0']['reference']
    born_ref = {'POSCAR': ref['poscar'], 'OUTCAR_ion_positions_cartesian_A': ref['outcar_positions_cartesian_A'],
                'OUTCAR_dipole_lines': ref['dipole_lines'], 'OUTCAR_energy_sigma0_eV': ref['energy_sigma0_eV'],
                'OUTCAR_pseudopotential_lines': P['born']['0_0']['outcar_header'],
                'INCAR': P['structures']['born_INCAR']}
    elk = {m: {t: {'elk.in': rec['elk_in'], 'INFO.OUT_final_moments_block': rec['moments_block'],
                   'INFO.OUT_scf_loops': rec['scf_loops'], 'INFO.OUT_messages': rec['messages'],
                   'INFO.OUT_last_rms_change_KS_potential': rec['last_rms_potential_change'],
                   'INFO.OUT_rms_target': rec['rms_target'],
                   'INFO.OUT_final_total_energy_Ha': rec['final_total_energy_Ha'],
                   'INFO.OUT_run_completed': rec['run_completed'],
                   'other_files_in_folder_last_lines': rec['other_files_in_folder']}
               for t, rec in P['elk'][m].items()} for m in MODES}
    vasp = {}
    for m in MODES:
        vasp[m] = {}
        for t, rec in P['vasp_moments'][m].items():
            out = {'POSCAR': rec['poscar'], 'OUTCAR_final_total_magnetization_line': rec['total_magnetization_line'],
                   'OUTCAR_electronic_steps': rec['electronic_steps'], 'OUTCAR_finished_normally': rec['finished_normally'],
                   'OUTCAR_last_energy_sigma0_eV': rec['last_energy_sigma0_eV']}
            for c in 'xyz':
                out[f'OUTCAR_last_ion_table_magnetization_{c}'] = rec[f'ion_magnetization_{c}']
                if f'ion_orbital_moment_{c}' in rec:
                    out[f'OUTCAR_last_ion_table_orbital_moment_{c}'] = rec[f'ion_orbital_moment_{c}']
            vasp[m][t] = out
    s = P['structures']
    structures = {'relaxation': {'POSCAR': s['relaxation_POSCAR'], 'CONTCAR': s['relaxation_CONTCAR'],
                                 'INCAR': s['relaxation_INCAR']},
                  'phonon': {'POSCAR': s['phonon_POSCAR'], 'INCAR': s['phonon_INCAR'], 'FORCE_SETS': s['FORCE_SETS'],
                             'bands.conf': s['bands_conf']}}
    fc = P['fc_eigen']
    return {'born': born, 'born_reference': born_ref, 'elk': elk, 'vasp': vasp,
            'vasp_moment_inputs': {'INCAR': P['vasp_moments']['INCAR'], 'KPOINTS': P['vasp_moments']['KPOINTS']},
            'phonon_gamma': P['phonon_gamma'], 'fc_eigen': {'Cr_ph_eigenval_pol.hdf5': fc['values'],
                                                            'Cr_ph_eigenvec_pol.hdf5': fc['vectors']},
            'structures': structures, 'qe': P['qe']}


def labels(tags):
    return sorted(tags, key=M.amp)


def catalog(data):
    return {
        'material': 'Cr2O3 (corundum, R-3c), antiferromagnetic, 10-atom rhombohedral cell (4 Cr, 6 O).',
        'source': 'Bousquet et al., J. Phys.: Condens. Matter (2024), doi 10.1088/1361-648X/ad1a59; Materials Cloud '
                  'Archive record xd06z-0pp29 (doi 10.24435/materialscloud:ek-fp), CC BY 4.0. Retrieved archived results, '
                  'no new calculations.',
        'archive_README': P['readme'].replace('https://doi.org/10.48550/arXiv.2309.02095', 'arXiv:2309.02095'),
        'budget': BUDGET,
        'requests': {
            'born_reference': {'cost': 1, 'request': {'kind': 'born_reference'},
                               'returns': 'Born_effective_charges/Disp_0_0: POSCAR, INCAR, OUTCAR ion positions, dipole '
                                          'lines, energy and pseudopotential header lines'},
            'born_displacement': {'cost': 1, 'request': {'kind': 'born_displacement', 'amplitude': '<folder tag>',
                                                         'atom': '<Cr1..Cr4, O1..O6>'},
                                  'returns': 'the three folders Disp_<tag>/Pos_<atom>_{x,y,z}: POSCAR, OUTCAR ion positions, '
                                             'dipole lines (LCALCPOL), energy',
                                  'amplitudes': labels(data['born']), 'atoms': M.NAMES},
            'elk_moment': {'cost': 1, 'request': {'kind': 'elk_moment', 'mode': '<Eu_1_x..Eu_4_x>', 'amplitude': '<tag>'},
                           'returns': 'Elk_calculations/Induced_magnetic_moments/<mode>/Amplitude_<tag>/Calc_force: elk.in, '
                                      'final Moments block of INFO.OUT, SCF information, other files in the folder',
                           'amplitudes': {m: labels(data['elk'][m]) for m in MODES}},
            'vasp_moment': {'cost': 1, 'request': {'kind': 'vasp_moment', 'mode': '<Eu_1_x..Eu_4_x>', 'amplitude': '<tag>'},
                            'returns': 'VASP_calculations/Induced_magnetic_moments/<mode>/Amplitude_<tag>: POSCAR, final '
                                       'total magnetization line, last per-ion tables, run status (INCAR and KPOINTS of '
                                       'these runs are included in every result)',
                            'amplitudes': {m: labels(data['vasp'][m]) for m in MODES}},
            'phonon_gamma': {'cost': 1, 'request': {'kind': 'phonon_gamma', 'masses': '<CrO|Cr|H>'},
                             'returns': 'Gamma-point block of band_<masses>.yaml (phonopy): masses, frequencies (THz), '
                                        'eigenvectors'},
            'force_constant_eigenmodes': {'cost': 1, 'request': {'kind': 'force_constant_eigenmodes'},
                                          'returns': 'contents of Cr_ph_eigenval_pol.hdf5 and Cr_ph_eigenvec_pol.hdf5 '
                                                     '(dataset names, values and 10x3 vectors as archived)'},
            'structure': {'cost': 1, 'request': {'kind': 'structure', 'key': '<relaxation|phonon>'},
                          'returns': 'relaxation: POSCAR, CONTCAR, INCAR; phonon: POSCAR, INCAR, FORCE_SETS, bands.conf'},
            'qe_inputs': {'cost': 1, 'request': {'kind': 'qe_inputs', 'folder': '<phonons_at_gamma|mode_E_u_1..4>'},
                          'returns': 'Quantum ESPRESSO / thermo_pw input files of that folder (the archive holds no QE '
                                     'outputs)'},
        },
        'amplitude_tags': "Folder tags as archived: '0_02' is +0.02 A, '-0_005' is -0.005 A, and so on.",
        'atom_order': 'POSCAR order: Cr1..Cr4 then O1..O6 (Elk lists the same atoms per species in the same order).',
    }


def main():
    data = lab_data()
    Z = M.born_charges()
    se = {m: M.elk_slope(m, ELK_TAGS) for m in MODES}
    a_ref, c_ref = M.alpha(Z, se)
    f_ref = M.fractions(c_ref)
    frac0 = P['born']['0_0']['reference']['poscar']['frac']
    elk_dom = M.domain(frac0, [row[2] for row in P['elk']['Eu_3_x']['0_01']['muffin_tin_atom_moments_muB'][:4]])
    vasp_dom = M.domain(frac0, [row[-1] for row in P['vasp_moments']['Eu_3_x']['0_025']['ion_magnetization_z'][:4]])
    assert (elk_dom, vasp_dom) == ('out_pointing', 'in_pointing')
    assert a_ref < 0
    ref = {'alpha_out_pointing_ps_m': a_ref, 'mode_fractions': f_ref,
           'dmdu_out_pointing_muB_per_A': {m: se[m] for m in ('Eu_3_x', 'Eu_4_x')},
           'born_Cr1': {'xx': Z[0, 0, 0], 'xy': Z[0, 0, 1], 'zz': Z[0, 2, 2]},
           'dataset_domain': {'VASP': vasp_dom, 'Elk': elk_dom},
           'dominant_mode': max(f_ref, key=f_ref.get), 'annealed_domain_parallel_EH': 'out_pointing'}

    def payload(alpha, fr, dmdu, born, dom=None, dominant=None, anneal='out_pointing'):
        return {'alpha_perp_ps_per_m': {'out_pointing': r(alpha), 'in_pointing': r(-alpha)},
                'mode_fractions': {m: r(fr[m]) for m in MODES},
                'dominant_mode': dominant or max(fr, key=fr.get),
                'dM_du_out_pointing_muB_per_A': {m: r(dmdu[m]) for m in ('Eu_3_x', 'Eu_4_x')},
                'born_charge_Cr1': {k: r(v) for k, v in born.items()},
                'dataset_domain': dom or {'VASP': vasp_dom, 'Elk': elk_dom},
                'domain_selected_by_parallel_field_annealing': anneal}

    born_ref = ref['born_Cr1']
    ab, rob = {}, {}
    # ---- robustness
    sv = {m: -M.vasp_slope(m, ['0_025', '-0_025']) for m in MODES}          # VASP in-pointing -> flip
    a, c = M.alpha(Z, sv)
    rob['vasp_flipped'] = payload(a, M.fractions(c), sv, born_ref)
    sw = {m: M.elk_slope(m, ['0_01', '-0_01', '0_015', '-0_015']) for m in MODES}
    Z5 = M.born_charges(('0_05', '-0_05'))
    a, c = M.alpha(Z5, sw)
    rob['elk_window_010_015_born_005'] = payload(a, M.fractions(c), sw, {'xx': Z5[0, 0, 0], 'xy': Z5[0, 0, 1], 'zz': Z5[0, 2, 2]})
    sm = {m: 0.5 * (se[m] + sv[m]) for m in MODES}
    a, c = M.alpha(Z, sm)
    rob['elk_vasp_mean'] = payload(a, M.fractions(c), sm, born_ref)
    Z1 = M.born_charges(('0_02',))
    a, c = M.alpha(Z1, se)
    rob['born_one_sided_002'] = payload(a, M.fractions(c), se, {'xx': Z1[0, 0, 0], 'xy': Z1[0, 0, 1], 'zz': Z1[0, 2, 2]})
    # ---- ablations
    # dp1 real-phonon (mass-weighted) eigenvectors used as frozen patterns with the frozen FC-mode responses
    fc = P['fc_eigen']
    g = P['phonon_gamma']['CrO']
    mass = np.repeat(g['masses_amu'], 3)
    Ey = np.array([np.array(fc['vectors'][m.replace('_x', '_y')]).ravel() for m in MODES])
    Ex = np.array([np.array(fc['vectors'][m]).ravel() for m in MODES])
    freqs = np.array(g['frequencies_THz'])
    pairs = []
    for n in range(len(freqs)):
        u = np.array(g['eigenvectors'][n]).ravel() / np.sqrt(mass)
        u /= np.linalg.norm(u)
        tr = np.array([np.tile([1.0, 0, 0], 10), np.tile([0, 1.0, 0], 10)]) / np.sqrt(10)
        if np.linalg.norm(np.r_[Ex @ u, Ey @ u, tr @ u]) > 0.99 and freqs[n] > 1:
            pairs.append(n)
    pairs = [tuple(pairs[i:i + 2]) for i in range(0, len(pairs), 2)]
    assert len(pairs) == 4
    vecs, stiff = {}, {}
    for m, (i, j) in zip(MODES, pairs):
        ui = np.array(g['eigenvectors'][i]).ravel() / np.sqrt(mass)
        uj = np.array(g['eigenvectors'][j]).ravel() / np.sqrt(mass)
        _, _, vt = np.linalg.svd(np.c_[Ey @ ui, Ey @ uj])
        u = vt[-1][0] * ui + vt[-1][1] * uj
        nu = np.linalg.norm(u)
        uhat = u / nu
        if (Ex @ uhat)[MODES.index(m)] < 0:
            uhat = -uhat
        vecs[m] = uhat.reshape(10, 3)
        stiff[m] = (2 * np.pi * freqs[i] * 1e12) ** 2 / nu ** 2 * M.AMU_EV
    a, c = M.alpha(Z, se, stiff=stiff, vecs=vecs)
    ab['phonon_eigenvectors_as_modes'] = payload(a, M.fractions(c), se, born_ref)
    thz = {m: 1.00794 * (2 * np.pi * fc['values'][m] * 1e12) ** 2 * M.AMU_EV for m in MODES}
    a, c = M.alpha(Z, se, stiff=thz)
    ab['hdf5_values_as_THz'] = payload(a, M.fractions(c), se, born_ref)
    # dp2 Elk domain read from the input field bfcmt (opposite to the converged moments)
    bf = []
    for l in P['elk']['Eu_3_x']['0_01']['elk_in_atoms_block'].split('\n'):
        try:
            v = [float(x) for x in l.split()[:6]]
        except ValueError:
            continue
        if len(v) == 6:
            bf.append(v[5])
    bf = bf[:4]
    elk_dom_input = M.domain(frac0, bf)
    assert elk_dom_input == 'in_pointing'
    a, c = M.alpha(Z, se)    # computed for what is believed to be the in-pointing domain
    ab['elk_domain_from_input_field'] = payload(-a, M.fractions(c), {m: -se[m] for m in se}, born_ref,
                                                dom={'VASP': vasp_dom, 'Elk': 'in_pointing'})
    # dp3 the axial Cr2-Cr3 pair (4.1 A, same column) taken as the nearest-neighbour pair
    d_long = {'VASP': M.domain(frac0, [row[-1] for row in P['vasp_moments']['Eu_3_x']['0_025']['ion_magnetization_z'][:4]], pair=(1, 2)),
              'Elk': M.domain(frac0, [row[2] for row in P['elk']['Eu_3_x']['0_01']['muffin_tin_atom_moments_muB'][:4]], pair=(1, 2))}
    assert d_long == {'VASP': 'out_pointing', 'Elk': 'in_pointing'}
    ab['axial_second_pair_as_neighbours'] = payload(-a_ref, f_ref, {m: -se[m] for m in se}, born_ref, dom=d_long)
    # dp4 both codes assumed to model the same domain and averaged
    sav = {m: 0.5 * (se[m] - sv[m]) for m in MODES}
    a, c = M.alpha(Z, sav)
    ab['codes_assumed_same_domain_averaged'] = payload(a, M.fractions(c), sav, born_ref,
                                                       dom={'VASP': 'out_pointing', 'Elk': 'out_pointing'})
    # dp5 Born charges from the mislabeled +-0.005 folders (POSCAR displacement taken at face value)
    Zb = M.born_charges(('0_005', '-0_005'))
    a, c = M.alpha(Zb, se)
    ab['born_from_mislabeled_0005_pair'] = payload(a, M.fractions(c), se, {'xx': Zb[0, 0, 0], 'xy': Zb[0, 0, 1], 'zz': Zb[0, 2, 2]})
    # dp6 Born tensor transposed (displacement index first) / species-averaged tensors
    Zt = M.born_charges(transpose=True)
    a, c = M.alpha(Zt, se)
    ab['born_tensor_transposed'] = payload(a, M.fractions(c), se, {'xx': Zt[0, 0, 0], 'xy': Zt[0, 0, 1], 'zz': Zt[0, 2, 2]})
    Zs = Z.copy()
    Zs[:4] = Z[:4].mean(0)
    Zs[4:] = Z[4:].mean(0)
    a, c = M.alpha(Zs, se)
    ab['born_species_averaged'] = payload(a, M.fractions(c), se, {'xx': Zs[0, 0, 0], 'xy': Zs[0, 0, 1], 'zz': Zs[0, 2, 2]})
    # dp11 one-sided differences against Disp_0_0 without the acoustic sum rule (reference dipole offset)
    Zn = M.born_charges(('0_02',), asr=False)
    a, c = M.alpha(Zn, se)
    ab['born_one_sided_no_asr'] = payload(a, M.fractions(c), se, {'xx': Zn[0, 0, 0], 'xy': Zn[0, 0, 1], 'zz': Zn[0, 2, 2]})
    # dp7 Eu_3 slope from the +-0.03 pair, which includes the prematurely 'converged' +0.03 run
    s3 = dict(se)
    s3['Eu_3_x'] = M.elk_slope('Eu_3_x', ['0_03', '-0_03'])
    a, c = M.alpha(Z, s3)
    ab['elk_eu3_from_003_pair'] = payload(a, M.fractions(c), s3, born_ref)
    # dp8 muffin-tin (sphere-projected) moments instead of the net cell moment
    smt = {m: M.elk_slope(m, ELK_TAGS, 'total_in_muffin_tins_muB') for m in MODES}
    a, c = M.alpha(Z, smt)
    ab['muffin_tin_moments'] = payload(a, M.fractions(c), smt, born_ref)
    svi = {m: -M.vasp_slope(m, ['0_025', '-0_025'], 'ion') for m in MODES}
    a, c = M.alpha(Z, svi)
    ab['vasp_ion_projected_moments'] = payload(a, M.fractions(c), svi, born_ref)
    # dp9 conventional hexagonal cell volume (3 x rhombohedral) in the normalisation
    ab['hexagonal_cell_volume'] = payload(a_ref / 3, f_ref, se, born_ref)
    # dp10 annealing: alpha_parallel taken to keep its 0 K sign up to T_N
    ab['annealing_ignores_sign_change'] = payload(a_ref, f_ref, se, born_ref, anneal='in_pointing')
    # hidden data
    hidden = {'reference': {k: (v if not isinstance(v, (float, np.floating)) else float(v)) for k, v in ref.items()},
              'thresholds': THRESH}
    hidden['reference']['born_Cr1'] = {k: float(v) for k, v in ref['born_Cr1'].items()}
    hidden['reference']['mode_fractions'] = {k: float(v) for k, v in f_ref.items()}
    hidden['reference']['dmdu_out_pointing_muB_per_A'] = {k: float(v) for k, v in ref['dmdu_out_pointing_muB_per_A'].items()}
    hidden['reference']['alpha_out_pointing_ps_m'] = float(a_ref)
    values = [float(a_ref), float(-a_ref), *[float(v) for v in f_ref.values()],
              *[float(v) for v in ref['dmdu_out_pointing_muB_per_A'].values()], *[float(v) for v in ref['born_Cr1'].values()]]
    lab = TASK/'environment/lab'
    (lab/'private').mkdir(parents=True, exist_ok=True)
    (lab/'data.json').write_text(json.dumps(data))
    (lab/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (lab/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog(data), indent=1) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': ab, 'robustness': rob}, indent=1) + '\n')
    # calibration printout
    print('reference', json.dumps({k: v for k, v in hidden['reference'].items()}, indent=1))
    for name, p in {**{'ROB ' + k: v for k, v in rob.items()}, **{'AB ' + k: v for k, v in ab.items()}}.items():
        print(f"{name:45s} alpha={p['alpha_perp_ps_per_m']['out_pointing']:+.4f} "
              f"f3={p['mode_fractions']['Eu_3_x']:.3f} f4={p['mode_fractions']['Eu_4_x']:.3f} "
              f"d3={p['dM_du_out_pointing_muB_per_A']['Eu_3_x']:+.4f} d4={p['dM_du_out_pointing_muB_per_A']['Eu_4_x']:+.4f} "
              f"Z={p['born_charge_Cr1']} dom={p['dataset_domain']} {p['dominant_mode']} {p['domain_selected_by_parallel_field_annealing']}")


if __name__ == '__main__':
    main()
