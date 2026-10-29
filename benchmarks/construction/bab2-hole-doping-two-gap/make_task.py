"""Generate the staged, lab, private and ablation data for bab2-hole-doping-two-gap from pool.json.

    bash benchmarks/construction/bab2-hole-doping-two-gap/run_in_image.sh make_task.py

Writes catalog.json, lab/data.json, private/hidden.json, private/hidden_values.json,
solution/ablation_values.json, the budget in both public.json files, and calibration.json (this folder).
Asserts the archive facts the task relies on and prints the calibration numbers.
"""
import json
import re
import numpy as np
import model as M

TASK, HERE, RUNS, SUPER = M.TASK, M.HERE, M.RUNS, M.SUPER
EL = M.EL
BUDGET = 34
COSTS = {'dfpt_phonons': 2, 'a2f': 1, 'iso_log': 1, 'aniso_gap': 2, 'eliashberg_iso': 1, 'supercell_phonons': 2}
MU, CUTOFF = 0.10, 500.0
THRESHOLDS = {'critical_doping_abs': 0.0025, 'iso_tc_rel': 0.03, 'aniso_tc_rel': 0.04, 'gap_ratio_rel': 0.10,
              'homogeneous_soft_mode_abs_cm1': 3.5, 'cs_soft_mode_abs_cm1': 1.5}


def header(folder):
    r = RUNS[folder]
    return {'run': f'sensitivity/{folder}', 'scf.in': r['scf_in'], 'scf.out (selected lines)': r['scf_out_lines'],
            'relax.out (final structure)': r['relaxed_structure']}


def lab_data():
    runs = {}
    for f, r in RUNS.items():
        rows = np.array(r['a2f']['rows'])
        smear = [float(v) for v in r['a2f']['footer'][r['a2f']['footer'].index(' Phonon smearing (meV)') + 1].split()[1:]]
        runs[f] = {
            'header': header(f),
            'dfpt': None if f in M.WITHHELD_DFPT else {
                'ph.in': r['ph_in'], 'bab2.dyn0': r['dfpt_grid'], 'bab2.dyn1 header': r['dfpt_header'],
                'q_points': r['dfpt'],
                'note': 'q in cartesian coordinates, units 2 pi/alat (alat = celldm(1) of the dyn header); '
                        'frequencies in cm-1 from the archived dynamical-matrix files, negative = imaginary'},
            'a2f': {'file': 'epw/bab2.a2f', 'first_line': r['a2f']['header'], 'rows': r['a2f']['rows'],
                    'footer': r['a2f']['footer']},
            'iso_log': {'iso.in': r['iso_in'], 'iso.out (isotropic Eliashberg section)': r['iso_log']},
            'aniso_in': r['aniso_in'], 'aniso_log_head': r['aniso_log_head'],
            'aniso': {T: {'imaginary': a['imaginary'], 'real': a['real'], 'log': a['log']} for T, a in r['aniso'].items()},
            'solver_a2f': {'w_meV': rows[:, 0].tolist(), 'columns': rows[:, 1:11].T.tolist(), 'smearing_meV': smear}}
    supercells = {}
    for name, s in SUPER.items():
        entry = {'run': name, 'scf input': s['scf_in'], 'scf.out (selected lines)': s['scf_out_lines'], 'ph.in': s['ph_in'],
                 'bab2.dyn0': s['dfpt_grid'], 'bab2.dyn1 header': s['dfpt_header'], 'q_points': s['dfpt'],
                 'note': 'q in cartesian coordinates, units 2 pi/alat of the supercell; frequencies in cm-1, negative = imaginary'}
        if 'relax_in' in s:
            entry['relax.in'] = s['relax_in']
            entry['relax.out (final structure)'] = s['relaxed_structure']
        supercells[name] = entry
    return {'costs': COSTS, 'runs': runs, 'supercells': supercells, 'withheld_dfpt': [f'sensitivity/{f}' for f in M.WITHHELD_DFPT]}


def catalog():
    runs = []
    for f in '0123456':
        runs.append({'run': f'sensitivity/{f}', 'archive_folder': f'bab2-analysis/sensitivity/{f}',
                     'purchasable': (['iso_log', 'a2f', 'aniso_gap', 'eliashberg_iso'] +
                                     ([] if f in M.WITHHELD_DFPT else ['dfpt_phonons'])),
                     'archived_anisotropic_temperatures_K': sorted(float(T) for T in RUNS[f]['aniso'])})
    runs += [{'run': '331', 'archive_folder': 'bab2-analysis/331', 'purchasable': ['supercell_phonons']},
             {'run': '331_cs', 'archive_folder': 'bab2-analysis/331_cs', 'purchasable': ['supercell_phonons']}]
    return {
        'title': 'Hole-doped BaB2: archived DFPT, EPW and Migdal-Eliashberg results (Bercx et al., PRX Energy 4, 033012 (2025))',
        'material': 'BaB2 in the AlB2 structure (P6/mmm): one formula unit (Ba at 1a, B at 2d) per primitive cell; hexagonal B layers as in MgB2.',
        'archive_readme_excerpt': POOL_README,
        'methods': ('Quantum ESPRESSO, PBE, 80 Ry cutoff, cold smearing 0.02 Ry, 12x12x8 k-points; holes added as a homogeneous '
                    'background charge with the cell relaxed at each charge; DFPT on a 6x6x4 q-grid (21 irreducible '
                    'points). EPW v5.8 with Wannier interpolation to 40x40x24 k and q grids; isotropic and anisotropic Migdal-Eliashberg '
                    'equations on the imaginary axis, continued to the real axis with Pade approximants. Supercells: 3x3x1 '
                    'cells (27 atoms), phonons at two wave vectors. The returned input files and log lines record each run\'s settings.'),
        'runs': runs,
        'requests': {
            'dfpt_phonons': {'example': {'kind': 'dfpt_phonons', 'run': 'sensitivity/0'}, 'cost': COSTS['dfpt_phonons'],
                             'returns': 'harmonic frequencies at the 21 irreducible DFPT q-points (archived bab2.dyn1-21) with ph.in and the run header'},
            'a2f': {'example': {'kind': 'a2f', 'run': 'sensitivity/0'}, 'cost': COSTS['a2f'],
                    'returns': 'the archived EPW file bab2.a2f (rows as archived, with its first line and footer) and the run header'},
            'iso_log': {'example': {'kind': 'iso_log', 'run': 'sensitivity/0'}, 'cost': COSTS['iso_log'],
                        'returns': "the archived isotropic Eliashberg input (iso.in) and the solver section of iso.out"},
            'aniso_gap': {'example': {'kind': 'aniso_gap', 'run': 'sensitivity/0', 'temperature_K': 15, 'axis': 'imaginary'},
                          'cost': COSTS['aniso_gap'],
                          'returns': ("the archived EPW gap-distribution file at one archived temperature, on the imaginary axis "
                                      "(bab2.imag_aniso_gap0_T) or continued to the real axis (bab2.pade_aniso_gap0_T), as archived "
                                      "(first line and two columns), with that temperature's log lines, epw_int.in and the run header")},
            'eliashberg_iso': {'example': {'kind': 'eliashberg_iso', 'run': 'sensitivity/0', 'mu_star': 0.13, 'cutoff_meV': 500,
                                           'temperatures_K': [15, 20, 25]},
                               'optional': {'a2f_column': '1-10 (smearing column of bab2.a2f; default 1)',
                                            'max_iter': 'default 500, at most 5000', 'conv_thr': 'default 1e-3, at least 1e-10',
                                            'mixing': 'default 0.4, in (0, 1]', 'initial_gap_meV': 'default 5'},
                               'cost': COSTS['eliashberg_iso'],
                               'returns': 'a new numpy solution of the isotropic Migdal-Eliashberg equations on the archived alpha2F of that run, at up to 8 temperatures (5-150 K)'},
            'supercell_phonons': {'example': {'kind': 'supercell_phonons', 'run': '331'}, 'cost': COSTS['supercell_phonons'],
                                  'returns': 'frequencies at the two archived supercell wave vectors with the scf/ph inputs, selected scf output lines (and the relaxation for 331_cs)'}},
        'withheld': ['dfpt_phonons of sensitivity/3', 'dfpt_phonons of sensitivity/5'],
        'eliashberg_iso_solver': (
            'Matsubara frequencies w_n = (2n+1) pi kB T, n = 0..N-1, N = int((cutoff/(pi kB T) - 1)/2) + 1; with '
            'lambda(m) = int 2 w a2F(w)/(w^2 + (2 m pi kB T)^2) dw: '
            'Z_n = 1 + (pi kB T/w_n) sum_m [lambda(n-m) - lambda(n+m+1)] w_m/sqrt(w_m^2 + D_m^2); '
            'Z_n D_n = pi kB T sum_m [lambda(n-m) + lambda(n+m+1) - 2 mu*] D_m/sqrt(w_m^2 + D_m^2) (the form solved by EPW). '
            'Linear mixing D <- mixing*D_new + (1-mixing)*D_old from a constant initial gap; residual = sum|D_new - D_old|/sum|D_new|; '
            'a run stops when the residual falls below conv_thr, when D(iw_0) falls below 1e-6 meV (reported as zero), or at max_iter. '
            'Returns D(iw_0) in meV, Z(iw_0), N, the iteration count, the final residual and whether it converged.'),
        'purchase_costs': COSTS, 'budget': BUDGET}


POOL_README = M.POOL['readme_bab2']


def tc_from_scan(folder, temps, mu=MU, **kw):
    """Tc as the first scanned temperature where the reported gap is below 0.01 meV (unconverged values taken at face value)."""
    w, a = M.a2f_arrays(folder)
    for T in temps:
        r = EL.solve(w, a, T, mu, CUTOFF, **kw)
        if r['gap_meV'] < 0.01:
            return float(T)
    return float('nan')


def tc_converged_extrapolation(folder, temps, mu=MU):
    """Oracle-style estimate: tight-tolerance runs, Delta^2 extrapolated from the two highest converged gapped temperatures."""
    w, a = M.a2f_arrays(folder)
    pts = []
    for T in temps:
        r = EL.solve(w, a, T, mu, CUTOFF, max_iter=5000, conv_thr=1e-6)
        if r['converged'] and r['gap_meV'] > 0:
            pts.append((T, r['gap_meV']))
    pts.sort()
    return M.zero_of_square([p[0] for p in pts], [p[1] for p in pts]), pts


def main():
    # --- archive facts the task relies on
    dop = {f: M.doping(f) for f in RUNS}
    assert dop == {'0': 0.02, '1': 0.04, '2': 0.06, '3': 0.08, '4': 0.1, '5': 0.065, '6': 0.0}, dop
    soft = {M.doping(f): M.soft_mode(RUNS[f]) for f in RUNS}
    assert all(soft[x][0] < 0 for x in (0.0, 0.02, 0.04, 0.06)) and all(soft[x][0] > 0 for x in (0.065, 0.08, 0.1))
    s331, scs = SUPER['331'], SUPER['331_cs']
    assert any('144.00' in l for l in s331['scf_out_lines']) and 'tot_charge' not in s331['scf_in']
    assert any('143.00' in l for l in scs['scf_out_lines']) and 'Cs' in scs['scf_in']
    assert any('11.3085' in l for l in s331['scf_out_lines']) and any('11.3085' in l for l in RUNS['6']['scf_out_lines'])
    a331 = sorted(s331['dfpt'][1]['frequencies_cm1'])[0]
    acs = sorted(scs['dfpt'][1]['frequencies_cm1'])
    gcs = sorted(scs['dfpt'][0]['frequencies_cm1'])
    assert abs(a331 - soft[0.0][0]) < 0.2            # the '331' supercell reproduces the undoped primitive cell's H mode
    iso8 = M.iso_archived('3')
    assert np.isnan(iso8[-1][1]) and not iso8[-1][2]  # archived x=0.08 isotropic run: NaN, not converged at 30 K
    # the lab solver, converged tightly at the archived settings (mu* 0.13, 0.5 eV cutoff), reproduces the archived EPW gaps
    dev = []
    for f in RUNS:
        w, a = M.a2f_arrays(f)
        for t, g, ok in M.iso_archived(f):
            if ok:
                dev.append(abs(EL.solve(w, a, t, 0.13, CUTOFF, max_iter=5000, conv_thr=1e-8)['gap_meV'] - g))
    assert max(dev) < 0.03, max(dev)
    dos = {M.doping(f): float(re.search(r'DOS\(states/spin/eV/Unit Cell\) =\s*([\d.E+-]+)', RUNS[f]['iso_log']).group(1)) for f in RUNS}

    # --- references
    xs = sorted(soft)
    w2 = [M.signed_square(soft[x][0]) for x in xs]
    xc_ref = M.crossing(xs, w2)
    rel = [x for x in xs if M.folder_of(x) not in M.WITHHELD_DFPT]
    w2rel = [M.signed_square(soft[x][0]) for x in rel]
    hom_ref = float(np.sqrt(w2[-1] + (1 / 9 - xs[-1]) * (w2[-1] - w2[-2]) / (xs[-1] - xs[-2])))
    cs_ref = float(acs[0])
    iso_ref = {k: EL.tc_linearized(*M.a2f_arrays(M.folder_of(x)), MU, CUTOFF) for k, x in (('0.065', 0.065), ('0.10', 0.1))}
    an = M.aniso_series('5')
    T = [r['T'] for r in an]
    aniso_ref = M.zero_of_square(T, [r['dmax'] for r in an])
    a25 = RUNS['5']['aniso']['25']
    small, large = M.clusters(a25['imaginary']['rows'], 25.0)
    ratio_ref = large / small

    # --- robustness (defensible alternatives)
    rob = {}
    rob['xc_w2_linear_released_0.06_0.10'] = M.crossing(rel, w2rel)
    rob['xc_w2_quadratic_released'] = M.poly_root(rel, w2rel, 2)[0]
    rob['xc_w2_linear_regression_released'] = M.poly_root(rel, w2rel, 1)[0]
    p2 = M.poly_root(rel, w2rel, 2)[1]
    rob['hom_w2_quadratic_released'] = float(np.sqrt(np.polyval(p2, 1 / 9)))
    rob['hom_w2_linear_released_0.06_0.10'] = float(np.sqrt(w2rel[-1] + (1 / 9 - 0.1) * (w2rel[-1] - w2rel[-2]) / 0.04))
    rob['hom_w2_linear_regression_released'] = float(np.sqrt(np.polyval(M.poly_root(rel, w2rel, 1)[1], 1 / 9)))
    rob['iso_converged_extrapolation'] = {}
    for k, f, temps in (('0.065', '5', [33, 34, 35, 35.5, 36, 36.4]), ('0.10', '4', [30, 31, 32, 32.5, 33, 33.3])):
        rob['iso_converged_extrapolation'][k] = tc_converged_extrapolation(f, temps)[0]
    rob['iso_a2f_column_10'] = {k: EL.tc_linearized(*M.a2f_arrays(M.folder_of(x), 10), MU, CUTOFF)
                                for k, x in (('0.065', 0.065), ('0.10', 0.1))}
    rob['aniso_large_cluster_squared'] = M.zero_of_square(T, [r['large'] for r in an])
    rob['aniso_sqrt_free_energy'] = M.zero_linear(T, list(np.sqrt(-np.array([r['F'] for r in an]))))
    p = np.polyfit(T[-3:], np.array([r['dmax'] for r in an][-3:]) ** 2, 1)
    rob['aniso_three_point_fit'] = float(-p[1] / p[0])
    rob['aniso_real_axis_max_squared'] = M.zero_of_square(T, [r['dmax'] for r in M.aniso_series('5', 'real')])
    med = M.clusters(a25['imaginary']['rows'], 25.0, 'median')
    unw = M.clusters(a25['imaginary']['rows'], 25.0, 'unweighted')
    rob['ratio_median'] = med[1] / med[0]
    rob['ratio_unweighted_bins'] = unw[1] / unw[0]

    # --- ablations (plausible wrong paths, computed from the archive)
    abl = {}
    lab = {M.README_ORDER[int(f)]: M.soft_mode(RUNS[f])[0] for f in RUNS if f not in M.WITHHELD_DFPT}
    lx = sorted(lab)
    abl['xc_readme_order'] = M.crossing(lx, [M.signed_square(lab[x]) for x in lx])
    readme_folder = {k: str(M.README_ORDER.index(x)) for k, x in (('0.065', 0.065), ('0.10', 0.1))}
    abl['iso_readme_order'] = {k: EL.tc_linearized(*M.a2f_arrays(f), MU, CUTOFF) for k, f in readme_folder.items()}
    abl['xc_omega_linear'] = M.crossing(rel, [soft[x][0] for x in rel])
    abl['hom_omega_linear'] = float(soft[0.1][0] + (1 / 9 - 0.1) * (soft[0.1][0] - soft[0.06][0]) / 0.04)
    abl['hom_331_as_homogeneous'] = float(a331)
    abl['cs_gamma_point'] = float(gcs[0])
    abl['iso_archived_mu013'] = {}
    for k, f in (('0.065', '5'), ('0.10', '4')):
        pts = [(t, g) for t, g, ok in M.iso_archived(f) if ok and g > 0]
        abl['iso_archived_mu013'][k] = M.zero_of_square([q[0] for q in pts], [q[1] for q in pts])
    abl['iso_allen_dynes'] = {k: M.allen_dynes(*M.a2f_arrays(M.folder_of(x)), MU)[0] for k, x in (('0.065', 0.065), ('0.10', 0.1))}
    abl['iso_unconverged_scan'] = {k: tc_from_scan(M.folder_of(x), range(25, 46)) for k, x in (('0.065', 0.065), ('0.10', 0.1))}
    abl['aniso_iso_reported'] = iso_ref['0.065']
    abl['aniso_linear_gap'] = M.zero_linear(T, [r['dmax'] for r in an])
    ps, pl = M.clusters(a25['real']['rows'], 25.0)
    abl['ratio_real_axis'] = pl / ps
    mode = M.clusters(a25['imaginary']['rows'], 25.0, 'mode')
    info = {'ratio_mode_convention_violation': mode[1] / mode[0],
            'aniso_last_archived_T': 60.0, 'paper_values': {'xc': 0.065, 'tc_iso_mu013': 32.0, 'tc_aniso': 61.6, 'gap_ratio': 11.54 / 1.50},
            'iso_mu013_linearized': {k: EL.tc_linearized(*M.a2f_arrays(M.folder_of(x)), 0.13, CUTOFF) for k, x in (('0.065', 0.065), ('0.10', 0.1))},
            'soft_modes': {str(x): soft[x][0] for x in xs}, 'supercell_A_331': a331, 'supercell_A_cs': acs[:2],
            'archived_iso_x008': iso8, 'solver_vs_archived_max_dev_meV': max(dev),
            'solver_vs_archived_median_dev_meV': float(np.median(dev)), 'dos_EF': dos,
            'lambda': {M.doping(f): EL.coupling(*M.a2f_arrays(f)) for f in RUNS}}

    hidden = {'targets': {'critical_hole_doping': xc_ref, 'iso_tc_K': iso_ref, 'aniso_tc_K': aniso_ref,
                          'gap_ratio_25K': ratio_ref,
                          'soft_mode_cm1': {'homogeneous_one_ninth': hom_ref, 'cs_substituted': cs_ref}},
              'thresholds': THRESHOLDS}
    values = [xc_ref, *iso_ref.values(), aniso_ref, ratio_ref, hom_ref, cs_ref]

    # --- write files
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog(), indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab_data(), separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    for path in (TASK/'environment/lab/public.json', TASK/'tests/public.json'):
        pub = json.loads(path.read_text())
        pub['solver_budget'] = BUDGET
        path.write_text(json.dumps(pub, indent=1) + '\n')
    calib = {'references': hidden['targets'], 'robustness': rob, 'ablations': abl, 'info': info, 'thresholds': THRESHOLDS}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    print(json.dumps(calib, indent=1, default=float))


if __name__ == '__main__':
    main()
