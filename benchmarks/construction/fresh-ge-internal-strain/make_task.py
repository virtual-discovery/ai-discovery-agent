"""Derive the hidden references, ablation values and calibration for fresh-ge-internal-strain.

    python3 -B benchmarks/construction/fresh-ge-internal-strain/make_task.py --stage-a   # writes reference_lattice.json
    python3 -B benchmarks/construction/fresh-ge-internal-strain/make_task.py             # everything else

Reads reference_runs.jsonl (lab results of reference_plan.PLAN made through task_lab.run by run_reference.py)
and applies the estimators of the task's own solution/solve.py. Writes environment/lab/private/hidden.json,
private/hidden_values.json, solution/ablation_values.json and the "calibration" block of design.json, and
prints the calibration table.
"""
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-ge-internal-strain'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
ROWS = {r['tag']: r for r in map(json.loads, (HERE/'reference_runs.jsonl').read_text().splitlines()) if r}
P6 = S.P_TARGET
A_EXP = 5.658
TOL = {'absolute': {'a0_A': 0.006, 'gamma_gap_at_a0_eV': 0.05, 'a_6GPa_A': 0.006, 'kleinman_zeta': 0.012,
                    'raman_frequency_cm-1': 1.5, 'raman_gruneisen': 0.05},
       'relative': {'C11_GPa': 0.02, 'C12_GPa': 0.03, 'C44_clamped_ion_GPa': 0.015, 'C44_relaxed_ion_GPa': 0.02}}
FIELDS = ['a0_A', 'gamma_gap_at_a0_eV', 'a_6GPa_A', 'C11_GPa', 'C12_GPa', 'C44_clamped_ion_GPa',
          'C44_relaxed_ion_GPa', 'kleinman_zeta', 'raman_frequency_cm-1', 'raman_gruneisen']


def res(tag):
    row = ROWS[tag]
    if 'result' not in row:
        raise SystemExit(f'reference calculation {tag} failed: {row.get("error")}')
    return row['result']


def eos(prefix_a, suffix):
    return [res(t) for t in ROWS if t.startswith('cub_') and t.endswith('_' + suffix)
            and t.split('_')[1].replace('.', '', 1).isdigit() and abs(float(t.split('_')[1]) - prefix_a) < 0.1]


def pick(lattices, suffix):
    return [res(f'cub_{a}_{suffix}') for a in lattices]


EOS0, EOS6 = (5.75, 5.77), (5.59, 5.61)          # the oracle's bracketing points (k12s0, 60 Ry)


def stage_a():
    a0 = S.lattice_at_pressure(pick(EOS0, 'k12s0_e60'), 0.0)
    a6 = S.lattice_at_pressure(pick(EOS6, 'k12s0_e60'), P6)
    out = {'A0': round(a0, 4), 'A6': round(a6, 4), 'A0_exact': a0, 'A6_exact': a6}
    (HERE/'reference_lattice.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))
    for suffix in ('k12s0_e45', 'k12s0_e35', 'k12s0_psl30', 'k12s1_e60'):
        try:
            print(suffix, 'a0', S.lattice_at_pressure(eos(5.76, suffix), 0.0),
                  'a6', S.lattice_at_pressure(eos(5.60, suffix), P6) if eos(5.60, suffix) else None)
        except Exception as exc:
            print(suffix, exc)


def cauchy_energy_c44(e4_tag, cub_tag, pressure_work):
    """Clamped C44 from the energy of the sheared cell; optionally removing the -P dV work of the pre-stress."""
    e, c = res(e4_tag), res(cub_tag)
    a = e['settings']['lattice_A']
    e4 = 2 * e['settings']['strain'][1][2]
    om = S.omega_bohr3(a)
    dE = e['total_energy_Ry'] - c['total_energy_Ry']
    if pressure_work:     # det(1 + eps) = 1 - (e4/2)^2 for a pure shear; dE includes +P |dV|
        dE -= c['pressure_GPa'] / S.RY_BOHR3_GPA * om * (e4 / 2) ** 2
    return 2 * dE / (om * e4 ** 2) * S.RY_BOHR3_GPA


def main():
    L = json.loads((HERE/'reference_lattice.json').read_text())
    A0, A6 = L['A0'], L['A6']
    tag6 = lambda t: t                                   # noqa: E731 (tags already encode a6)
    # ---------------- reference (the oracle's estimators on the reference runs)
    ref = {}
    ref['a0_A'] = S.lattice_at_pressure(pick(EOS0, 'k12s0_e60'), 0.0)
    ref['gamma_gap_at_a0_eV'] = S.interpolate([(a, S.gamma_gap(res(f'cub_{a}_k12s0_e60'))) for a in EOS0], ref['a0_A'])
    ref['a_6GPa_A'] = S.lattice_at_pressure(pick(EOS6, 'k12s0_e60'), P6)
    ref['C11_GPa'], ref['C12_GPa'] = S.normal_constants(res('exx_0.01_a6_k12s0'), res('exx_-0.01_a6_k12s0'))
    I = S.internal_strain(res('e4_a6_k12s0'), res('ph100_a6_k12s0'))
    ref['C44_clamped_ion_GPa'], ref['C44_relaxed_ion_GPa'], ref['kleinman_zeta'] = I['C44_clamped'], I['C44_relaxed'], I['zeta']
    ref['raman_frequency_cm-1'] = S.frequency_cm(I['Phi_Ry_per_bohr2'])
    ref['raman_gruneisen'] = S.gruneisen(res('ph100_a6mh_k12s0'), res('ph100_a6ph_k12s0'))

    def fields_from(e4, ph, extra=None):
        X = S.internal_strain(res(e4), res(ph))
        out = {'C44_clamped_ion_GPa': X['C44_clamped'], 'C44_relaxed_ion_GPa': X['C44_relaxed'],
               'kleinman_zeta': X['zeta'], 'raman_frequency_cm-1': S.frequency_cm(X['Phi_Ry_per_bohr2'])}
        out.update(extra or {})
        return out

    # ---------------- robustness: defensible alternatives that must pass
    Rb = {}
    Rb['gamma_k16'] = fields_from('e4_a6_k16s0', 'ph100_a6_k16s0')
    Rb['cutoff_60_Ry'] = fields_from('e4_a6_k12s0_e60', 'ph100_a6_k12s0_e60')
    Rb['small_smearing_mv002'] = fields_from('e4_a6_k12s0_mv002', 'ph100_a6_k12s0_mv002')
    scan = S.relaxed_from_scan(res('e4_a6_k12s0'), res('e4u_a6_k12s0'))
    Rb['displacement_scan'] = {'C44_relaxed_ion_GPa': scan['C44_relaxed'], 'kleinman_zeta': scan['zeta']}
    rel = res('relax_e4_a6_k12s0_1e-5')
    cell = rel['cell_vectors_A']
    pos = [[sum(f[i] * cell[i][j] for i in range(3)) for j in range(3)] for f in rel['final_positions_crystal']]
    b = [sum(c[j] for c in cell) * 0.25 for j in range(3)]           # (1+eps)(a/4)(1,1,1) = (a1+a2+a3)/4
    ux = pos[1][0] - pos[0][0] - b[0]
    e4 = 2 * rel['settings']['strain'][1][2]
    Rb['relaxation_1e-5'] = {'C44_relaxed_ion_GPa': S.cauchy(rel)[1][2] / e4, 'kleinman_zeta': -4 * ux / (A6 * e4)}
    # Lambda from the stress of the displaced cell (Maxwell relation) instead of the force of the sheared cell
    ph = res('ph100_a6_k12s0')
    lam_m = S.cauchy(ph)[1][2] / S.RY_BOHR3_GPA * S.omega_bohr3(A6) / (ph['settings']['displacement_A'][0] / S.BOHR_A)
    phi = I['Phi_Ry_per_bohr2']
    Rb['maxwell_lambda'] = {'C44_relaxed_ion_GPa': I['C44_clamped'] - lam_m ** 2 / (S.omega_bohr3(A6) * phi) * S.RY_BOHR3_GPA,
                            'kleinman_zeta': 4 * lam_m / (A6 / S.BOHR_A * phi)}
    phi111 = (res('ph111_-0.02_a6_k12s0')['forces_Ry_per_bohr'][1][0] - res('ph111_0.02_a6_k12s0')['forces_Ry_per_bohr'][1][0]) / (2 * 0.02 / S.BOHR_A)
    Rb['central_111_displacement'] = {'raman_frequency_cm-1': S.frequency_cm(phi111)}
    phb = res('ph100b_a6_k12s0')
    Rb['displacement_0.03_along_x'] = {'raman_frequency_cm-1': S.frequency_cm(-phb['forces_Ry_per_bohr'][1][0] / (0.03 / S.BOHR_A))}
    Rb['normal_strain_0.005'] = dict(zip(('C11_GPa', 'C12_GPa'), S.normal_constants(res('exx_0.005_a6_k12s0'), res('exx_-0.005_a6_k12s0'))))
    Rb['normal_strain_shifted_k8'] = dict(zip(('C11_GPa', 'C12_GPa'), S.normal_constants(res('exx_0.01_a6_k8s1'), res('exx_-0.01_a6_k8s1'))))
    Rb['gruneisen_compression_side_1pct'] = {'raman_gruneisen': S.gruneisen(res('ph100_a6m_k12s0'), res('ph100_a6_k12s0'))}
    Rb['gruneisen_compression_side_half_pct'] = {'raman_gruneisen': S.gruneisen(res('ph100_a6mh_k12s0'), res('ph100_a6_k12s0'))}
    Rb['gruneisen_expansion_side_half_pct'] = {'raman_gruneisen': S.gruneisen(res('ph100_a6_k12s0'), res('ph100_a6ph_k12s0'))}
    Rb['pressures_at_45_Ry'] = {'a0_A': S.lattice_at_pressure(eos(5.76, 'k12s0_e45'), 0.0),
                                'a_6GPa_A': S.lattice_at_pressure(eos(5.60, 'k12s0_e45'), P6)}
    Rb['five_point_pressure_fits'] = {'a0_A': S.lattice_at_pressure(eos(5.76, 'k12s0_e60'), 0.0),
                                      'a_6GPa_A': S.lattice_at_pressure(eos(5.60, 'k12s0_e60'), P6)}
    Rb['gap_direct_at_a0_k8'] = {'gamma_gap_at_a0_eV': S.gamma_gap(res('cub_a0_k8s0_e60'))}
    Rb['shifted_mesh_pressures'] = {'a0_A': S.lattice_at_pressure(eos(5.76, 'k12s1_e60') + [res('cub_5.76_k12s0_e60')], 0.0)}
    Rb['gap_smeared_k12'] = {'gamma_gap_at_a0_eV': S.gamma_gap(res('cub_a0_k12s0_e60_mv002'))}
    Rb['energy_route_with_pressure_work'] = {'C44_clamped_ion_GPa': cauchy_energy_c44('e4_a6_k12s0', 'cub_a6_k12s0', True)}

    # ---------------- ablations: realistic wrong paths computed from the reference runs
    A = {}
    sh = res('cub_a0_k8s1_e45')
    A['gap_from_shifted_mesh_band_edges'] = {'gamma_gap_at_a0_eV': sh['lowest_unoccupied_level_eV'] - sh['highest_occupied_level_eV']}
    g = res('cub_a0_k8s0_e60')
    A['gap_from_gamma_mesh_band_edges'] = {'gamma_gap_at_a0_eV': g['lowest_unoccupied_level_eV'] - g['highest_occupied_level_eV']}
    A['pressures_at_35_Ry'] = {'a0_A': S.lattice_at_pressure(eos(5.76, 'k12s0_e35'), 0.0),
                               'a_6GPa_A': S.lattice_at_pressure(eos(5.60, 'k12s0_e35'), P6)}
    A['frozen_3d_core_pseudopotential'] = dict(fields_from('e4_a6_k12s0_psl30', 'ph100_a6_k12s0_psl30'),
                                               a0_A=S.lattice_at_pressure(eos(5.76, 'k12s0_psl30'), 0.0),
                                               a_6GPa_A=S.lattice_at_pressure(eos(5.60, 'k12s0_psl30'), P6))
    A['shifted_k8_mesh'] = fields_from('e4_a6_k8s1', 'ph100_a6_k8s1')
    INFO = {'shifted_k12_mesh (not graded as robustness or ablation)': fields_from('e4_a6_k12s1', 'ph100_a6_k12s1')}
    rs = res('relax_e4_a6_k8s1_1e-4')
    cell = rs['cell_vectors_A']
    pos = [[sum(f[i] * cell[i][j] for i in range(3)) for j in range(3)] for f in rs['final_positions_crystal']]
    b = [sum(c[j] for c in cell) * 0.25 for j in range(3)]
    A['shifted_k8_relaxation'] = {'C44_relaxed_ion_GPa': S.cauchy(rs)[1][2] / 0.01,
                                  'kleinman_zeta': -4 * (pos[1][0] - pos[0][0] - b[0]) / (A6 * 0.01)}
    A['gruneisen_from_1pct_steps'] = {'raman_gruneisen': S.gruneisen(res('ph100_a6m_k12s0'), res('ph100_a6p_k12s0'))}
    A['gruneisen_from_1pct_steps_k8'] = {'raman_gruneisen': S.gruneisen(res('ph100_a6m_k8s0'), res('ph100_a6p_k8s0'))}
    A['wide_smearing_mv002'] = dict(fields_from('e4_a6_k8s1_mv02', 'ph100_a6_k8s1_mv02'),
                                    raman_gruneisen=S.gruneisen(res('ph100_a6m_k8s1_mv02'), res('ph100_a6p_k8s1_mv02')))
    A['smearing_0.01_gamma_k12'] = {'raman_frequency_cm-1': S.frequency_cm(-res('ph100_a6_k12s0_mv01')['forces_Ry_per_bohr'][1][0] / (0.01 / S.BOHR_A))}
    A['clamped_reported_as_relaxed'] = {'C44_relaxed_ion_GPa': I['C44_clamped'], 'kleinman_zeta': 0.0}
    loose = res('relax_e4s_a6_k12s0_1e-3')
    A['loose_relaxation_at_small_strain'] = {'C44_relaxed_ion_GPa': S.cauchy(loose)[1][2] / 0.002,
                                             'kleinman_zeta': 0.0 if len(loose['ionic_steps']) == 1 else None}
    for d in (0.02, 0.05, -0.05):
        r = res(f'ph111_{d}_a6_k12s0')
        A[f'one_sided_111_displacement_{d}'] = {'raman_frequency_cm-1': S.frequency_cm(-r['forces_Ry_per_bohr'][1][0] / (d / S.BOHR_A))}
    A['energy_route_ignoring_prestress'] = {'C44_clamped_ion_GPa': cauchy_energy_c44('e4_a6_k12s0', 'cub_a6_k12s0', False)}
    A['atomic_mass_in_frequency'] = {'raman_frequency_cm-1': ref['raman_frequency_cm-1'] / math.sqrt(2)}
    A['tensor_shear_strain_as_engineering'] = {'C44_clamped_ion_GPa': 2 * ref['C44_clamped_ion_GPa'],
                                               'C44_relaxed_ion_GPa': 2 * ref['C44_relaxed_ion_GPa'],
                                               'kleinman_zeta': 2 * ref['kleinman_zeta']}
    one = res('exx_0.02_a6_k12s0')
    c0 = S.cauchy(res('cub_a6_k12s0'))
    A['one_sided_normal_strain_0.02'] = {'C11_GPa': (S.cauchy(one)[0][0] - c0[0][0]) / 0.02,
                                         'C12_GPa': (S.cauchy(one)[1][1] - c0[1][1]) / 0.02}
    A['experimental_lattice_parameter'] = {'a0_A': A_EXP}
    for name in list(A):
        A[name] = {k: v for k, v in A[name].items() if v is not None}

    hidden = {'reference': ref, 'tolerances': TOL, 'provenance': {
        'runs': 'benchmarks/construction/fresh-ge-internal-strain/reference_runs.jsonl (task_lab.run, QE 6.7, 2 MPI ranks)',
        'lattice_parameters_of_stage_b': {'A0': A0, 'A6': A6},
        'Lambda_Ry_per_bohr': I['Lambda_Ry_per_bohr'], 'Phi_Ry_per_bohr2': I['Phi_Ry_per_bohr2']}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [ref[f] for f in FIELDS] + [I['Lambda_Ry_per_bohr'], I['Phi_Ry_per_bohr2']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': Rb}, indent=1) + '\n')

    # ---------------- calibration table and checks
    calibration, problems = {}, []
    for f in FIELDS:
        rel = f in TOL['relative']
        tol = TOL['relative'][f] if rel else TOL['absolute'][f]
        err = (lambda v: v / ref[f] - 1) if rel else (lambda v: v - ref[f])
        rob = {n: v[f] for n, v in Rb.items() if f in v}
        abl = {n: v[f] for n, v in A.items() if f in v}
        calibration[f] = {'reference': ref[f], 'tolerance': f'within {tol:.1%}' if rel else f'+-{tol}',
                          'robustness': rob, 'ablations': abl,
                          'max_robust_error': max([abs(err(v)) for v in rob.values()], default=None),
                          'nearest_ablation_error': min([abs(err(v)) for v in abl.values()], default=None)}
        problems += [f'robustness {n} fails {f}: {v}' for n, v in rob.items() if abs(err(v)) > tol]
    for n, v in A.items():
        if all((abs(v[f] / ref[f] - 1) <= TOL['relative'][f]) if f in TOL['relative'] else (abs(v[f] - ref[f]) <= TOL['absolute'][f]) for f in v):
            problems.append(f'ablation {n} passes every field it changes')
    design_path = TASK/'design.json'
    design = json.loads(design_path.read_text())
    design['calibration'] = json.loads(json.dumps(calibration, default=float))
    design_path.write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps(calibration, indent=1, default=float))
    design['calibration']['informational'] = json.loads(json.dumps(INFO, default=float))
    design_path.write_text(json.dumps(design, indent=1) + '\n')
    print('INFO', json.dumps(INFO, indent=1, default=float))
    print('PROBLEMS:', problems or 'none')


if __name__ == '__main__':
    stage_a() if '--stage-a' in sys.argv else main()
