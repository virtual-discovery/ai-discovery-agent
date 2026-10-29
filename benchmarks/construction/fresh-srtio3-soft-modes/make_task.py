"""Write the hidden references, ablation/robustness values and calibration of fresh-srtio3-soft-modes.

    .venv/bin/python -B benchmarks/construction/fresh-srtio3-soft-modes/make_task.py

Reads reference_runs.jsonl (task_lab.run results), applies the estimators of solution/solve.py and writes
environment/lab/private/hidden.json, private/hidden_values.json, solution/ablation_values.json, calibration.json
and design.json -> calibration (every variation graded with task_lab.grade).
"""
import copy
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-srtio3-soft-modes'
sys.path.insert(0, str(HERE))
import reference_plan as P  # noqa: E402

S = P.S
A0 = P.A0
R = {t: r['result'] for t, r in P.rows().items()}


def fp(a, k, u=0.01):
    return {n: R[P.fp_tag(a, k, u, n)] for n in S.COLS_Z}


def low(K):
    return float(S.modes(K, S.Z_MASSES)[0][0])


def lo_max(K, d):
    return float(S.modes(K, S.Z_MASSES, S.z_nac(d))[0][-1])


def polar_depth(pat, k, amps=S.POLAR_TI_AMPS):
    c = R[f'cub_a{A0}_k{k}']
    de = [(R[f'polar_pat{pat}_k{k}_t{t}']['total_energy_Ry'] - c['total_energy_Ry']) * S.RY_EV * 1000 for t in amps]
    return S.well(amps, de)[0]


def tilt(angles, ref_tag=f'tilt_a{A0}_k4_t0'):
    ref = R[ref_tag]
    return S.tilt_analysis(ref, {t: R[f'tilt_a{A0}_k4_t{t}'] for t in angles}, A0)


def strained_low(c, k):
    return float(S.modes(S.x_block({n: R[f'strainfp_c{c}_k{k}_{n}x'] for n in S.COLS_X}), S.X_MASSES)[0][0])


def stress_zero(tags):
    cs = [R[t] for t in tags]
    return float(np.interp(0, [-r['stress_GPa'][2][2] for r in cs], [r['settings']['c_A'] for r in cs]))


def uncoupled_lo(K, d):
    """Wrong path: LO of the highest TO mode alone, from its own mode effective charge (no mixing of modes)."""
    w, v = S.modes(K, S.Z_MASSES)
    zq, eps, vol = S.z_nac(d)
    e = v[:, -1] / np.sqrt(S.Z_MASSES)
    zm = float(np.dot(zq, e))
    add = 4 * np.pi * 2.0 / (vol / S.BOHR_A ** 3) * zm ** 2 / eps
    return float(math.sqrt((w[-1] / S.RY_CM) ** 2 + add) * S.RY_CM)


def main():
    eos60 = [R[f'eos_a{a}_e60_k6'] for a in S.EOS_A]
    eos75 = [R[f'eos_a{a}_e75_k6'] for a in S.EOS_A]
    a0_75, a0_60 = S.a0_from_eos(eos75), S.a0_from_eos(eos60)
    assert round(a0_60, 3) == A0, a0_60
    K4, K6, K8 = (S.z_block(fp(A0, k)) for k in (4, 6, 8))
    K6u4 = S.z_block(fp(A0, 6, 0.04), 0.04)
    Kexp = S.z_block(fp(P.A_EXP, 6))
    d4, d6, d8 = (R[f'diel_a{A0}_k{k}'] for k in (4, 6, 8))
    tilt_all = tilt((1, 2, 3, 3.5, 4, 4.5, 5, 6))
    tilt_or = tilt(S.TILT_DEG)
    c75 = S.c_relaxed([R[f'strain_c{c}_k6_e75'] for c in S.C_SCAN])
    c60 = S.c_relaxed([R[f'strain_c{c}_k6'] for c in S.C_SCAN])
    c_or = round(c60, 3)
    ref = {
        'a0_A': a0_75,
        'gamma_soft_mode_cm-1': low(K8),
        'gamma_LO_highest_cm-1': lo_max(K8, d8),
        'eps_inf': d8['dielectric_tensor_eps_inf'][0][0],
        'born_charge_Ti': float(S.neutral_charges(d8)[1][0][0]),
        'polar_well_depth_meV_per_fu': polar_depth(8, 8),
        'R_mode_cm-1': tilt_all[2], 'tilt_angle_deg': tilt_all[0], 'tilt_well_depth_meV_per_fu': tilt_all[1],
        'strained_c_A': c75,
        'strained_inplane_polar_mode_cm-1': strained_low(c_or, 8),
    }
    ref['larger_energy_gain'] = 'polar' if ref['polar_well_depth_meV_per_fu'] > ref['tilt_well_depth_meV_per_fu'] else 'tilt'
    tol = {'absolute': {'a0_A': 0.004, 'gamma_soft_mode_cm-1': 10.0, 'gamma_LO_highest_cm-1': 4.0, 'eps_inf': 0.12,
                        'born_charge_Ti': 0.1, 'polar_well_depth_meV_per_fu': 0.5, 'R_mode_cm-1': 8.0,
                        'tilt_angle_deg': 0.3, 'tilt_well_depth_meV_per_fu': 0.5, 'strained_c_A': 0.005,
                        'strained_inplane_polar_mode_cm-1': 10.0},
           'relative': {}}
    # ---- oracle-path values (what solve.py submits from the same runs)
    oracle = {
        'a0_A': a0_60, 'gamma_soft_mode_cm-1': low(K8), 'gamma_LO_highest_cm-1': lo_max(K8, d8),
        'eps_inf': d8['dielectric_tensor_eps_inf'][0][0], 'born_charge_Ti': float(S.neutral_charges(d8)[1][0][0]),
        'polar_well_depth_meV_per_fu': polar_depth(8, 6), 'R_mode_cm-1': tilt_or[2], 'tilt_angle_deg': tilt_or[0],
        'tilt_well_depth_meV_per_fu': tilt_or[1], 'strained_c_A': c60, 'strained_inplane_polar_mode_cm-1': strained_low(c_or, 6)}
    # ---- robustness: defensible alternatives
    eos75_p = float(np.interp(0, [r['pressure_GPa'] for r in eos75][::-1], list(S.EOS_A)[::-1]))
    robust = {
        'a0_from_zero_pressure_75Ry': {'a0_A': eos75_p},
        'soft_mode_and_LO_from_6x6x6': {'gamma_soft_mode_cm-1': low(K6), 'gamma_LO_highest_cm-1': lo_max(K6, d6),
                                        'eps_inf': d6['dielectric_tensor_eps_inf'][0][0],
                                        'born_charge_Ti': float(S.neutral_charges(d6)[1][0][0])},
        'soft_mode_0.04A_displacements': {'gamma_soft_mode_cm-1': low(K6u4)},
        'LO_8x8x8_force_constants_6x6x6_tensors': {'gamma_LO_highest_cm-1': lo_max(K8, d6)},
        'polar_well_6x6x6_pattern_and_mesh': {'polar_well_depth_meV_per_fu': polar_depth(6, 6)},
        'polar_well_six_amplitudes_8x8x8': {'polar_well_depth_meV_per_fu': polar_depth(8, 8, (0.01, 0.02, 0.03, 0.04, 0.08, 0.12))},
        'polar_well_lowest_sampled_energy_8x8x8': {'polar_well_depth_meV_per_fu': -min(
            (R[f'polar_pat8_k8_t{t}']['total_energy_Ry'] - R[f'cub_a{A0}_k8']['total_energy_Ry']) * S.RY_EV * 1000
            for t in S.POLAR_TI_AMPS)},
        'tilt_angles_3_4_5': dict(zip(('tilt_angle_deg', 'tilt_well_depth_meV_per_fu', 'R_mode_cm-1'), tilt((3, 4, 5)))),
        'R_mode_from_1_2_deg': {'R_mode_cm-1': tilt((1, 2))[2]},
        'strained_c_from_stress_zero_60Ry': {'strained_c_A': stress_zero([f'strain_c{c}_k6' for c in S.C_SCAN])},
        'strained_c_from_stress_zero_75Ry': {'strained_c_A': stress_zero([f'strain_c{c}_k6_e75' for c in S.C_SCAN])},
    }
    # ---- ablations (realistic wrong paths)
    cub6 = R[f'cub_a{A0}_k6']
    offset = (R[f'tilt_a{A0}_k4_t0']['total_energy_Ry'] / 2 - cub6['total_energy_Ry']) * S.RY_EV * 1000
    Kti = S.z_block(fp(A0, 6))[1, 1]
    ti_only = float(np.sign(Kti) * math.sqrt(abs(Kti) / (S.MASS['Ti'] * S.AMU_RY)) * S.RY_CM)
    k4_depth = polar_depth(4, 4)
    abl = {
        'k4_force_constants_and_polar_well': ('dp1', {
            'gamma_soft_mode_cm-1': low(K4), 'gamma_LO_highest_cm-1': lo_max(K4, d4),
            'polar_well_depth_meV_per_fu': k4_depth,
            'larger_energy_gain': 'polar' if k4_depth > oracle['tilt_well_depth_meV_per_fu'] else 'tilt'}),
        'ti_only_slater_displacement': ('dp2', {'gamma_soft_mode_cm-1': ti_only, 'polar_well_depth_meV_per_fu': 0.0,
                                                'larger_energy_gain': 'tilt'}),
        'experimental_lattice_parameter': ('dp3', {'a0_A': P.A_EXP, 'gamma_soft_mode_cm-1': low(Kexp)}),
        'k4_dielectric_tensors': ('dp4', {'eps_inf': d4['dielectric_tensor_eps_inf'][0][0],
                                          'born_charge_Ti': float(S.neutral_charges(d4)[1][0][0]),
                                          'gamma_LO_highest_cm-1': lo_max(K8, d4)}),
        'polar_well_from_coarse_large_amplitudes': ('dp8', {'polar_well_depth_meV_per_fu': polar_depth(8, 6, (0.04, 0.08, 0.12))}),
        'analytic_frequency_as_LO': ('dp5', {'gamma_LO_highest_cm-1': float(S.modes(K8, S.Z_MASSES)[0][-1])}),
        'uncoupled_mode_by_mode_LO': ('dp5', {'gamma_LO_highest_cm-1': uncoupled_lo(K8, d8)}),
        'tilt_referenced_to_5atom_cell_other_mesh': ('dp6', {'tilt_well_depth_meV_per_fu': tilt_or[1] - offset}),
        'tilt_energies_per_10atom_cell': ('dp6', {'tilt_well_depth_meV_per_fu': 2 * tilt_or[1],
                                                  'R_mode_cm-1': math.sqrt(2) * tilt_or[2]}),
        'no_poisson_relaxation': ('dp7', {'strained_c_A': A0,
                                          'strained_inplane_polar_mode_cm-1': strained_low(A0, 6)}),
        'cubic_soft_mode_for_strained_cell': ('dp7', {'strained_c_A': A0,
                                                      'strained_inplane_polar_mode_cm-1': oracle['gamma_soft_mode_cm-1']}),
    }
    # ---- write private and solution files
    private = {'reference': ref, 'tolerances': tol,
               'provenance': {'runs': 'benchmarks/construction/fresh-srtio3-soft-modes/reference_runs.jsonl '
                                      '(task_lab.run, QE 6.7, 2 MPI ranks)',
                              'definitions': 'a0: cubic fit of 75 Ry 6x6x6 EOS; Gamma: 0.01 A frozen phonons 8x8x8 with '
                                             'ASR, dielectric 8x8x8; polar well: 8x8x8 pattern and mesh; tilt: 4x4x4 '
                                             'rotated_10, 8 angles, same-cell reference; strain: 75 Ry energy minimum '
                                             'in c, 8x8x8 in-plane frozen phonons'}}
    lab = TASK/'environment/lab/private'
    (lab/'hidden.json').write_text(json.dumps(private, indent=1) + '\n')
    hv = [v for k, v in ref.items() if k != 'larger_energy_gain']
    (lab/'hidden_values.json').write_text(json.dumps(hv, indent=1) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps(
        {'ablations': {n: {'decision_point': dp, 'changes': ch} for n, (dp, ch) in abl.items()},
         'robustness': robust, 'oracle_path': oracle}, indent=1) + '\n')
    # ---- calibration: grade each variation on top of the oracle-path payload
    sys.path.insert(0, str(TASK/'environment/lab'))
    import task_lab
    base = dict(oracle, larger_energy_gain=ref['larger_energy_gain'], evidence_ids=[1])
    trace = [{'id': 1, 'status': 'complete'}]
    ok, rep = task_lab.grade(copy.deepcopy(base), trace, [])
    print('oracle-path payload passes:', ok, {k: round(v, 4) for k, v in rep['errors'].items()})
    calib = {}
    for f in task_lab.NUMERIC:
        calib[f] = {'reference': ref[f], 'tolerance': tol['absolute'][f], 'oracle_path': oracle[f],
                    'robustness': {n: ch[f] for n, ch in robust.items() if f in ch},
                    'ablations': {n: ch[f] for n, (_, ch) in abl.items() if f in ch}}
        errs = [abs(v - ref[f]) for v in list(calib[f]['robustness'].values()) + [oracle[f]]]
        calib[f]['max_robustness_error'] = max(errs)
        ab = [abs(v - ref[f]) for v in calib[f]['ablations'].values()]
        calib[f]['nearest_ablation_error'] = min(ab) if ab else None
    calib['larger_energy_gain'] = {'reference': ref['larger_energy_gain'],
                                   'ablations': {n: ch['larger_energy_gain'] for n, (_, ch) in abl.items()
                                                 if 'larger_energy_gain' in ch}}
    for n, ch in robust.items():
        ok, rep = task_lab.grade(dict(copy.deepcopy(base), **ch), trace, [])
        print('robust', n, ok)
    for n, (dp, ch) in abl.items():
        ok, rep = task_lab.grade(dict(copy.deepcopy(base), **ch), trace, [])
        print('ablation', dp, n, ok, [k for k, v in rep['checks'].items() if not v])
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    design = json.loads((TASK/'design.json').read_text())
    design['calibration'] = {f: {k: (round(v, 5) if isinstance(v, float) else
                                     ({kk: round(vv, 5) if isinstance(vv, float) else vv for kk, vv in v.items()}
                                      if isinstance(v, dict) else v))
                                 for k, v in c.items()} for f, c in calib.items()}
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in ref.items()}, indent=1))
    print('offset 10-atom k4 vs 5-atom k6 (meV/fu):', round(offset, 3), 'c60', round(c60, 4))


if __name__ == '__main__':
    main()
