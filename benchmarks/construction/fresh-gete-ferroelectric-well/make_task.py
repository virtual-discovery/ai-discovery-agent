"""Derive the hidden references, ablation and robustness values of fresh-gete-ferroelectric-well from the reference runs.

    python3 -B make_task.py --stage-a     # after run_reference.py --stage A: writes reference_stage_a.json
    python3 -B make_task.py               # after stage B: writes private/hidden*.json, solution/ablation_values.json,
                                          # calibration.json and design.json -> calibration

The estimators are those of solution/solve.py (imported from the task package without writing bytecode).
"""
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-gete-ferroelectric-well'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
RY = S.RY_MEV
S2 = math.sqrt(2)

RUNS = {}
for line in (HERE/'reference_runs.jsonl').read_text().splitlines():
    row = json.loads(line)
    if 'result' in row:
        RUNS[row['tag']] = row['result']
    else:
        print('failed reference run:', row['tag'], row.get('error'))


def E(tag):
    return RUNS[tag]['total_energy_Ry'] * RY


def stage_a():
    out = {}
    for k in (12, 16, 20):
        eos = [RUNS[f'A_cub{a:.2f}_k{k}'] for a in (5.96, 5.98, 6.00, 6.02, 6.04, 6.06)]
        out[f'cubic_a0_k{k}_six_points'] = S.zero_pressure_a(eos)
    oracle_eos = [RUNS[f'A_cub{a:.2f}_k16'] for a in S.EOS_A]
    out['cubic_a0_A'] = round(S.zero_pressure_a(oracle_eos), 3)
    grid_tags = {k: [f'A_r3m_g{i}_k{k}' for i in range(10)] for k in (8, 12, 16)}
    for k, tags in grid_tags.items():
        pmin, H = S.model_minimum([RUNS[t] for t in tags])
        out[f'r3m_model_min_all10_k{k}'] = [float(v) for v in pmin]
    # the oracle's seven k12 points (grid indices of solve.R3M_POINTS)
    index = {(0, 0, 0): 0, (1, 0, 0): 1, (-1, 0, 0): 2, (0, 1, 0): 3, (0, -1, 0): 4, (0, 0, 1): 5, (0, 0, -1): 6,
             (1, 1, 0): 7, (1, 0, 1): 8, (0, 1, 1): 9}
    oracle = [RUNS[f'A_r3m_g{index[p]}_k12'] for p in S.R3M_POINTS]
    for res, req in zip(oracle, S.r3m_grid_requests()):
        assert all(res['settings'][k] == req[k] for k in req), 'oracle grid request differs from the reference grid'
    pmin, _ = S.model_minimum(oracle)
    out['r3m_min'] = list(S.rounded_minimum(pmin))
    out['r3m_min_volume_A3'] = S.volume(out['r3m_min'][0], out['r3m_min'][1])
    (HERE/'reference_stage_a.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


FIELDS = ['cubic_a0_A', 'r3m_a_rh_A', 'r3m_alpha_deg', 'r3m_x_Te', 'well_depth_meV_per_fu',
          'clamped_well_depth_meV_per_fu', 'soft_mode_critical_pressure_GPa', 'soc_well_depth_change_meV_per_fu']
TOL_ABS = {'cubic_a0_A': 0.008, 'r3m_a_rh_A': 0.005, 'r3m_alpha_deg': 0.1, 'r3m_x_Te': 0.0008,
           'soft_mode_critical_pressure_GPa': 0.6, 'soc_well_depth_change_meV_per_fu': 0.6}
TOL_REL = {'well_depth_meV_per_fu': 0.04, 'clamped_well_depth_meV_per_fu': 0.05}


def grid(k):
    return [RUNS[f'A_r3m_g{i}_k{k}'] for i in range(10)]


def soft(prefix_suffix, alats=(5.80, 5.75, 5.70, 5.65)):
    return S.soft_mode_pressure([RUNS[f'B_soft{a:.2f}_{prefix_suffix}'] for a in alats])[0]


def clamped(k, xs=(0.475, 0.480, 0.485, 0.490), cubic=None):
    return S.clamped_depth(RUNS[cubic or f'B_cub_a0_k{k}'], [RUNS[f'B_clamp_x{x:.3f}_k{k}'] for x in xs])[0]


def stage_b():
    A = json.loads((HERE/'reference_stage_a.json').read_text())
    a0o = A['cubic_a0_A']
    # ---- reference (densest data)
    a0_ref = A['cubic_a0_k20_six_points']
    _, H16 = S.model_minimum(grid(16))
    p24, e24 = S.newton_refine(RUNS['B_r3m_min_k24'], H16)
    depth_ref = E('B_cub_a0_k24') - e24
    clamped_ref = clamped(20)
    pc_ref = soft('dx0.005_k24_smearing0.002')
    soc = lambda k: ((E(f'B_pslibrary_fullrel_cub_a0_k{k}') - E(f'B_pslibrary_fullrel_r3m_min_k{k}'))      # noqa: E731
                     - (E(f'B_pslibrary_cub_a0_k{k}') - E(f'B_pslibrary_r3m_min_k{k}')))
    ref = {'cubic_a0_A': a0_ref, 'r3m_a_rh_A': p24[0], 'r3m_alpha_deg': p24[1], 'r3m_x_Te': p24[2],
           'well_depth_meV_per_fu': depth_ref, 'clamped_well_depth_meV_per_fu': clamped_ref,
           'soft_mode_critical_pressure_GPa': pc_ref, 'soc_well_depth_change_meV_per_fu': soc(16)}
    ref = {k: float(v) for k, v in ref.items()}
    # ---- robustness (defensible alternatives)
    _, H12o = S.model_minimum([RUNS[f'A_r3m_g{i}_k12'] for i in range(7)])
    p20o, e20o = S.newton_refine(RUNS['B_r3m_min_k20'], H12o)
    p16g, _ = S.model_minimum(grid(16))
    p20, e20 = S.newton_refine(RUNS['B_r3m_min_k20'], H16)
    rob = {
        'oracle_route': {'cubic_a0_A': a0o, 'r3m_a_rh_A': p20o[0], 'r3m_alpha_deg': p20o[1], 'r3m_x_Te': p20o[2],
                         'well_depth_meV_per_fu': E('B_cub_a0_k24') - e20o,
                         'clamped_well_depth_meV_per_fu': clamped(20, (0.480, 0.485)),
                         'soft_mode_critical_pressure_GPa': soft('dx0.005_k20_smearing0.005', (5.80, 5.75, 5.70)),
                         'soc_well_depth_change_meV_per_fu': soc(12)},
        'k16_grid_model_minimum_structure': {'r3m_a_rh_A': p16g[0], 'r3m_alpha_deg': p16g[1], 'r3m_x_Te': p16g[2]},
        'k20_both_phases': {'well_depth_meV_per_fu': E('B_cub_a0_k20') - e20},
        'eos_k16_six_points': {'cubic_a0_A': A['cubic_a0_k16_six_points']},
        'clamped_k16_four_points': {'clamped_well_depth_meV_per_fu': clamped(16)},
        'soft_mode_k20_mv0.002': {'soft_mode_critical_pressure_GPa': soft('dx0.005_k20_smearing0.002')},
        'soft_mode_k20_mv0.005_all_points': {'soft_mode_critical_pressure_GPa': soft('dx0.005_k20_smearing0.005')},
    }
    # ---- ablations (realistic wrong paths, from the same reference calculations)
    p12, e12 = S.newton_refine(RUNS['B_r3m_min_k12'], H12o)
    # quadratic model through one-sided steps from the starting guess, minimum taken without a further SCF there
    one_sided = [RUNS[f'A_r3m_g{i}_k12'] for i in (0, 1, 3, 5, 7, 8, 9)]
    p1s, _ = S.model_minimum(one_sided)
    centre, e0, b, Hs = S.quadratic_model(one_sided)
    u = (p1s - centre) / S.SCALE
    e1s = e0 + b @ u + 0.5 * u @ Hs @ u + (E('A_r3m_g0_k20') - E('A_r3m_g0_k12'))
    exp_x = [RUNS[f'A_exp_cell_x{x:.3f}_k16'] for x in (0.470, 0.475, 0.480)]
    xs = np.array([r['settings']['x_Te'] for r in exp_x]); es = np.array([r['total_energy_Ry'] * RY for r in exp_x])
    c = np.polyfit(xs, es, 2); x_exp = -c[1] / (2 * c[0]); e_exp = np.polyval(c, x_exp)
    abl = {
        'clamped_reported_as_well_depth': {'well_depth_meV_per_fu': clamped(20)},
        'cubic_at_r3m_volume': {'well_depth_meV_per_fu': E('B_cub_at_r3m_volume_k20') - e20},
        'x_relaxed_in_experimental_cell': {'well_depth_meV_per_fu': E('B_cub_a0_k16') - e_exp, 'r3m_a_rh_A': 4.28,
                                            'r3m_alpha_deg': 58.36, 'r3m_x_Te': float(x_exp)},
        'k12_mesh_both_phases': {'well_depth_meV_per_fu': E('B_cub_a0_k12') - e12,
                                 'clamped_well_depth_meV_per_fu': clamped(12)},
        'one_sided_quadratic_model': {'r3m_a_rh_A': p1s[0], 'r3m_alpha_deg': p1s[1], 'r3m_x_Te': p1s[2],
                                      'well_depth_meV_per_fu': E('B_cub_a0_k24') - e1s},
        'per_atom_energies': {'well_depth_meV_per_fu': depth_ref / 2, 'clamped_well_depth_meV_per_fu': clamped_ref / 2},
        'fixed_occupations_soft_mode_k16': {'soft_mode_critical_pressure_GPa': soft('dx0.005_k16_fixed')},
        'soft_mode_k12_mv0.01': {'soft_mode_critical_pressure_GPa': soft('dx0.005_k12_smearing0.01')},
        'soc_compared_with_sssp': {'soc_well_depth_change_meV_per_fu':
                                   (E('B_pslibrary_fullrel_cub_a0_k12') - E('B_pslibrary_fullrel_r3m_min_k12'))
                                   - (E('B_cub_a0_k12') - E('B_r3m_min_k12'))},
    }
    # large displacement (x = 0.48): force-based "force constant" includes the quartic term
    big = [RUNS[f'B_soft{a:.2f}_dx0.02_k20_smearing0.005'] for a in (5.80, 5.75, 5.70)]
    # the force at x = 0.48 includes the quartic term: the apparent force constant is positive at every pressure
    # bought, and its linear extrapolation to zero from the two lowest pressures gives the wrong path's estimate
    pts = sorted((r['pressure_GPa'], S.gradient(r)[2] / (r['settings']['x_Te'] - 0.5)) for r in big)
    (p1, k1), (p2, k2) = pts[:2]
    abl['soft_mode_from_dx0.02'] = {'soft_mode_critical_pressure_GPa': p1 - k1 * (p2 - p1) / (k2 - k1)}
    to_float = lambda d: {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in d.items()}  # noqa: E731
    rob = {k: to_float(v) for k, v in rob.items()}
    abl = {k: to_float(v) for k, v in abl.items()}
    report = {'reference': ref, 'robustness': rob, 'ablations': abl}
    print(json.dumps(report, indent=1))

    def err(field, value):
        return value / ref[field] - 1 if field in TOL_REL else value - ref[field]

    def passes(field, value):
        return abs(err(field, value)) <= (TOL_REL.get(field) or TOL_ABS[field])
    calib = {}
    for f in FIELDS:
        rv = {n: d[f] for n, d in rob.items() if f in d}
        av = {n: d[f] for n, d in abl.items() if f in d and d[f] is not None}
        calib[f] = {'reference': ref[f], 'tolerance': f'+-{TOL_ABS[f]}' if f in TOL_ABS else f'within {100 * TOL_REL[f]:.1f}%',
                    'robustness': rv, 'ablations': av,
                    'max_robust_error': max([abs(err(f, v)) for v in rv.values()], default=None),
                    'robust_all_pass': all(passes(f, v) for v in rv.values()),
                    'ablation_errors': {n: err(f, v) for n, v in av.items()},
                    'failing_ablations': sorted(n for n, v in av.items() if not passes(f, v))}
    for n, d in abl.items():
        fails = [f for f in FIELDS if f in d and d[f] is not None and not passes(f, d[f])]
        print('ablation', n, 'fails on', fails)
    for f, c in calib.items():
        print(f, 'ref', round(c['reference'], 5), 'max robust err', c['max_robust_error'], 'all robust pass', c['robust_all_pass'])
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    hidden = {'reference': ref, 'tolerances': {'absolute': TOL_ABS, 'relative': TOL_REL},
              'provenance': {'runs': 'benchmarks/construction/fresh-gete-ferroelectric-well/reference_runs.jsonl '
                                     '(task_lab.run, QE 6.7, 2 MPI ranks)'}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps([ref[f] for f in FIELDS]) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    design = json.loads((TASK/'design.json').read_text())
    design['calibration'] = calib
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')


if __name__ == '__main__':
    if '--stage-a' in sys.argv:
        stage_a()
    else:
        stage_b()
