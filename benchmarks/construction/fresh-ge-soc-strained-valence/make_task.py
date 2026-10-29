"""Derive the hidden targets, ablation and robustness values and the calibration table for fresh-ge-soc-strained-valence.

    .venv/bin/python -B benchmarks/construction/fresh-ge-soc-strained-valence/make_task.py

Reads reference_runs.jsonl (task_lab.run in the runtime image; stages A/B from reference_plan.py, stage C = the
oracle's own chain from oracle_chain.py), applies the estimators of solution/solve.py, and writes
environment/lab/private/hidden.json, private/hidden_values.json, solution/ablation_values.json and design.json ->
calibration (plus calibration.json here).
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-ge-soc-strained-valence'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

ROWS = {}
for line in (HERE/'reference_runs.jsonl').read_text().splitlines():
    if line.strip():
        r = json.loads(line)
        if 'result' in r:
            ROWS[r['tag']] = r['result']
LAT = json.loads((HERE/'reference_lattice.json').read_text())
A0 = LAT['A0']
A_EXP = 5.658
TOL = {'a0_A': 0.005, 'delta0_at_a0_eV': 0.005, 's_level_at_a0_eV': 0.035, 'c_over_a_at_minus1pct': 0.0008,
       'hh_minus_lh_at_minus1pct_eV': 0.005, 'so_minus_hh_at_minus1pct_eV': 0.008,
       'b_shear_deformation_potential_eV': 0.12, 's_so_crossing_strain_percent': 0.1}


def pairs(e):
    return [(e[i] + e[i + 1]) / 2 for i in range(0, len(e), 2)]


def p_root(tags, deg=None):
    a = np.array([ROWS[t]['settings']['a_par_A'] for t in tags])
    p = np.array([ROWS[t]['pressure_GPa'] for t in tags])
    c = np.polyfit(a, p, deg if deg is not None else len(a) - 1)
    r = np.roots(c).real
    return float(r[np.argmin(abs(r - 5.77))])


def relaxed(e, mode='soc'):
    """Out-of-plane strain with sigma_zz = 0 at in-plane strain e (linear in the two scan points) and, for
    spin-orbit scans, the labelled Gamma levels interpolated to it."""
    cz = LAT['scan_c'][str(e)]
    tags = [f'B_scan_{mode}_e{e}_c{c}' for c in cz]
    ez = S.interp_root(cz, [ROWS[t]['stress_GPa'][2][2] for t in tags])
    if mode != 'soc':
        return ez, None, tags
    Ls = [S.gamma_levels(ROWS[t]) for t in tags]
    f = (ez - cz[0]) / (cz[1] - cz[0])
    return ez, {k: Ls[0][k] + f * (Ls[1][k] - Ls[0][k]) for k in Ls[0]}, tags


def scalar_levels(result):
    """Scalar-relativistic Gamma levels (eV): s-like singlet (odd parity) and the p-like states (even)."""
    levs = result['band_kpoints'][0]['symmetry']['levels']
    out = {'s': None, 'p': []}
    for lev in levs:
        i, j = lev['bands']
        if 12 <= i <= 15:
            name = lev['label']
            if 'A_2u' in name or 'B_1u' in name:
                out['s'] = lev['energy_eV']
            else:
                out['p'] += [lev['energy_eV']] * (j - i + 1)
    return out


def bp6(D0):
    def split(Q):
        r = np.sqrt(D0 * D0 + 2 * D0 * Q + 9 * Q * Q)
        return -Q, 0.5 * (Q - D0 + r), 0.5 * (Q - D0 - r)
    return split


def b_from_6x6(dE, x, D0):
    split = bp6(D0)
    Q = brentq(lambda q: split(q)[0] - split(q)[1] - dE, -0.3, 0.3)
    return -2 * Q / x


def main():
    ref, abl, rob = {}, {}, {}
    # ---------------- a0
    ref['a0_A'] = p_root([f'A_soc_k10_ec80_{a}' for a in (5.75, 5.77, 5.79)])
    rob['a0_spin_orbit_60Ry_k8'] = {'a0_A': p_root(['A_soc_k8_ec60_5.75', 'A_soc_k8_ec60_5.79'])}
    rob['a0_scalar_relativistic_80Ry'] = {'a0_A': p_root([f'A_scalar_k10_ec80_{a}' for a in (5.75, 5.77, 5.79)])}
    abl['pressures_at_40_Ry'] = {'a0_A': p_root(['A_soc_k8_ec40_5.75', 'A_soc_k8_ec40_5.79'])}
    abl['experimental_lattice_parameter'] = {'a0_A': A_EXP}
    # ---------------- levels at a0
    L0 = S.gamma_levels(ROWS['B_a0_bands_k12'])
    ref['delta0_at_a0_eV'] = L0['hh'] - L0['so']
    ref['s_level_at_a0_eV'] = L0['s'] - L0['hh']
    for tag, name in (('B_a0_bands_k8', 'levels_k8_60Ry'), ('B_a0_k6_ec50', 'levels_k6_50Ry_mv001')):
        L = S.gamma_levels(ROWS[tag])
        rob[name] = {'delta0_at_a0_eV': L['hh'] - L['so'], 's_level_at_a0_eV': L['s'] - L['hh']}
    pf = pairs(ROWS['B_a0_gamma_fixed']['mesh_point_nearest_gamma']['band_energies_eV'])  # Gamma, fixed occupations
    rob['levels_fixed_occupations_by_degeneracy'] = {'delta0_at_a0_eV': pf[13] - pf[12], 's_level_at_a0_eV': pf[11] - pf[13]}
    e0 = ROWS['B_a0_bands_k12']['band_kpoints'][0]['band_energies_eV']
    abl['first_level_above_gamma8_as_s'] = {'s_level_at_a0_eV': e0[30] - L0['hh']}
    abl['band_count_identification'] = {'delta0_at_a0_eV': e0[27] - e0[23]}   # top occupied (28) minus band 24
    ps = pairs(ROWS['B_a0_shift1']['mesh_point_nearest_gamma']['band_energies_eV'])
    abl['shifted_mesh_nearest_point'] = {'delta0_at_a0_eV': ps[13] - ps[12], 's_level_at_a0_eV': ps[11] - ps[14]}
    sc = scalar_levels(ROWS['B_a0_scalar_bands'])
    abl['spin_orbit_off_s_level'] = {'s_level_at_a0_eV': sc['s'] - max(sc['p'])}
    # ---------------- -1 % relaxed
    e1 = -0.01
    cz = LAT['scan_c'][str(e1)]
    La, Lb = [S.gamma_levels(ROWS[f'B_scan_soc_e{e1}_c{c}']) for c in cz]
    slope = {k: (Lb[k] - La[k]) / (cz[1] - cz[0]) for k in La}          # level shifts per unit eps_zz (60 Ry k8)
    s60 = [ROWS[f'B_scan_soc_e{e1}_c{c}']['stress_GPa'][2][2] for c in cz]
    dsdz = (s60[1] - s60[0]) / (cz[1] - cz[0])                          # d sigma_zz / d eps_zz (GPa)

    def at(tag, ez_target=None):
        """sigma_zz = 0 strain from one run plus the scan slope, and that run's levels moved to it."""
        r = ROWS[tag]
        cz0 = r['settings']['c_A'] / A0 - 1
        ez = cz0 - r['stress_GPa'][2][2] / dsdz if ez_target is None else ez_target
        L = S.gamma_levels(r)
        return ez, {k: L[k] + slope[k] * (ez - cz0) for k in L}
    # reference: the most converged run (Gamma-k10, 80 Ry) for both the stress zero and the levels
    ez1, L1 = at('B_m1_k10_ec80')
    ref['c_over_a_at_minus1pct'] = (1 + ez1) / (1 + e1)
    ref['hh_minus_lh_at_minus1pct_eV'] = L1['hh'] - L1['lh']
    ref['so_minus_hh_at_minus1pct_eV'] = L1['so'] - L1['hh']
    put = lambda ez, L: {'c_over_a_at_minus1pct': (1 + ez) / (1 + e1), 'hh_minus_lh_at_minus1pct_eV': L['hh'] - L['lh'],   # noqa: E731
                         'so_minus_hh_at_minus1pct_eV': L['so'] - L['hh']}
    ez6, L6 = relaxed(e1)[:2]
    rob['relaxation_and_levels_60Ry_k8'] = put(ez6, L6)
    rob['relaxation_and_levels_80Ry_k8'] = put(*at('B_m1_k8_ec80'))
    ezs, _, stags = relaxed(e1, 'scalar')
    rob['scalar_relativistic_relaxation_60Ry_k8'] = put(ezs, {k: La[k] + slope[k] * (ezs - cz[0]) for k in La})
    Lc0, Lcub = S.gamma_levels(ROWS['B_m1_c_a0']), S.gamma_levels(ROWS['B_m1_cubic'])
    abl['c_fixed_at_a0'] = {'c_over_a_at_minus1pct': 1 / (1 + e1), 'hh_minus_lh_at_minus1pct_eV': Lc0['hh'] - Lc0['lh'],
                            'so_minus_hh_at_minus1pct_eV': Lc0['so'] - Lc0['hh']}
    abl['cubic_compression'] = {'c_over_a_at_minus1pct': 1.0, 'hh_minus_lh_at_minus1pct_eV': Lcub['hh'] - Lcub['lh'],
                                'so_minus_hh_at_minus1pct_eV': Lcub['so'] - Lcub['hh']}
    sa, sb = [scalar_levels(ROWS[t]) for t in stags]
    f = (ezs - LAT['scan_c'][str(e1)][0]) / (LAT['scan_c'][str(e1)][1] - LAT['scan_c'][str(e1)][0])
    trip = [x + f * (y - x) for x, y in zip(sorted(sa['p'], reverse=True), sorted(sb['p'], reverse=True))]
    abl['spin_orbit_off_splitting'] = {'hh_minus_lh_at_minus1pct_eV': trip[0] - trip[2]}
    abl['koster_label_for_heavy_hole'] = {'hh_minus_lh_at_minus1pct_eV': L1['lh'] - L1['hh'],
                                          'so_minus_hh_at_minus1pct_eV': L1['so'] - L1['lh']}
    # ---------------- b
    split = {}
    for e in (0.005, -0.005):
        ez, L, _ = relaxed(e)
        split[e] = (L['hh'] - L['lh'], 2 * e - 2 * ez)
    ref['b_shear_deformation_potential_eV'] = S.b_symmetric(*split[0.005], *split[-0.005])
    x1 = 2 * e1 - 2 * ez1
    D0 = ref['delta0_at_a0_eV']
    rob['b_6x6_bir_pikus_from_minus1pct'] = {'b_shear_deformation_potential_eV': b_from_6x6(ref['hh_minus_lh_at_minus1pct_eV'], x1, D0)}
    rob['b_6x6_bir_pikus_from_pm0.5pct'] = {'b_shear_deformation_potential_eV': float(np.mean([b_from_6x6(d, x, D0) for d, x in split.values()]))}
    X = np.array([split[0.005][1], split[-0.005][1], x1])
    Y = np.array([split[0.005][0], split[-0.005][0], ref['hh_minus_lh_at_minus1pct_eV']])
    rob['b_quadratic_fit_three_strains'] = {'b_shear_deformation_potential_eV': float(np.polyfit(X, Y, 2)[1])}
    abl['linear_bir_pikus_at_minus1pct'] = {'b_shear_deformation_potential_eV': ref['hh_minus_lh_at_minus1pct_eV'] / x1}
    abl['linear_formula_at_minus0.5pct'] = {'b_shear_deformation_potential_eV': split[-0.005][0] / split[-0.005][1]}
    abl['koster_label_for_heavy_hole']['b_shear_deformation_potential_eV'] = -ref['b_shear_deformation_potential_eV']
    abl['c_fixed_at_a0']['b_shear_deformation_potential_eV'] = (Lc0['hh'] - Lc0['lh']) / (2 * e1)   # eps_zz taken as 0
    # ---------------- crossing
    track = {0.0: L0['s'] - L0['so']}
    track[0.0] = S.gamma_levels(ROWS['B_a0_bands_k8'])['s'] - S.gamma_levels(ROWS['B_a0_bands_k8'])['so']
    for e in (-0.01, -0.02, -0.025, -0.03):
        track[e] = relaxed(e)[1]['s'] - relaxed(e)[1]['so']
    xs = np.array(sorted(track))
    ys = np.array([track[x] for x in xs])
    near = sorted(track, key=lambda x: abs(track[x]))[:3]
    c2 = np.polyfit(near, [track[x] for x in near], 2)
    r = np.roots(c2).real
    ref['s_so_crossing_strain_percent'] = 100 * float(r[np.argmin(abs(r + 0.02))])
    rob['crossing_linear_between_bracketing_strains'] = {'s_so_crossing_strain_percent': 100 * S.interp_root([-0.02, -0.01], [track[-0.02], track[-0.01]])}
    rob['crossing_quadratic_0_to_-2pct'] = {'s_so_crossing_strain_percent': 100 * float(
        [x for x in np.roots(np.polyfit([0, -0.01, -0.02], [track[0], track[-0.01], track[-0.02]], 2)).real if -0.03 < x < -0.01][0])}
    abl['crossing_linear_extrapolation_from_0_and_-1pct'] = {'s_so_crossing_strain_percent': 100 * S.interp_root([0.0, -0.01], [track[0.0], track[-0.01]])}
    ec = ref['s_so_crossing_strain_percent'] / 100
    # strain reference = experimental lattice parameter: the -1 % state is then eps_par = 0.99 A_EXP/A0 - 1 of this protocol
    ex = 0.99 * A_EXP / A0 - 1
    e_lo, e_hi = -0.03, -0.025
    (za, La_), (zb, Lb_) = [relaxed(e)[:2] for e in (e_lo, e_hi)]
    g = (ex - e_lo) / (e_hi - e_lo)
    Lx = {k: La_[k] + g * (Lb_[k] - La_[k]) for k in La_}
    zx = za + g * (zb - za)
    abl['strain_from_experimental_lattice_parameter'] = {
        'c_over_a_at_minus1pct': (1 + zx) / (1 + ex), 'hh_minus_lh_at_minus1pct_eV': Lx['hh'] - Lx['lh'],
        'so_minus_hh_at_minus1pct_eV': Lx['so'] - Lx['hh'], 's_so_crossing_strain_percent': 100 * ((1 + ec) * A0 / A_EXP - 1)}
    # ---------------- outputs
    ref = {k: float(v) for k, v in ref.items()}
    order = [('dp1', 'first_level_above_gamma8_as_s'), ('dp1', 'band_count_identification'),
             ('dp2', 'spin_orbit_off_s_level'), ('dp2', 'spin_orbit_off_splitting'),
             ('dp3', 'strain_from_experimental_lattice_parameter'), ('dp3', 'experimental_lattice_parameter'),
             ('dp4', 'c_fixed_at_a0'), ('dp4', 'cubic_compression'),
             ('dp5', 'linear_bir_pikus_at_minus1pct'), ('dp5', 'linear_formula_at_minus0.5pct'),
             ('dp6', 'shifted_mesh_nearest_point'), ('dp7', 'koster_label_for_heavy_hole'),
             ('dp8', 'crossing_linear_extrapolation_from_0_and_-1pct'),
             ('dp9', 'pressures_at_40_Ry')]
    to_f = lambda d: {k: float(v) for k, v in d.items()}   # noqa: E731
    abl = {k: to_f(v) for k, v in abl.items()}
    rob = {k: to_f(v) for k, v in rob.items()}
    (TASK/'solution/ablation_values.json').write_text(json.dumps(
        {'ablation_order': order, 'ablations': abl, 'robustness': rob,
         'source': 'benchmarks/construction/fresh-ge-soc-strained-valence/make_task.py from reference_runs.jsonl'}, indent=1) + '\n')
    private = {'reference': ref, 'tolerances': TOL,
               'provenance': {'runs': 'benchmarks/construction/fresh-ge-soc-strained-valence/reference_runs.jsonl '
                                      '(task_lab.run, QE 6.7, 2 MPI ranks)', 'A0_of_stage_B': A0}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(private, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(list(ref.values())) + '\n')
    cal = {}
    for f_ in ref:
        rr = {n: v[f_] for n, v in rob.items() if f_ in v}
        aa = {n: v[f_] for n, v in abl.items() if f_ in v}
        cal[f_] = {'reference': ref[f_], 'tolerance': f'+-{TOL[f_]}', 'robustness': rr, 'ablations': aa,
                   'max_robust_error': max([abs(v - ref[f_]) for v in rr.values()], default=None),
                   'nearest_ablation_error': min([abs(v - ref[f_]) for v in aa.values()], default=None)}
    cal['informational'] = {'D_minus_eperp_over_epar': {str(e): -relaxed(e)[0] / e for e in LAT['scan_strains']},
                            's_minus_so_track_eV': {str(k): float(v) for k, v in track.items()},
                            'b_one_sided_split': {str(k): [float(v[0]), float(v[1])] for k, v in split.items()}}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1) + '\n')
    design = json.loads((TASK/'design.json').read_text())
    design['calibration'] = cal
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    for f_, c in cal.items():
        if f_ != 'informational':
            print(f"{f_:36s} ref {c['reference']:.5f} tol {c['tolerance']:8s} robust {c['max_robust_error']} nearest abl {c['nearest_ablation_error']}")
    print(json.dumps(cal['informational'], indent=1))
    bad = [(n, f_, v) for n, d in abl.items() for f_, v in d.items() if abs(v - ref[f_]) <= TOL[f_]]
    print('ablation fields inside tolerance (fine only if another field of the same ablation fails):', bad)
    badr = [(n, f_, v) for n, d in rob.items() for f_, v in d.items() if abs(v - ref[f_]) > TOL[f_]]
    print('robustness fields outside tolerance:', badr)


if __name__ == '__main__':
    main()
