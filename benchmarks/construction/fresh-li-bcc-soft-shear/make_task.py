"""Build the private reference, ablation and robustness values for fresh-li-bcc-soft-shear.

    python3 -B benchmarks/construction/fresh-li-bcc-soft-shear/make_task.py

Reads the reference calculations made by run_batch.sh (each one a call of the task's own task_lab.run, stored under
~/.cache/ai-theorist/fresh-li-bcc-soft-shear/runs/<name>/result.json), copies their results into
reference_results.json (kept in this folder), and writes
  environment/lab/private/hidden.json, environment/lab/private/hidden_values.json,
  solution/ablation_values.json, calibration.json
and prints the calibration table. The analysis functions are the oracle's own (solution/solve.py).
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-li-bcc-soft-shear'
RUNS = Path.home()/'.cache/ai-theorist/fresh-li-bcc-soft-shear/runs'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
spec = importlib.util.spec_from_file_location('refjobs', HERE/'reference_jobs.py')
J = importlib.util.module_from_spec(spec)
spec.loader.exec_module(J)

A0, AC, AX = J.A0, J.AC, J.AX
KEEP = ('lattice_A', 'strain_voigt', 'cell_volume_A3', 'kmesh', 'smearing', 'ecutwfc_Ry', 'free_energy_Ry',
        'internal_energy_Ry', 'stress_GPa', 'pressure_GPa', 'fermi_energy_eV', 'irreducible_k_points',
        'symmetry_operations', 'scf_iterations', 'wall_seconds')
TOL = {'absolute': {'a0_A': 0.005, 'Cprime_GPa': 0.065, 'C44_GPa': 0.42, 'compressed_pressure_GPa': 0.05,
                    'compressed_Bprime_GPa': 0.09, 'compressed_B44_GPa': 0.56},
       'relative': {'bulk_modulus_GPa': 0.03, 'zener_anisotropy': 0.06}}
BOOLEANS = ['mechanically_stable_at_a0']
FIELDS = ['a0_A', 'bulk_modulus_GPa', 'Cprime_GPa', 'C44_GPa', 'zener_anisotropy', 'mechanically_stable_at_a0',
          'compressed_pressure_GPa', 'compressed_Bprime_GPa', 'compressed_B44_GPa']


def bulk_modulus_at(lattice, energy_ry, a):
    """B = V d2F/dV2 at lattice parameter a from the third-order Birch-Murnaghan polynomial in V^(-2/3)."""
    v = np.array(lattice, float) ** 3 / 2.0
    c = np.polyfit(v ** (-2.0 / 3.0), energy_ry, 3)
    vv = a ** 3 / 2.0
    x = vv ** (-2.0 / 3.0)
    d2 = np.polyval(np.polyder(c, 2), x) * (2.0 / 3.0 * vv ** (-5.0 / 3.0)) ** 2 \
        + np.polyval(np.polyder(c, 1), x) * (10.0 / 9.0 * vv ** (-8.0 / 3.0))
    return float(vv * d2 * S.RY_EV * S.EVA3_GPA)


def load_results():
    names = [n for jobs in J.batches().values() for n, _ in jobs]
    store = HERE/'reference_results.json'
    have = json.loads(store.read_text()) if store.exists() else {}
    for name in names:
        path = RUNS/name/'result.json'
        if path.exists():
            r = json.loads(path.read_text())
            have[name] = {k: r[k] for k in KEEP}
    missing = [n for n in names if n not in have]
    if missing:
        sys.exit('missing reference runs: %s' % missing)
    store.write_text(json.dumps({n: have[n] for n in names}, indent=0) + '\n')
    return have


def main():
    R = load_results()
    c = lambda a, k, s: R[f'c_a{a}_k{k}_s{s}']
    tet = lambda a, k, s, d=0.01: (R[f't_a{a}_k{k}_s{s}_d+{d:g}'], R[f't_a{a}_k{k}_s{s}_d-{d:g}'])
    tri = lambda a, k, s, e=0.01: (R[f'r_a{a}_k{k}_s{s}_e+{e:g}'], R[f'r_a{a}_k{k}_s{s}_e-{e:g}'])
    cp = lambda a, k, s, d=0.01: S.cprime_stress(*tet(a, k, s, d))
    c44 = lambda a, k, s, e=0.01: S.c44_stress(*tri(a, k, s, e))

    # ---------------- reference
    eos_a = [3.39, 3.41, 3.43, 3.45, 3.47, 3.49]
    a0_fit, b0 = S.birch_murnaghan(eos_a, [c(a, 48, 0.01)['free_energy_Ry'] for a in eos_a])
    assert abs(a0_fit - A0) < 1e-4, (a0_fit, A0)
    p = [c(a, 48, 0.01)['pressure_GPa'] for a in eos_a]
    a0_p0 = [float(r.real) for r in np.roots(np.polyfit(eos_a, p, 3)) if abs(r.imag) < 1e-9 and 3.3 < r.real < 3.6][0]
    cp_set = {f'k{k}/{s}': cp(A0, k, s) for k, s in [(48, 0.01), (56, 0.01), (48, 0.005), (56, 0.005)]}
    c44_set = {f'k{k}/{s}': c44(A0, k, s) for k in (48, 56, 64) for s in (0.01, 0.005)}
    pc_set = {f'k{k}/{s}': c(AC, k, s)['pressure_GPa'] for k, s in [(48, 0.01), (56, 0.01), (64, 0.01), (48, 0.005), (56, 0.005)]}
    bpc_set = {f'k{k}/{s}': cp(AC, k, s) for k, s in [(48, 0.01), (56, 0.01), (64, 0.01), (48, 0.005), (56, 0.005)]}
    b44c_set = {f'k{k}/{s}': c44(AC, k, s) for k in (48, 56, 64) for s in (0.01, 0.005)}
    mean = lambda d: float(np.mean(list(d.values())))
    ref = {'a0_A': round(A0, 4), 'bulk_modulus_GPa': b0, 'Cprime_GPa': mean(cp_set), 'C44_GPa': mean(c44_set),
           'mechanically_stable_at_a0': True, 'compressed_pressure_GPa': mean(pc_set),
           'compressed_Bprime_GPa': mean(bpc_set), 'compressed_B44_GPa': mean(b44c_set)}
    ref['zener_anisotropy'] = ref['C44_GPa'] / ref['Cprime_GPa']
    ref = {f: ref[f] for f in FIELDS}

    def with_(**changes):
        out = dict(ref)
        out.update(changes)
        if 'Cprime_GPa' in changes or 'C44_GPa' in changes:
            out['zener_anisotropy'] = out['C44_GPa'] / out['Cprime_GPa']
            out['mechanically_stable_at_a0'] = bool(out['Cprime_GPa'] > 0 and out['C44_GPa'] > 0)
        return out

    # ---------------- robustness (defensible alternatives; must pass)
    rob = {}
    a0_32, b0_32 = S.birch_murnaghan([3.39, 3.41, 3.43, 3.45, 3.47],
                                     [c(a, 32, 0.01)['free_energy_Ry'] for a in [3.39, 3.41, 3.43, 3.45, 3.47]])
    rob['oracle_like_k48_both_smearings'] = with_(a0_A=round(a0_32, 4), bulk_modulus_GPa=b0_32,
                                                 Cprime_GPa=float(np.mean([cp(A0, 48, 0.01), cp(A0, 48, 0.005)])),
                                                 C44_GPa=c44(A0, 48, 0.01), compressed_pressure_GPa=c(AC, 48, 0.01)['pressure_GPa'],
                                                 compressed_Bprime_GPa=cp(AC, 48, 0.01), compressed_B44_GPa=c44(AC, 48, 0.01))
    rob['k56_smearing_0.01'] = with_(Cprime_GPa=cp(A0, 56, 0.01), C44_GPa=c44(A0, 56, 0.01),
                                     compressed_pressure_GPa=c(AC, 56, 0.01)['pressure_GPa'],
                                     compressed_Bprime_GPa=cp(AC, 56, 0.01), compressed_B44_GPa=c44(AC, 56, 0.01))
    rob['k56_smearing_0.005'] = with_(Cprime_GPa=cp(A0, 56, 0.005), C44_GPa=c44(A0, 56, 0.005),
                                      compressed_pressure_GPa=c(AC, 56, 0.005)['pressure_GPa'],
                                      compressed_Bprime_GPa=cp(AC, 56, 0.005), compressed_B44_GPa=c44(AC, 56, 0.005))
    rob['k64_smearing_0.01'] = with_(C44_GPa=c44(A0, 64, 0.01), compressed_Bprime_GPa=cp(AC, 64, 0.01),
                                     compressed_B44_GPa=c44(AC, 64, 0.01))
    rob['smearing_0.005_k48'] = with_(Cprime_GPa=cp(A0, 48, 0.005), C44_GPa=c44(A0, 48, 0.005),
                                      compressed_pressure_GPa=c(AC, 48, 0.005)['pressure_GPa'],
                                      compressed_Bprime_GPa=cp(AC, 48, 0.005), compressed_B44_GPa=c44(AC, 48, 0.005))
    rob['energy_route_volume_conserving'] = with_(
        Cprime_GPa=S.cprime_energy(c(A0, 48, 0.01), *tet(A0, 48, 0.01)),
        compressed_Bprime_GPa=S.cprime_energy(c(AC, 48, 0.01), *tet(AC, 48, 0.01)),
        compressed_B44_GPa=S.c44_energy_linear(c(AC, 48, 0.01), *tri(AC, 48, 0.01)) - c(AC, 48, 0.01)['pressure_GPa'] / 2)
    rob['half_percent_strains'] = with_(Cprime_GPa=cp(A0, 48, 0.01, 0.005))
    rob['two_percent_strains'] = with_(Cprime_GPa=cp(A0, 48, 0.01, 0.02), C44_GPa=c44(A0, 48, 0.01, 0.02))
    rob['three_percent_tetragonal_strains'] = with_(Cprime_GPa=cp(A0, 48, 0.01, 0.03))
    rob['a0_from_zero_pressure'] = with_(a0_A=a0_p0)
    # cutoff: the 60 Ry - 40 Ry change of C' at the same mesh (k40, 0.01 Ry) added to the reference
    cp60 = S.cprime_stress(R[f't_a{A0}_k40_s0.01_d+0.01_e60'], R[f't_a{A0}_k40_s0.01_d-0.01_e60'])
    rob['cutoff_60Ry_correction'] = with_(Cprime_GPa=ref['Cprime_GPa'] + cp60 - cp(A0, 40, 0.01))
    print('C prime k40/0.01: 40 Ry', cp(A0, 40, 0.01), '60 Ry', cp60,
          '; pressure at a0 k40: 40 Ry', c(A0, 40, 0.01)['pressure_GPa'], '60 Ry', R[f'c_a{A0}_k40_s0.01_e60']['pressure_GPa'])
    # C' from a uniaxial strain e3 = +-1 %: stress route (sigma_zz - sigma_xx = 2 C' e3) and C11 (energy) minus B.
    up, um = R[f'u_a{A0}_k48_s0.01_e+0.01'], R[f'u_a{A0}_k48_s0.01_e-0.01']
    c0 = c(A0, 48, 0.01)
    su, sm = S.cauchy(up), S.cauchy(um)
    rob['uniaxial_strain_stress_route'] = with_(Cprime_GPa=float(((su[2][2] - su[0][0]) - (sm[2][2] - sm[0][0])) / (2 * 0.02)))
    c11_e = (up['free_energy_Ry'] + um['free_energy_Ry'] - 2 * c0['free_energy_Ry']) * S.RY_EV / c0['cell_volume_A3'] * S.EVA3_GPA / 0.01 ** 2
    rob['cprime_from_c11_energy_minus_b'] = with_(Cprime_GPa=float(0.75 * (c11_e - b0)))

    # ---------------- ablations (realistic wrong paths; must fail)
    abl, dps = {}, []
    # dp1: the level at which total energies and pressures are converged (k16, 0.02 Ry) used for the moduli.
    a0_16, b0_16 = S.birch_murnaghan([3.39, 3.41, 3.43, 3.45, 3.47],
                                     [c(a, 16, 0.02)['free_energy_Ry'] for a in [3.39, 3.41, 3.43, 3.45, 3.47]])
    abl['total_energy_converged_mesh'] = with_(a0_A=round(a0_16, 4), bulk_modulus_GPa=b0_16,
                                               Cprime_GPa=cp(A0, 16, 0.02), C44_GPa=c44(A0, 16, 0.02),
                                               compressed_pressure_GPa=c(AC, 16, 0.02)['pressure_GPa'],
                                               compressed_Bprime_GPa=cp(AC, 16, 0.02), compressed_B44_GPa=c44(AC, 16, 0.02))
    dps.append(('dp1', 'total_energy_converged_mesh'))
    # dp2: small smearing on coarse meshes (noise; can change the sign of C').
    abl['small_smearing_coarse_mesh'] = with_(Cprime_GPa=cp(A0, 16, 0.005), C44_GPa=c44(A0, 16, 0.005),
                                              compressed_Bprime_GPa=cp(AC, 24, 0.005), compressed_B44_GPa=c44(AC, 24, 0.005))
    dps.append(('dp2', 'small_smearing_coarse_mesh'))
    abl['small_smearing_k24'] = with_(Cprime_GPa=cp(A0, 24, 0.005), compressed_Bprime_GPa=cp(AC, 24, 0.005))
    dps.append(('dp2', 'small_smearing_k24'))
    # dp3: the largest smearing width, k-converged but biased.
    abl['largest_smearing_k_converged'] = with_(Cprime_GPa=cp(A0, 48, 0.02), C44_GPa=c44(A0, 32, 0.02),
                                                compressed_Bprime_GPa=cp(AC, 32, 0.02), compressed_B44_GPa=c44(AC, 32, 0.02))
    dps.append(('dp3', 'largest_smearing_k_converged'))
    # dp4: one-sided strain (single strained cell; the Bain path's third-order term does not cancel).
    tp, tm = tet(A0, 48, 0.01)
    rp, rm = tri(A0, 48, 0.01)
    shear0 = lambda r: S.shear_stress(r) / r['strain_voigt'][3]
    abl['one_sided_plus_strain'] = with_(Cprime_GPa=S.cprime_one_sided(tp), C44_GPa=shear0(rp),
                                         compressed_Bprime_GPa=S.cprime_one_sided(tet(AC, 48, 0.01)[0]),
                                         compressed_B44_GPa=shear0(tri(AC, 48, 0.01)[0]))
    dps.append(('dp4', 'one_sided_plus_strain'))
    abl['one_sided_minus_strain'] = with_(Cprime_GPa=S.cprime_one_sided(tm), C44_GPa=shear0(rm),
                                          compressed_Bprime_GPa=S.cprime_one_sided(tet(AC, 48, 0.01)[1]),
                                          compressed_B44_GPa=shear0(tri(AC, 48, 0.01)[1]))
    dps.append(('dp4', 'one_sided_minus_strain'))
    # dp5: strain about the experimental lattice parameter (78 K) instead of the protocol's a0.
    abl['experimental_lattice'] = with_(a0_A=AX, Cprime_GPa=cp(AX, 48, 0.01), C44_GPa=c44(AX, 48, 0.01))
    dps.append(('dp5', 'experimental_lattice'))
    # dp6: energy curvature of volume-changing strains in the compressed cell, without the -P dV work.
    cc = c(AC, 48, 0.01)
    lin = (R[f'tl_a{AC}_k48_s0.01_d+0.01'], R[f'tl_a{AC}_k48_s0.01_d-0.01'])
    d2 = (lin[0]['free_energy_Ry'] + lin[1]['free_energy_Ry'] - 2 * cc['free_energy_Ry']) * S.RY_EV / cc['cell_volume_A3'] * S.EVA3_GPA
    bp_lin = float(d2 / (2 * 6 * 0.01 ** 2))
    abl['energy_curvature_under_pressure'] = with_(compressed_Bprime_GPa=bp_lin,
                                                   compressed_B44_GPa=S.c44_energy_linear(cc, *tri(AC, 48, 0.01)))
    dps.append(('dp6', 'energy_curvature_under_pressure'))
    pc = cc['pressure_GPa']
    abl['lagrangian_coefficients_under_pressure'] = with_(compressed_Bprime_GPa=ref['compressed_Bprime_GPa'] + pc,
                                                          compressed_B44_GPa=ref['compressed_B44_GPa'] + pc)
    dps.append(('dp6', 'lagrangian_coefficients_under_pressure'))
    # dp7: tensor shear strain read as engineering (or vice versa): energy route gives C44/4.
    abl['tensor_shear_convention'] = with_(C44_GPa=ref['C44_GPa'] / 4, compressed_B44_GPa=ref['compressed_B44_GPa'] / 4)
    dps.append(('dp7', 'tensor_shear_convention'))
    abl['trigonal_single_component'] = with_(C44_GPa=ref['C44_GPa'] * 3, compressed_B44_GPa=ref['compressed_B44_GPa'] * 3)
    dps.append(('dp7', 'trigonal_single_component'))
    # dp8: compressed-cell shear moduli not computed: a0 values reused, or scaled with the bulk modulus.
    abl['compressed_reuses_a0_moduli'] = with_(compressed_Bprime_GPa=ref['Cprime_GPa'], compressed_B44_GPa=ref['C44_GPa'])
    dps.append(('dp8', 'compressed_reuses_a0_moduli'))
    # B at a = 3.30 A from the a0 equation of state, B(P) = B0 + B0' P with B0' from the same fit
    b_c = bulk_modulus_at(eos_a, [c(a, 48, 0.01)['free_energy_Ry'] for a in eos_a], AC)
    print('B at the compressed cell (EOS extrapolation):', b_c)
    abl['compressed_scaled_with_bulk_modulus'] = with_(compressed_Bprime_GPa=ref['Cprime_GPa'] * b_c / b0,
                                                       compressed_B44_GPa=ref['C44_GPa'] * b_c / b0)
    dps.append(('dp8', 'compressed_scaled_with_bulk_modulus'))
    hidden = {'reference': ref, 'tolerances': TOL, 'booleans': {'mechanically_stable_at_a0': True},
              'provenance': {'a0_birch_murnaghan_A': a0_fit, 'a0_zero_pressure_A': a0_p0,
                             'Cprime_set_GPa': cp_set, 'C44_set_GPa': c44_set, 'compressed_pressure_set_GPa': pc_set,
                             'compressed_Bprime_set_GPa': bpc_set, 'compressed_B44_set_GPa': b44c_set,
                             'runs': 'benchmarks/construction/fresh-li-bcc-soft-shear/reference_results.json'}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [v for v in ref.values() if not isinstance(v, bool)] + [a0_fit, a0_p0] + \
        [v for d in (cp_set, c44_set, pc_set, bpc_set, b44c_set) for v in d.values()]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob,
                                                                  'ablation_decision_points': dps}, indent=1) + '\n')

    # ---------------- calibration table
    numeric = [f for f in FIELDS if f not in BOOLEANS]

    def err(sub, f):
        e = sub[f] - ref[f]
        return (e / abs(ref[f]) if f in TOL['relative'] else e), sub[f]

    table = {}
    for f in numeric:
        tol = TOL['relative'].get(f, TOL['absolute'].get(f))
        row = {'reference': ref[f], 'tolerance': ('relative ' if f in TOL['relative'] else 'absolute ') + str(tol),
               'robustness': {}, 'ablations_failing_this_field': {}}
        for name, sub in rob.items():
            e, v = err(sub, f)
            if abs(e) > 1e-12:
                row['robustness'][name] = {'value': v, 'error': e}
        for name, sub in abl.items():
            e, v = err(sub, f)
            if abs(e) > tol:
                row['ablations_failing_this_field'][name] = {'value': v, 'error': e}
        worst = max([abs(x['error']) for x in row['robustness'].values()], default=0.0)
        nearest = min([abs(x['error']) for x in row['ablations_failing_this_field'].values()], default=None)
        row['max_robust_error'], row['nearest_failing_ablation_error'] = worst, nearest
        table[f] = row
        flag = '' if worst <= tol else '  <-- ROBUSTNESS FAILS'
        print(f"{f:26s} ref={ref[f]:.5g} tol={tol} max_robust={worst:.4g} nearest_abl={nearest}{flag}")
    table['mechanically_stable_at_a0'] = {'reference': True, 'tolerance': 'exact',
                                          'ablations_failing_this_field': [n for n, s in abl.items()
                                                                           if s['mechanically_stable_at_a0'] is not True]}
    for name, sub in abl.items():
        failing = [f for f in numeric if abs(err(sub, f)[0]) > TOL['relative'].get(f, TOL['absolute'].get(f))]
        failing += [b for b in BOOLEANS if sub[b] is not True]
        print(f'ablation {name:40s} fails {failing}')
    print(json.dumps({'reference': ref, 'sets': hidden['provenance']}, indent=1, default=float))
    (HERE/'calibration.json').write_text(json.dumps(table, indent=1) + '\n')


if __name__ == '__main__':
    main()
