"""Derive the hidden references, ablation values and calibration for fresh-nio-exchange-dftu.

    python3 -B benchmarks/construction/fresh-nio-exchange-dftu/make_task.py

Reads reference_runs.jsonl (lab results of reference_plan.PLAN, made by run_reference.py through
task_lab.run) and applies the estimators of the task's own solution/solve.py. Writes
environment/lab/private/hidden.json, private/hidden_values.json, solution/ablation_values.json,
calibration.json and the "calibration" block of design.json, and prints the calibration table.
"""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-nio-exchange-dftu'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
ROWS = {r['tag']: r for r in map(json.loads, (HERE/'reference_runs.jsonl').read_text().splitlines())}
TOL = {'relative': {'J1_meV': 0.03, 'J2_meV': 0.02, 'neel_temperature_mean_field_K': 0.02,
                    'pbe_fm_minus_afm2_meV_per_Ni': 0.02},
       'absolute': {'fm_gap_closing_U_eV': 0.1}}
GRADED_U = ('4.0', '8.0')


def res(tag):
    row = ROWS[tag]
    if 'result' not in row:
        raise SystemExit(f'reference calculation {tag} failed: {row.get("error")}')
    return row['result']


def E(tag):
    return S.energy(res(tag))


def pair_rh(U, k, v=''):
    return res(f'U{U}_rh_fm_k{k}{v}'), res(f'U{U}_rh_af2_k{k}{v}')


def pair_te(U, k, v=''):
    return res(f'U{U}_te_fm_k{k}{v}'), res(f'U{U}_te_af1_k{k}{v}')


def J(U, rh='777', te='10108', v=''):
    return S.exchange(*pair_rh(U, rh, v), *pair_te(U, te, v))


def fields(j, keys=GRADED_U):
    """J1, J2 and T_N fields from {U: (J1, J2)}."""
    return {'J1_meV': {u: j[u][0] for u in keys}, 'J2_meV': {u: j[u][1] for u in keys},
            'neel_temperature_mean_field_K': {u: S.neel_mean_field(j[u][1]) for u in keys}}


def gaps(k=12):
    return [(float(t.split('_')[0][1:]), r['result']['band_gap_eV']) for t, r in ROWS.items()
            if t.endswith(f'_pr_fm_k{k}{k}{k}') and 'result' in r]


def main():
    j_ref = {u: J(u) for u in GRADED_U}
    de0_ref = E('U0.0_pr_fm_k161616') - E('U0.0_rh_af2_k777')
    uc_ref = S.gap_closing(gaps(12))
    ref = {**fields(j_ref), 'pbe_fm_minus_afm2_meV_per_Ni': de0_ref, 'fm_gap_closing_U_eV': uc_ref}
    # Sanity of the protocol: same states from other starts, and gapped at the graded U.
    for u in GRADED_U:
        assert abs(E(f'U{u}_rh_fm_k444_m1') - E(f'U{u}_rh_fm_k444')) < 0.05
        assert abs(E(f'U{u}_rh_af2_k444_m1') - E(f'U{u}_rh_af2_k444')) < 0.05
        assert abs(E(f'U{u}_te_af1_k886_m01') - E(f'U{u}_te_af1_k886')) < 0.05
        assert all(res(t)['band_gap_eV'] > 0.5 for t in (f'U{u}_rh_fm_k777', f'U{u}_te_af1_k10108'))

    # ---- ablations: realistic wrong paths, computed from the reference calculations ----
    A, order = {}, []

    def abl(dp, name, value):
        A[name] = value
        order.append((dp, name))

    # dp1: FM and AFM-II only (J1 undetermined, set to zero)
    abl('dp1', 'fm_afm2_only_J1_zero', fields({u: (0.0, (E(f'U{u}_rh_fm_k777') - E(f'U{u}_rh_af2_k777')) / 6)
                                               for u in GRADED_U}))
    # dp2: one FM reference (rhombohedral 4^3, where FM - AFM-II is converged) also used for AFM-I
    j = {}
    for u in GRADED_U:
        fm = E(f'U{u}_rh_fm_k444')
        j1 = (fm - E(f'U{u}_te_af1_k10108')) / 8
        j[u] = (j1, (fm - E(f'U{u}_rh_af2_k444')) / 6 - j1)
    abl('dp2', 'cross_cell_fm_reference', fields(j))
    # dp3: coarse rhombohedral mesh for FM - AFM-II
    abl('dp3', 'coarse_rhombohedral_mesh_k3', fields({u: J(u, rh='333') for u in GRADED_U}))
    # dp4: tetragonal mesh that looks converged at large U reused at small U (6x6x4)
    abl('dp4', 'tetragonal_mesh_664', fields({u: J(u, te='664') for u in GRADED_U}))
    abl('dp4', 'tetragonal_mesh_444', fields({u: J(u, te='444') for u in GRADED_U}))
    # dp5: energies per cell (2 Ni) instead of per Ni
    abl('dp5', 'energies_per_cell', {f: {u: 2 * v for u, v in ref[f].items()} for f in fields(j_ref)})
    # dp6: Neel temperature conventions
    abl('dp6', 'classical_S2_neel', {'neel_temperature_mean_field_K': {u: ref['neel_temperature_mean_field_K'][u] / 2
                                                                       for u in GRADED_U}})
    abl('dp6', 'curie_weiss_as_neel', {'neel_temperature_mean_field_K': {
        u: S.S * (S.S + 1) / 3 * abs(12 * j_ref[u][0] + 6 * j_ref[u][1]) / S.KB_MEV for u in GRADED_U}})
    # dp7: metal-style smearing at U = 4 eV (FM gap 0.77 eV): Fermi-Dirac 0.02 Ry or Gaussian 0.05 Ry close the FM gap
    for name, v in (('fermi_dirac_0p02_at_U4', '_fd0.02'), ('gaussian_0p05_at_U4', '_gaussian0.05')):
        abl('dp7', name, fields({'4.0': S.exchange(*pair_rh('4.0', '666', v), *pair_te('4.0', '886', v))}, keys=('4.0',)))
    # dp8: plain PBE metallic FM: smearing (Fermi-Dirac; Gaussian free energy F without the zero-smearing estimate)
    af2 = E('U0.0_rh_af2_k777')
    abl('dp8', 'pbe_fermi_dirac_0p02', {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k121212_fd0.02') - af2})
    g = res('U0.0_pr_fm_k121212_gaussian0.02')
    abl('dp8', 'pbe_gaussian_0p02_free_energy', {'pbe_fm_minus_afm2_meV_per_Ni':
                                                 g['free_energy_Ry_per_cell'] * S.RY_MEV - af2})
    # dp9: plain PBE with meshes adequate for the gapped U > 0 states (8^3 primitive FM, 4^3 rhombohedral AFM-II)
    abl('dp9', 'pbe_coarse_meshes', {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k888') - E('U0.0_rh_af2_k444')})
    abl('dp9', 'pbe_primitive_fm_k8', {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k888') - af2})
    # dp10: small starting moment: metastable low-moment FM solution
    abl('dp10', 'pbe_small_starting_moment', {'pbe_fm_minus_afm2_meV_per_Ni':
                                              E('U0.0_rh_fm_k666_m01') - E('U0.0_rh_af2_k666')})
    # dp11: FM gap closing from the graded U values only, or from a coarse U bracket
    g4, g8 = res('U4.0_pr_fm_k121212')['band_gap_eV'], res('U8.0_pr_fm_k121212')['band_gap_eV']
    abl('dp11', 'gap_closing_linear_from_U4_U8', {'fm_gap_closing_U_eV': 4.0 - g4 * 4.0 / (g8 - g4)})
    abl('dp11', 'gap_closing_coarse_bracket', {'fm_gap_closing_U_eV': 3.0})
    # dp12: shifted (kshift = 1) meshes miss the FM band extrema, so the gap on the mesh stays open below U_c
    shifted = [(u, res(f'U{u}_pr_fm_k121212_s1')['band_gap_eV']) for u in ('2.0', '2.2', '2.5')]
    abl('dp12', 'gap_closing_shifted_mesh', {'fm_gap_closing_U_eV': S.gap_closing([(float(u), g) for u, g in shifted])})

    # ---- robustness: defensible alternatives that must pass ----
    Rb = {}
    Rb['oracle_meshes_rh666_te886'] = fields({u: J(u, rh='666', te='886') for u in GRADED_U})
    Rb['rh444_te888'] = fields({u: J(u, rh='444', te='888') for u in GRADED_U})
    j = {}
    for u in GRADED_U:   # converged primitive-cell FM as the common reference for both differences
        fm = E(f'U{u}_pr_fm_k121212')
        j1 = (fm - E(f'U{u}_te_af1_k10108')) / 8
        j[u] = (j1, (fm - E(f'U{u}_rh_af2_k777')) / 6 - j1)
    Rb['converged_primitive_fm_reference'] = fields(j)
    Rb['cutoff_75_Ry'] = fields({u: J(u, rh='444', te='886', v='_e75') for u in GRADED_U})
    Rb['cutoff_40_Ry'] = fields({u: J(u, rh='444', te='886', v='_e40') for u in GRADED_U})
    Rb['pbe_k12_gauss0p002'] = {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k121212_g002') - af2}
    Rb['pbe_k16_gauss0p002'] = {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k161616_g002') - af2}
    Rb['pbe_k12_mv0p005'] = {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k121212_mv0.005') - af2}
    Rb['pbe_afm2_k6'] = {'pbe_fm_minus_afm2_meV_per_Ni': E('U0.0_pr_fm_k161616') - E('U0.0_rh_af2_k666')}
    Rb['gap_closing_k16_points'] = {'fm_gap_closing_U_eV': S.gap_closing(gaps(16) + [(2.5, dict(gaps(12))[2.5])])}
    Rb['gap_closing_line_from_3_and_4'] = {'fm_gap_closing_U_eV': S.gap_closing(
        [(u, g) for u, g in gaps(12) if u in (3.0, 3.5)])}
    Rb['gap_closing_oracle_points'] = {'fm_gap_closing_U_eV': S.gap_closing(
        [(u, g) for u, g in gaps(12) if u in (2.0, 2.2, 2.5)] + [(4.0, res('U4.0_rh_fm_k666')['band_gap_eV'])])}

    hidden = {'reference': ref, 'tolerances': TOL,
              'provenance': 'derived from reference_runs.jsonl (task_lab.run, QE 6.7, 2 MPI ranks) by make_task.py'}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [v for f in ('J1_meV', 'J2_meV', 'neel_temperature_mean_field_K') for v in ref[f].values()]
    values += [ref['pbe_fm_minus_afm2_meV_per_Ni'], ref['fm_gap_closing_U_eV']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablation_order': order, 'ablations': A,
                                                                  'robustness': Rb}, indent=1) + '\n')

    # ---- calibration table ----
    def flat(d):
        out = {}
        for f, v in d.items():
            if isinstance(v, dict):
                out.update({f'{f}[{k}]': x for k, x in v.items()})
            else:
                out[f] = v
        return out
    fref = flat(ref)
    cal = {}
    for key, value in fref.items():
        field = key.split('[')[0]
        rel = field in TOL['relative']
        err = (lambda x: x / value - 1) if rel else (lambda x: x - value)
        cal[key] = {'reference': value,
                    'tolerance': f"{'+-%g%%' % (100 * TOL['relative'][field]) if rel else '+-%g' % TOL['absolute'][field]}",
                    'robustness': {n: round(err(flat(c)[key]), 5) for n, c in Rb.items() if key in flat(c)},
                    'ablations': {n: round(err(flat(c)[key]), 5) for n, c in A.items() if key in flat(c)}}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1) + '\n')
    design = json.loads((TASK/'design.json').read_text())
    design['calibration'] = cal
    (TASK/'design.json').write_text(json.dumps(design, indent=1, ensure_ascii=False) + '\n')
    for key, c in cal.items():
        worst_rob = max([abs(v) for v in c['robustness'].values()], default=0)
        near_abl = min([abs(v) for v in c['ablations'].values()], default=None)
        print(f"{key:45s} ref={c['reference']:.5g} tol={c['tolerance']:>8s} max|robust|={worst_rob:.4f} "
              f"nearest|ablation|={near_abl}")
        for n, v in c['ablations'].items():
            print(f'      abl {n}: {v:+.4f}')
        for n, v in c['robustness'].items():
            print(f'      rob {n}: {v:+.4f}')


if __name__ == '__main__':
    main()
