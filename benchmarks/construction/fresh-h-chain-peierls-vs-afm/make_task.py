"""Derive the hidden references, ablation values and calibration for fresh-h-chain-peierls-vs-afm.

    python3 -B benchmarks/construction/fresh-h-chain-peierls-vs-afm/make_task.py

Reads reference_runs.jsonl (lab results of reference_plan.PLAN, made by run_reference.py through
task_lab.run) and applies the estimators of the task's own solution/solve.py. Writes
environment/lab/private/hidden.json, private/hidden_values.json, solution/ablation_values.json and the
"calibration" block of design.json, and prints the calibration table.
"""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-h-chain-peierls-vs-afm'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
ROWS = {r['tag']: r for r in map(json.loads, (HERE/'reference_runs.jsonl').read_text().splitlines())}
TOL = {'absolute': {'delta_eq': 0.01, 'dimerization_energy_meV_per_atom': 10.0, 'critical_dimerization': 0.015,
                    'uniform_chain_band_gap_eV': 0.04},
       'relative': {'uniform_magnetic_energy_gain_meV_per_atom': 0.10}}


def res(tag):
    row = ROWS[tag]
    if 'result' not in row:
        raise SystemExit(f'reference calculation {tag} failed: {row.get("error")}')
    return row['result']


def E(tag):
    return S.energy(res(tag))


def equilibrium(R, deltas):
    return S.equilibrium([(d, res(f'R{R}_d{d}_nm')) for d in deltas], R)


def boundary_points(R):
    pts = [(0.0, res(f'R{R}_d0_afm'))]
    pts += [(float(t.split('_d')[1].split('_')[0]), r['result']) for t, r in ROWS.items()
            if t.startswith(f'R{R}_d') and t.endswith('_afm') and '_d0_' not in t and 'result' in r]
    return sorted(pts, key=lambda p: p[0])


def label(R, dimer_tags, afm_dimer_tag):
    uniform = min(E(f'R{R}_d0_nm'), E(f'R{R}_d0_afm'))
    dimer = min(E(t) for t in dimer_tags)
    if dimer < uniform:
        return 'dimerized_magnetic' if S.is_magnetic(res(afm_dimer_tag)) else 'dimerized_nonmagnetic'
    return 'uniform_antiferromagnetic' if E(f'R{R}_d0_afm') < E(f'R{R}_d0_nm') else 'uniform_nonmagnetic'


def per_key(fn, keys):
    return {k: fn(k) for k in keys}


def main():
    eq = {'1.25': equilibrium(1.25, (0.39, 0.4, 0.41)), '1.75': equilibrium(1.75, (0.56, 0.57, 0.58))}
    uniform = {k: min(E(f'R{k}_d0_nm'), E(f'R{k}_d0_afm')) for k in eq}
    crit = {k: S.critical_delta(boundary_points(float(k))) for k in ('1.75', '2.0')}
    ref = {
        'delta_eq': {k: eq[k][0] for k in eq},
        'dimerization_energy_meV_per_atom': {k: uniform[k] - eq[k][1] for k in eq},
        'uniform_magnetic_energy_gain_meV_per_atom': {k: E(f'R{k}_d0_nm_mv002_k96') - E(f'R{k}_d0_afm_mv002_k96')
                                                      for k in ('1.75', '2.0')},
        'critical_dimerization': crit,
        'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_afm_mv002_k192')['band_gap_eV']},
        'ground_state': {'1.0': label(1.0, [f'R1.0_d{d}_nm' for d in (0.2, 0.25, 0.3)], 'R1.0_d0.25_afm'),
                         '2.0': label(2.0, [f'R2.0_d{d}_nm' for d in (0.6, 0.625, 0.65)], 'R2.0_d0.625_afm')}}
    # The lowest uniform state at 1.25 A must be the (weak) SDW in the converged limit.
    assert E('R1.25_d0_afm_mv002_k192') < E('R1.25_d0_nm_mv002_k192') and S.is_magnetic(res('R1.25_d0_afm_mv002_k192'))

    # ---- ablations: realistic wrong paths, computed from the reference calculations ----
    A = {}
    A['ferromagnetic_start_only'] = {
        'uniform_magnetic_energy_gain_meV_per_atom': {k: E(f'R{k}_d0_nm') - E(f'R{k}_d0_fm') for k in ('1.75', '2.0')},
        'critical_dimerization': {k: 0.0 for k in ('1.75', '2.0')},
        'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_fm')['band_gap_eV']},
        'dimerization_energy_meV_per_atom': {k: min(E(f'R{k}_d0_nm'), E(f'R{k}_d0_fm')) - eq[k][1] for k in eq}}
    A['one_atom_cell_spin_polarized'] = {
        'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_fm_1atom')['band_gap_eV']},
        'uniform_magnetic_energy_gain_meV_per_atom': {'1.75': E('R1.75_d0_nm') - E('R1.75_d0_fm_1atom'),
                                                      '2.0': E('R2.0_d0_nm') - E('R2.0_d0_fm')}}
    A['large_smearing_mv002'] = {'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_afm_mv02')['band_gap_eV']}}
    A['fermi_dirac_0005'] = {'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_afm_fd005_k96')['band_gap_eV']}}
    A['coarse_k_sdw_gap'] = {'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_afm')['band_gap_eV']}}
    grid = {'1.75': (0.1, 0.2), '2.0': (0.2, 0.3)}       # last magnetic / first nonmagnetic on a 0.1 grid
    A['coarse_boundary_grid'] = {'critical_dimerization': {k: 0.5 * sum(grid[k]) for k in grid}}

    def linear_m(k):
        (d1, m1), (d2, m2) = [(d, S.site_moment(res(f'R{k}_d{d}_afm' if d else f'R{k}_d0_afm')))
                              for d in (grid[k][0] - 0.1, grid[k][0])]
        return d2 + (d2 - d1) * m2 / (m1 - m2)
    A['linear_moment_extrapolation'] = {'critical_dimerization': per_key(linear_m, grid)}
    A['nonmagnetic_uniform_reference'] = {'dimerization_energy_meV_per_atom':
                                          {k: E(f'R{k}_d0_nm') - eq[k][1] for k in eq}}
    A['afm_wins_at_large_spacing'] = {'ground_state': {'1.0': 'dimerized_nonmagnetic', '2.0': 'uniform_antiferromagnetic'}}
    gmin = {'1.25': 0.4, '1.75': 0.6}
    A['grid_minimum_dimerization'] = {'delta_eq': dict(gmin),
                                      'dimerization_energy_meV_per_atom':
                                      {k: uniform[k] - E(f'R{k}_d{gmin[k]}_nm') for k in gmin}}
    A['energies_per_cell'] = {f: {k: 2 * v for k, v in ref[f].items()}
                              for f in ('dimerization_energy_meV_per_atom', 'uniform_magnetic_energy_gain_meV_per_atom')}

    # ---- robustness: defensible alternatives that must pass ----
    Rb = {}
    Rb['standard_settings_magnetic_gain'] = {'uniform_magnetic_energy_gain_meV_per_atom':
                                             {k: E(f'R{k}_d0_nm') - E(f'R{k}_d0_afm') for k in ('1.75', '2.0')}}

    def bisection(k):
        pts = boundary_points(float(k))
        last = max(d for d, r in pts if S.is_magnetic(r))
        first = min(d for d, r in pts if not S.is_magnetic(r) and d > last)
        return 0.5 * (last + first)
    Rb['boundary_bisection_midpoint'] = {'critical_dimerization': per_key(bisection, ('1.75', '2.0'))}
    Rb['boundary_converged_settings'] = {'critical_dimerization': {
        '1.75': S.critical_delta([(0.0, res('R1.75_d0_afm_mv002_k96')), (0.1, res('R1.75_d0.1_afm')),
                                  (0.11, res('R1.75_d0.11_afm_mv002_k96')), (0.12, res('R1.75_d0.12_afm'))]),
        '2.0': S.critical_delta([(0.0, res('R2.0_d0_afm_mv002_k96')), (0.19, res('R2.0_d0.19_afm')),
                                 (0.2, res('R2.0_d0.2_afm_mv002_k96')), (0.21, res('R2.0_d0.21_afm'))])}}
    Rb['gap_k96_mv0002'] = {'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_afm_mv002_k96')['band_gap_eV']}}
    Rb['gap_gaussian_k96'] = {'uniform_chain_band_gap_eV': {'1.25': res('R1.25_d0_afm_g005_k96')['band_gap_eV']}}
    Rb['gap_4atom_cell_k96'] = {'uniform_chain_band_gap_eV': {'1.25': res('R1.25_4atom_afm_mv002_k96')['band_gap_eV']}}

    def vertex(R, deltas):
        a, b, c = S.parabola([(d, E(f'R{R}_d{d}_nm')) for d in deltas])
        return -b / (2 * c), a - b * b / (4 * c)
    vx = {'1.25': vertex(1.25, (0.39, 0.4, 0.41)), '1.75': vertex(1.75, (0.56, 0.57, 0.58))}
    Rb['energy_parabola_minimum'] = {'delta_eq': {k: vx[k][0] for k in vx},
                                     'dimerization_energy_meV_per_atom': {k: uniform[k] - vx[k][1] for k in vx}}
    Rb['cutoff_60_Ry'] = {'dimerization_energy_meV_per_atom': {
        '1.25': ref['dimerization_energy_meV_per_atom']['1.25'],
        '1.75': E('R1.75_d0_afm_ec60') - E('R1.75_d0.57_nm_ec60')}}
    Rb['converged_uniform_reference'] = {'dimerization_energy_meV_per_atom': {
        '1.25': min(E('R1.25_d0_nm_mv002_k192'), E('R1.25_d0_afm_mv002_k192')) - eq['1.25'][1],
        '1.75': min(E('R1.75_d0_nm_mv002_k96'), E('R1.75_d0_afm_mv002_k96')) - eq['1.75'][1]}}

    hidden = {'reference': ref, 'tolerances': TOL,
              'provenance': 'derived from reference_runs.jsonl (task_lab.run, QE 6.7, 2 MPI ranks) by make_task.py'}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [v for f, d in ref.items() if f != 'ground_state' for v in d.values()]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': Rb}, indent=1) + '\n')

    # ---- calibration table ----
    calibration = {}
    for field, keys in ((f, list(d)) for f, d in ref.items() if f != 'ground_state'):
        for k in keys:
            tol = TOL['relative'].get(field)
            entry = {'reference': ref[field][k],
                     'tolerance': f'within {tol:.0%}' if tol else f"+-{TOL['absolute'][field]}",
                     'robustness': {n: v[field][k] for n, v in Rb.items() if field in v and k in v[field]},
                     'ablations': {n: v[field][k] for n, v in A.items() if field in v and k in v[field]}}
            calibration[f'{field}[{k}]'] = entry
    calibration['ground_state'] = {'reference': ref['ground_state'],
                                   'ablations': {'afm_wins_at_large_spacing': A['afm_wins_at_large_spacing']['ground_state']}}
    design_path = TASK/'design.json'
    design = json.loads(design_path.read_text())
    design['calibration'] = json.loads(json.dumps(calibration, default=float))
    design_path.write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps(calibration, indent=1, default=float))


if __name__ == '__main__':
    main()
