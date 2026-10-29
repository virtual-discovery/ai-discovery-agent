"""Write the catalog, lab data and private targets of mpea-vacancy-diffusion-kmc.

    python3 -B make_task.py          (needs pool.json from extract.py and reference.json from reference.py)
"""
import collections
import json
import numpy as np
import kmclib as K
import common as C

TASK = K.TASK
BUDGET = 60
HOPS_PER_PURCHASE = 20
TOLERANCES = {'barrier_rmse_eV': 0.23, 'log10_eta': 0.3, 'log10_tracer_ratio_B': 0.3, 'log10_eta_thermo_ideal_B': 0.3,
              'log10_eta_kra_simple_B': 0.3, 'activation_energy_eV': 0.2}
LIMITS = {'temperature_K': [1000, 3500], 'cell': [4, 14], 'steps': [100, 200000], 'trajectories': [1, 400],
          'equilibration_sweeps': [0, 500], 'checkpoints': [0, 10], 'work_per_unit': 2000000, 'sweep_site_per_step': 8,
          'max_cost': 8}
OUTPUT_DEFS = {
    'time_ps': 'kMC time of the trajectory in ps (attempt frequency nu = 1 THz)',
    'vacancy_r2_a2': 'squared net displacement of the vacancy over the trajectory, in units of a^2 (a = cubic lattice parameter)',
    'sum_r2_a2': 'sum over all atoms of each element of the squared net displacement, a^2',
    'hops': 'number of vacancy exchanges performed by atoms of each element',
    'mean_encountered_barrier_eV': 'mean barrier of all 8 candidate hops, averaged over steps',
    'mean_accepted_barrier_eV': 'mean barrier of the hops performed',
    'equilibration_acceptance': 'acceptance ratio of the equilibration swaps',
    'checkpoints': 'optional intermediate (step, time_ps, vacancy_r2_a2, sum_r2_a2) records',
    'nn_pair_counts_after_equilibration': 'ordered nearest-neighbour pair counts between elements, summed over trajectories, '
                                          'right after equilibration (before the kMC hops)'}


def hop_id(n):
    return 'H%04d' % n


def composition(h):
    c = collections.Counter(x for x in h['occupation_initial'] if x != 'Va')
    return {e: c[e] for e in K.EL if c[e]}


def category(h):
    c = composition(h)
    return 'pure' if len(c) == 1 else ('dilute' if min(c.values()) == 1 and len(c) == 2 else 'alloy')


def main():
    hidden = set(C.hidden_indices().tolist())
    lab_hops, avail, predict, hidden_barriers = {}, [], [], {}
    for n, h in enumerate(K.POOL):
        rec = {'id': hop_id(n), 'occupation_initial': h['occupation_initial'],
               'vacancy_site_initial': h['vacancy_site_initial'], 'vacancy_site_final': h['vacancy_site_final'],
               'hopping_species': h['hopping_species']}
        if n in hidden:
            predict.append({**rec, 'composition': composition(h)})
            e = h['energy_eV']
            hidden_barriers[hop_id(n)] = e['transition_state'] - e['end_state_1']
        else:
            lab_hops[hop_id(n)] = {**rec, 'energy_eV': h['energy_eV'], 'archived_cell_edge_A': h['archived_cell_edge_A']}
            avail.append({'id': hop_id(n), 'category': category(h), 'hopping_species': h['hopping_species'],
                          'composition': composition(h)})
    orbits = [{'index': i + 1, 'distances_a': list(key), 'n_sites': len(off),
               'offsets_half_a': K.half_units(off).tolist()} for i, (key, off) in enumerate(K.ORBITS)]
    shells = [{'index': s + 1, 'distance_a': float(K.SHELL_R[s]), 'n_sites': len(off),
               'offsets_half_a': K.half_units(off).tolist()} for s, off in enumerate(K.SHELLS)]
    data = {'hops': lab_hops, 'orbits': [o['offsets_half_a'] for o in orbits], 'shells': [s['offsets_half_a'] for s in shells],
            'kmc_limits': LIMITS, 'hops_per_purchase': HOPS_PER_PURCHASE, 'kmc_output_definitions': OUTPUT_DEFS}
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')))
    catalog = {
        'material': 'bcc V-Cr-Nb-Mo-Ta-W alloys; single-vacancy nearest-neighbour hops released by Lee and Natarajan (2026)',
        'elements': K.EL,
        'supercell': {
            'description': '4 x 4 x 4 conventional bcc supercell, 128 sites (127 atoms and one vacancy) in every archived calculation',
            'site_index': 'site ((i*4 + j)*4 + k)*2 + b sits at (i + b/2, j + b/2, k + b/2) a, with i, j, k in 0-3 and b in {0, 1}',
            'occupation': 'occupation_initial[site] is the element on each site in end state 1; "Va" marks the vacancy'},
        'hop_records': {
            'vacancy_site_initial': 'vacancy site in end state 1',
            'vacancy_site_final': 'site of the hopping atom in end state 1 = vacancy site in end state 2',
            'energy_eV': 'DFT total energies (eV) of the whole cell: end_state_1, end_state_2 (hopping atom moved into the '
                         'vacancy) and transition_state; PBE single points on MACE-relaxed geometries and MACE CI-NEB saddles',
            'archived_cell_edge_A': 'cell edge of the stored structures, as archived'},
        'hops_available': avail,
        'hops_to_predict': predict,
        'compositions': {c: C.COMPOSITIONS[c] for c in 'ABC'},
        'temperature_K': C.T_ETA,
        'arrhenius_temperatures_K': list(C.T_ARRHENIUS),
        'purchase': {'hops_per_purchase': HOPS_PER_PURCHASE, 'hops_cost': 1,
                     'kmc_cost': 'max(1, ceil(trajectories * (steps + equilibration_sweeps * N / 8) / 2e6)), N = 2 cell^3 sites; '
                                 'at most 8 units per run',
                     'budget': BUDGET},
        'kmc_service': {
            'description': 'Rejection-free kinetic Monte Carlo of one vacancy in a periodic cell x cell x cell conventional '
                           'bcc supercell (N = 2 cell^3 sites, N - 1 atoms). Each trajectory starts from a random alloy '
                           'with the requested composition (atom counts by largest-remainder rounding of N - 1), optionally '
                           'equilibrated by canonical Monte Carlo (swaps of any two different occupants, the vacancy '
                           'included, Metropolis with the model end-state energies; one sweep = N attempts), then runs '
                           '`steps` vacancy hops with rates nu exp(-E_b / k_B T), nu = 1 THz.',
            'barrier': 'E_b = KRA + (E_after - E_before) / 2 for the exchange of the vacancy with the neighbouring atom X; '
                       'KRA = K0[X] + sum over orbits o of mean over the sites of o of W[X][o][element on site] + '
                       'sum_ij Q[X][i][j] c_i c_j, where c is the element composition of orbits 1-3 together (14 sites). '
                       'End-state energy E = sum over shells s and over bonds of V[s][a][b]; with end_state null all end '
                       'states have equal energy.',
            'orbits': orbits,
            'orbit_convention': 'orbit sites for a hop from the vacancy at r_v to the atom at r_v + (sx, sy, sz) a/2 '
                                '(s = +-1) are r_v + (sx*ox, sy*oy, sz*oz) a/2 for the listed offsets (ox, oy, oz); '
                                'orbits are classes of sites with the same sorted distances to the two hop sites',
            'pair_shells': shells,
            'model_format': {'kra': {'orbits': 'number n of leading orbits used (0-12)', 'K0': '{element: eV}',
                                     'W': '{element: n rows of 6 values (eV) in element order}',
                                     'Q': 'optional {element: 6 x 6 matrix (eV)}, needs orbits >= 3'},
                             'end_state': 'null, or {"shells": 1-4, "V": shells x 7 x 7 symmetric matrices (eV), order '
                                          'V Cr Nb Mo Ta W Va}'},
            'request': {'kind': 'kmc', 'composition': {'Nb': 0.5, 'Ta': 0.25, 'W': 0.25}, 'temperature_K': 2000,
                        'cell': 10, 'steps': 20000, 'trajectories': 40, 'equilibration_sweeps': 50, 'checkpoints': 0,
                        'seed': 1, 'model': {'kra': {'orbits': 0, 'K0': {e: 1.0 for e in K.EL}}, 'end_state': None}},
            'limits': {k: v for k, v in LIMITS.items() if k not in ('work_per_unit', 'sweep_site_per_step')},
            'outputs': OUTPUT_DEFS}}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    ref = json.loads((K.HERE/'reference.json').read_text())['targets']
    targets = {k: ref[k] for k in ('log10_eta_vacancy', 'fastest_tracer_B', 'log10_tracer_ratio_B', 'activation_energy_C_eV',
                                   'log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B')}
    private = {'hidden_barriers_eV': hidden_barriers, 'targets': targets, 'tolerances': TOLERANCES,
               'reference_info': {'slowest_tracer_B': ref['slowest_tracer_B_info'], 'fastest_tracer_C': ref['fastest_tracer_C']}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(private, indent=1) + '\n')
    values = list(hidden_barriers.values()) + list(targets['log10_eta_vacancy'].values()) + [
        targets[k] for k in ('log10_tracer_ratio_B', 'activation_energy_C_eV', 'log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B')]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'environment/lab/public.json').write_text(json.dumps({
        'task_id': 'mpea-vacancy-diffusion-kmc', 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
        'requests': {'hops': {'kind': 'hops', 'ids': ['<up to 20 ids from hops_available>'], 'cost': 1},
                     'kmc': {'kind': 'kmc', 'cost': 'see catalog purchase.kmc_cost'}},
        'catalog': '/app/catalog.json'}, indent=1) + '\n')
    (TASK/'tests/public.json').write_text(json.dumps({'task_id': 'mpea-vacancy-diffusion-kmc', 'solver_budget': BUDGET,
                                                      'min_report_chars': 400}, indent=1) + '\n')
    print('available', len(avail), 'to predict', len(predict), 'targets', json.dumps(targets))


if __name__ == '__main__':
    main()
