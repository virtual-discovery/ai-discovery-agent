"""Write the staged, lab, private and ablation files of ti-codoping-exchange-tc.

    python3 -B make_task.py            (after fetch.py; reads reference.json and ablations.json from this folder)

Inputs: the archived Jij_table.csv, Bi2Te3.cif and Tcs_Bi2Te3.txt (cache), reference.json (reference Monte Carlo
ordering temperatures, written by reference.py) and ablations.json (wrong-path values, written by reference.py).
"""
import csv
import json
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/ti-codoping-exchange-tc'
BUDGET = 60
COMPOSITIONS = {
    'V5Cr5': {'V': 0.05, 'Cr': 0.05},
    'Cr5Mo5': {'Cr': 0.05, 'Mo': 0.05},
    'Cr5Mn5': {'Cr': 0.05, 'Mn': 0.05},
    'V2Cr2': {'V': 0.02, 'Cr': 0.02},
    'Cr5Co5': {'Cr': 0.05, 'Co': 0.05},
}
CANDIDATES = ['V', 'Nb', 'Mo', 'Mn']
MC = {'cost_per_16_samples': {'2': 1, '3': 1, '4': 2, '5': 4, '6': 8},
      'max_samples': {'2': 32, '3': 32, '4': 32, '5': 32, '6': 16},
      'max_temperatures': 8, 'equilibration_sweeps': 1000, 'measurement_sweeps': 3000, 'overrelaxation': 2,
      'threads': 4}
COLUMNS = {
    'i_shell': 'index of the impurity-dimer configuration (neighbour shell), 0-16',
    'R': 'impurity-impurity distance (Angstrom)',
    'Rx, Ry, Rz': 'one connecting vector of the shell (Angstrom), as extracted by the authors',
    'ilayer2': 'Bi layer of impurity 2 (impurity 1 always substitutes Bi1): 3 = Bi1 sublattice, 4 = Bi2 sublattice',
    'Zimp1, Zimp2': 'atomic numbers of impurity 1 and impurity 2 (Zimp1 <= Zimp2)',
    'num_neighbors': 'number of equivalent Bi sites in this shell',
    'm1, m2': 'spin moments (muB) of impurity 1 and impurity 2 read from the dimer calculation',
    'm01, m02': 'spin moments (muB) of impurity 1 and impurity 2 as isolated single impurities',
    'J': 'isotropic exchange constant J_ij (meV) of the dimer',
    'J_data_uuid': 'AiiDA uuid of the JijData node (provenance)',
}


def read_rows():
    rows = {}
    for r in csv.DictReader(open(M.CACHE/'Jij_table.csv')):
        z1, z2 = int(float(r['Zimp1'])), int(float(r['Zimp2']))
        key = f'{M.NAMEOF[z1]}-{M.NAMEOF[z2]}'
        row = {'i_shell': int(r['i_shell']), 'R': float(r['R']), 'Rx': float(r['Rx']), 'Ry': float(r['Ry']),
               'Rz': float(r['Rz']), 'ilayer2': int(float(r['ilayer2'])), 'Zimp1': z1, 'Zimp2': z2,
               'num_neighbors': int(float(r['num_neighbors'])), 'm1': float(r['m1']), 'm2': float(r['m2']),
               'm01': float(r['m01']), 'm02': float(r['m02']), 'J': float(r['J']), 'J_data_uuid': r['J_data_uuid']}
        rows.setdefault(key, []).append(row)
    for v in rows.values():
        v.sort(key=lambda x: x['i_shell'])
    return rows


def check_archive_convention():
    """The archived mean-field matrix equals (1/3k_B) sum_shells N J (pair-once Hamiltonian, paper Eq. 2-3), with the
    authors' (T_AA + T_BB + T_AB)/3 co-doping rule. Its non-magnetic mask is permuted (notebook bug: the mask is built
    from a dict in AiiDA group order): Ru (0.225 muB) is blanked and Zr (3e-8 muB) is kept."""
    T = np.loadtxt(M.CACHE/'Tcs_Bi2Te3.txt')
    N = [s[2] for s in M.shells()]
    mom = M.single_moments()
    blank = {a for i, a in enumerate(M.NAMES) if not np.any(T[i])}
    assert 'Ru' in blank and 'Zr' not in blank and mom['Ru'] > 0.2 and mom['Zr'] < 1e-6, blank

    def t(a, b):
        return sum((j or 0.0) * n for j, n in zip(M.pair_J(a, b), N)) / 3 / M.KB
    worst = 0.0
    for i, a in enumerate(M.NAMES):
        for j, b in enumerate(M.NAMES):
            if a not in blank and b not in blank:
                ref = t(a, a) if i == j else (t(a, b) + t(a, a) + t(b, b)) / 3
                worst = max(worst, abs(T[i, j] - ref))
    assert worst < 1e-3, worst
    return worst


def check_m2_column(rows):
    """Archive defect: m2 is impurity 2's moment only in shell 0; beyond, it is ~0.01 muB for every pair."""
    for key, rs in rows.items():
        for r in rs:
            if r['m02'] > 0.5:
                if r['i_shell'] == 0:
                    assert abs(r['m2'] / r['m02'] - 1) < 0.25, (key, r)
                else:
                    assert abs(r['m2']) < 0.02, (key, r)


def cif_block():
    lines = (M.CACHE/'Bi2Te3.cif').read_text().split('\n')
    return [l for l in lines if l.strip() and not l.startswith('#')]


def main():
    rows = read_rows()
    worst = check_archive_convention()
    check_m2_column(rows)
    ref = json.loads((HERE/'reference.json').read_text())
    abl = json.loads((HERE/'ablations.json').read_text())
    species = {n: M.ZOF[n] for n in M.NAMES}
    catalog = {
        'material': 'Bi2Te3 with substitutional 3d/4d transition-metal dopants on Bi sites',
        'source': 'R. Mozumder, J. Wasmer, D. Antognini Silva, S. Bluegel, P. Ruessmann, High-throughput magnetic co-doping '
                  'and design of exchange interactions in topological insulators, Materials Cloud Archive 2024 '
                  '(doi:10.24435/materialscloud:b7-6k), arXiv:2407.04413; KKR impurity embedding (JuKKR/aiida-kkr), '
                  'spin-orbit coupling included, exchange from infinitesimal rotations in impurity-dimer calculations',
        'host_structure_cif': cif_block(),
        'host_hexagonal_A': {'a': M.A_HEX, 'c': M.C_HEX,
                             'note': 'Bi1 at z = 0.4 and Bi2 at z = 0.6 of the rhombohedral cell; the two Bi layers of a '
                                     'quintuple layer are related by inversion'},
        'dopants': species,
        'archive': {
            'file': 'Jij_table.csv (2077 rows)',
            'columns': COLUMNS,
            'authors_hamiltonian': 'H = -(1/2) sum_{i,j} J_ij S_i . S_j with S_i = m_i/|m_i| the unit vector along the '
                                   'moment of impurity i (paper Eq. 2)',
            'authors_note': 'single impurities with spin moment below 1e-3 muB are treated as non-magnetic (paper, Sec. III D)',
        },
        'purchases': {
            'pair': {'request': {'kind': 'pair', 'species': ['Cr', 'V'], 'shells': '0-6'},
                     'shell_windows': {'0-6': 'shells 0-6 (R = 4.38-7.83 A)', '7-16': 'shells 7-16 (R = 8.77-12.17 A)'},
                     'cost': 1,
                     'returns': 'every archived row of that dopant pair (either order) within the window, verbatim; '
                                'calculations that are not in the archive are simply absent'},
            'monte_carlo': {
                'request': {'kind': 'monte_carlo', 'spins': {'V': 0.05, 'Cr': 0.05},
                            'exchange_meV': {'V-V': [0.0] * 17, 'V-Cr': [0.0] * 17, 'Cr-Cr': [0.0] * 17},
                            'size': 4, 'samples': 16, 'temperatures_K': [30, 20, 15, 10], 'seed': 0},
                'model': 'classical unit-vector spins on randomly chosen Bi sites (both Bi layers), H = -sum_{i<j} J_ij e_i.e_j '
                         'with exchange_meV[pair][k] applied to every pair of spins at the distance of archive shell k '
                         '(k = 0-16) and nothing beyond; keys are species pairs ordered by atomic number',
                'spins': 'concentration of each spin-carrying species as a fraction of all Bi sites (one or two species, '
                         '0.005-0.11 each, at most 0.11 in total); non-magnetic dopants need not be listed',
                'supercell': 'periodic orthorhombic cell of (7 size) x (4 size) x (size) units of (a, sqrt(3) a, c), '
                             '336 size^3 Bi sites; each species gets round(concentration x sites) spins',
                'samples': 'independent random dopant placements; seed selects a different set of placements',
                'sweeps': f"per temperature {MC['equilibration_sweeps']} equilibration + {MC['measurement_sweeps']} measurement "
                          f"sweeps (heat bath + {MC['overrelaxation']} over-relaxation); temperatures are run from high to low, "
                          'each continuing from the previous configuration',
                'cost_per_16_samples_by_size': MC['cost_per_16_samples'],
                'max_samples_by_size': MC['max_samples'],
                'max_temperatures': MC['max_temperatures'],
                'note': 'runs take up to a few minutes; poll the job'},
        },
        'budget': BUDGET,
        'compositions': COMPOSITIONS,
        'codopant_candidates_with_Cr': CANDIDATES,
    }
    mc_defs = {'m': 'order parameter vector (1/n) sum e_i over the spins of the set (n = number of spins in the set); '
                    'for "A-minus-B", (1/n_total)(sum_A e_i - sum_B e_i)',
               'm2_kmin_per_sample': '<|m(k)|^2> of the total magnetization, m(k) = (1/n) sum_i e_i exp(i k.r_i), averaged over '
                                     'the three smallest nonzero wave vectors k = 2 pi / box_A[d] along x, y and z',
               'm2_per_sample, m4_per_sample, abs_m_per_sample': 'thermal averages <|m|^2>, <|m|^4>, <|m|> per sample '
                                                                   '(rows) and temperature (columns)',
               'energy_meV_per_spin': 'sample-averaged thermal energy per spin'}
    data = {'species': species, 'windows': {'0-6': [0, 6], '7-16': [7, 16]}, 'rows': rows, 'monte_carlo': MC,
            'mc_output_definitions': mc_defs}
    hidden = {'targets': ref['targets'], 'tolerances': ref['tolerances'], 'best_codopant_for_Cr': ref['best_codopant'],
              'reference_details': ref['details']}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [float(x) for x in ref['hidden_numbers']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps(abl, indent=1) + '\n')
    for p in (TASK/'environment/lab/public.json', TASK/'tests/public.json'):
        pub = json.loads(p.read_text())
        pub['solver_budget'] = BUDGET
        p.write_text(json.dumps(pub, indent=1) + '\n')
    print('archive mean-field matrix reproduced to', worst, 'K/c; m2 defect asserted; files written')


if __name__ == '__main__':
    main()
