"""Generate the staged, lab, private and ablation data for hydrogen-plasma-facing-materials.

    python3 -B make_task.py [CACHE]     (run in ai-theorist-benchmark-runtime:20261007 from this folder)

Reads CACHE/pool.json (extract.py) and CACHE/raw/../README.md, writes into benchmarks/tasks/<id>/:
environment/agent/catalog.json, environment/lab/data.json, environment/lab/private/hidden.json,
environment/lab/private/hidden_values.json and solution/ablation_values.json, and prints the calibration.
"""
import json
import sys
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/hydrogen-plasma-facing-materials'
CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else M.CACHE
P = M.load_pool(CACHE)
BUDGET = 61

H_SET = ['C-h-194', 'HfB2-h-191', 'Mo-c-229', 'Re-h-194', 'TaC-c-225', 'TiN-c-225', 'WC-h-187']
SBE_SET = {'W-c-229': ['W'], 'WC-h-187': ['W', 'C'], 'TiB2-h-191': ['B', 'Ti']}
COH_SET = ['Re-h-194', 'W-c-229', 'WC-h-187']
# TiB2: the archived Ti atom is a byte-identical copy of the TiN Ti atom (60/480 Ry) while every TiB2 slab,
# supercell and unit cell uses 35/280 Ry. cutoff_check/ shows the Ti atom energy changes by only a few meV
# between 35 and 60 Ry, so the archived combination is numerically sound; a null (not determined) answer is
# also accepted for this one entry because the mismatch cannot be quantified from the archive alone.
UNDETERMINED = set()
NULL_ALSO_ACCEPTED = ['TiB2-h-191/Ti']
THRESHOLDS = {'h_formation_tol_eV': 0.03, 'site_energy_tol_eV': 0.03, 'distance_tol_A': 0.04,
              'sbe_tol_eV': 0.06, 'cohesive_tol_eV': 0.03}
FORMULA = {'BC2N-h-160': 'BC2N', 'BN-c-216': 'BN', 'C-c-227': 'C (diamond)', 'C-h-194': 'C (graphite)',
           'C3N4-h-176': 'C3N4', 'Cr-c-229': 'Cr', 'HfB2-h-191': 'HfB2', 'HfC-c-225': 'HfC', 'Mo-c-229': 'Mo',
           'Re-h-194': 'Re', 'Ta-c-229': 'Ta', 'TaC-c-225': 'TaC', 'TaN-h-187': 'TaN', 'TaW-c-221': 'TaW (ordered)',
           'TiB2-h-191': 'TiB2', 'TiN-c-225': 'TiN', 'VB2-h-191': 'VB2', 'W-c-229': 'W', 'W2C-o-60': 'W2C',
           'WC-h-187': 'WC', 'ZrB2-h-191': 'ZrB2'}


def r(x, n=4):
    return None if x is None else round(float(x), n)


def lowest_termination(mat):
    terms = M.terminations(P, mat)
    return min(terms, key=lambda t: M.surface_energy(P, mat, t))


def min_sbe(mat, term, el, atom_id=None):
    labs = [l for l in M.vacancies(P, mat, term) if l.split('_')[0] == el]
    return min(M.sbe(P, mat, term, l, atom_id) for l in labs)


def mean_sbe(mat, terms, el, atom_id=None):
    vals = [M.sbe(P, mat, t, l, atom_id) for t in terms for l in M.vacancies(P, mat, t) if l.split('_')[0] == el]
    return float(np.mean(vals))


def preferred(mat, site=None, initial=False):
    cl = M.distinct_sites(P, mat)
    site = site or cl[0]['rep_site']
    d, el = M.h_environment(P[f'H-IFE/{mat}/site-{site}'], initial=initial)[0]
    return {'nearest_host_element': el, 'nearest_host_distance_A': r(d, 4)}


def reference():
    ref = {'h_formation_energy_eV': {}, 'distinct_site_energies_eV': {}, 'preferred_site': {}, 'preferred_site_runs': {},
           'surface_binding_energy_eV': {}, 'cohesive_energy_eV_per_atom': {}, 'lowest_termination': {}}
    for mat in H_SET:
        cl = M.distinct_sites(P, mat)
        ref['h_formation_energy_eV'][mat] = r(cl[0]['rep_energy'])
        ref['distinct_site_energies_eV'][mat] = [r(c['rep_energy']) for c in cl]
        ref['preferred_site'][mat] = preferred(mat)
        ref['preferred_site_runs'][mat] = cl[0]['sites']
    for mat, els in SBE_SET.items():
        t = lowest_termination(mat)
        ref['lowest_termination'][mat] = t
        ref['surface_binding_energy_eV'][mat] = {el: None if (mat, el) in UNDETERMINED else r(min_sbe(mat, t, el))
                                                 for el in els}
    for mat in COH_SET:
        ref['cohesive_energy_eV_per_atom'][mat] = None if (mat, 'cohesive') in UNDETERMINED else r(M.cohesive(P, mat))
    return ref


def ablations(ref):
    A, info = {}, {}
    # dp1: isolated (spin-polarized) H atom as the reservoir instead of 1/2 H2
    A['atom_reference'] = {
        'h_formation_energy_eV': {m: r(M.distinct_sites(P, m, ref='H')[0]['rep_energy']) for m in H_SET},
        'distinct_site_energies_eV': {m: [r(c['rep_energy']) for c in M.distinct_sites(P, m, ref='H')] for m in H_SET}}
    # dp2: the 'textbook' site (tetrahedral in the hcp metal, metal-coordinated void in the carbide)
    A['textbook_site'] = {'h_formation_energy_eV': {'Re-h-194': r(M.hife(P, 'Re-h-194', '17')),
                                                    'TaC-c-225': r(M.hife(P, 'TaC-c-225', '1'))},
                          'preferred_site': {'Re-h-194': preferred('Re-h-194', '17'), 'TaC-c-225': preferred('TaC-c-225', '1')}}
    # dp3: every archived relaxation counted as a distinct site; or energy-only binning at 0.1 eV
    A['archived_runs_as_sites'] = {'distinct_site_energies_eV': {
        m: sorted(r(M.hife(P, m, s)) for s in M.site_ids(P, m)) for m in H_SET}}
    A['energy_binning_0p1eV'] = {'distinct_site_energies_eV': {
        m: [r(c['rep_energy']) for c in M.distinct_sites(P, m, e_tol=0.1, r_tol=None)] for m in H_SET}}
    # dp4: one W atom reused across materials (the 30 Ry W-set atom for WC, as in the paper's table; or the reverse)
    w30, w45 = 'SBE/W-c-229/atom/W', 'SBE/WC-h-187/atom/W'
    tw = ref['lowest_termination']['WC-h-187']
    A['w_atom_from_w_set'] = {'surface_binding_energy_eV': {'WC-h-187': {'W': r(min_sbe('WC-h-187', tw, 'W', w30))}},
                              'cohesive_energy_eV_per_atom': {'WC-h-187': r(M.cohesive(P, 'WC-h-187', atom_ids={'W': w30}))}}
    A['w_atom_from_wc_set'] = {'surface_binding_energy_eV': {'W-c-229': {'W': r(min_sbe('W-c-229', '110', 'W', w45))}},
                               'cohesive_energy_eV_per_atom': {'W-c-229': r(M.cohesive(P, 'W-c-229', atom_ids={'W': w45}))}}
    info['TiB2_cohesive_mixed_vs_60Ry_route'] = [r(M.cohesive(P, 'TiB2-h-191')), r(
        (M.energy(P['SBE/TiB2-h-191/atom/Ti']) + 2 * M.energy(P['SBE/BN-c-216/atom/B'])
         - M.energy(P['H-IFE/TiB2-h-191/supercell']) / 27) / 3 * M.RY)]
    # dp5: wrong termination; the paper's average over all sites and terminations; site mean
    other = {m: [t for t in M.terminations(P, m) if t != ref['lowest_termination'][m]] for m in SBE_SET}
    A['other_termination'] = {'surface_binding_energy_eV': {
        m: {el: (None if (m, el) in UNDETERMINED else r(min_sbe(m, other[m][0], el))) for el in els}
        for m, els in SBE_SET.items() if other[m]}}
    A['paper_protocol_mean'] = {'surface_binding_energy_eV': {
        'WC-h-187': {'W': r(mean_sbe('WC-h-187', M.terminations(P, 'WC-h-187'), 'W', w30)),
                     'C': r(mean_sbe('WC-h-187', M.terminations(P, 'WC-h-187'), 'C'))},
        'TiB2-h-191': {'B': r(mean_sbe('TiB2-h-191', M.terminations(P, 'TiB2-h-191'), 'B')),
                       'Ti': r(mean_sbe('TiB2-h-191', M.terminations(P, 'TiB2-h-191'), 'Ti'))}}}
    A['site_mean_lowest_termination'] = {'surface_binding_energy_eV': {
        m: {el: (None if (m, el) in UNDETERMINED else r(mean_sbe(m, [ref['lowest_termination'][m]], el))) for el in els}
        for m, els in SBE_SET.items()}}
    # dp6: cohesive energy with the H-IFE supercell (60/480 Ry, already bought for the H work) as the bulk
    A['cohesive_bulk_from_hife_supercell'] = {'cohesive_energy_eV_per_atom': {
        m: r(M.cohesive(P, m, bulk_id=f'H-IFE/{m}/supercell')) for m in ('Re-h-194', 'WC-h-187')}}
    info['Mo_cohesive_consistent_vs_hife_supercell'] = [r(M.cohesive(P, 'Mo-c-229')),
                                                        r(M.cohesive(P, 'Mo-c-229', bulk_id='H-IFE/Mo-c-229/supercell'))]
    info['W_hife_supercell_mixing'] = r(M.cohesive(P, 'W-c-229', bulk_id='H-IFE/W-c-229/supercell'))
    # dp8 (competing explanation): bulk cohesive energy per atom (sublimation proxy) used as the SBE of every element
    A['cohesive_energy_as_sbe'] = {'surface_binding_energy_eV': {
        m: {el: r(M.cohesive(P, m)) for el in els} for m, els in SBE_SET.items()}}
    # dp7: the site read from the starting (input) coordinates instead of the relaxed ones
    A['starting_geometry'] = {'preferred_site': {m: preferred(m, initial=True) for m in H_SET}}
    return A, info


def robustness(ref):
    R = {}
    R['duplicate_cluster_mean'] = {
        'h_formation_energy_eV': {m: r(np.mean(M.distinct_sites(P, m)[0]['energies'])) for m in H_SET},
        'distinct_site_energies_eV': {m: [r(np.mean(c['energies'])) for c in M.distinct_sites(P, m)] for m in H_SET}}
    t = ref['lowest_termination']['WC-h-187']
    R['carbon_atom_at_60Ry'] = {'surface_binding_energy_eV': {
        'WC-h-187': {'C': r(min_sbe('WC-h-187', t, 'C', 'SBE/C3N4-h-176/atom/C'))}}}
    R['distance_from_another_duplicate'] = {'preferred_site': {
        m: preferred(m, M.distinct_sites(P, m)[0]['sites'][-1]) for m in H_SET}}
    R['rydberg_13p6057'] = {'h_formation_energy_eV': {m: r(v * 13.6057 / M.RY) for m, v in ref['h_formation_energy_eV'].items()},
                            'cohesive_energy_eV_per_atom': {m: (None if v is None else r(v * 13.6057 / M.RY))
                                                            for m, v in ref['cohesive_energy_eV_per_atom'].items()}}
    return R


def catalog():
    readme = (CACHE/'README.md').read_text()
    mats = sorted({k.split('/')[1] for k in P if k.startswith('H-IFE/') and k.count('/') == 2})
    entries = {}
    for m in mats:
        sites = [{'id': f'H-IFE/{m}/site-{s}', 'generation_number': int(s)} for s in M.site_ids(P, m)]
        slabs = []
        for k in sorted(P):
            if k.startswith(f'SBE/{m}/slab/'):
                lab = k.split('/')[3]
                slabs.append({'id': k, 'miller_index': lab.split('_')[0], 'termination': lab.split('_')[1] if '_' in lab else None})
        vac = []
        for t in M.terminations(P, m):
            for lab in M.vacancies(P, m, t):
                el, idx, depth = lab.split('_')
                vac.append({'id': f'SBE/{m}/vacancy/{t}/{lab}', 'slab_supercell': f'SBE/{m}/slab_supercell/{t}',
                            'removed_element': el, 'site_number': int(idx), 'depth_below_top_atom_A': float(depth)})
        facet = M.terminations(P, m)[0].split('_')[0]
        entries[m] = {
            'formula': FORMULA[m], 'space_group_number': int(m.split('-')[-1]),
            'hydrogen_interstitial': {'pristine_supercell': f'H-IFE/{m}/supercell', 'interstitial_runs': sites},
            'surface_binding': {
                'bulk_unit_cell': f'SBE/{m}/unitcell',
                'isolated_atoms': {k.split('/')[-1]: k for k in sorted(P) if k.startswith(f'SBE/{m}/atom/')},
                'slabs': slabs, 'selected_facet': facet,
                'slab_supercells': [{'id': f'SBE/{m}/slab_supercell/{t}', 'termination': t} for t in M.terminations(P, m)],
                'surface_vacancies': vac}}
    return {
        'source': 'A. Fedrigucci, N. Marzari, P. Ricci, Materials Cloud Archive record gpxwh-dx645 '
                  '(DOI 10.24435/materialscloud:9k-09, CC BY 4.0); Quantum ESPRESSO, PBEsol, SSSP efficiency pseudopotentials',
        'units': 'Energies are reported as printed by Quantum ESPRESSO, in Ry (1 Ry = 13.605693122994 eV). '
                 'Cells in Angstrom, positions in crystal (fractional) coordinates, forces in Ry/bohr, pressure in kbar.',
        'what_a_purchase_returns': 'The released run: the archive folder, the run header (code version, start time, MPI '
                                   'processes, wall time), the input settings (calculation type, cutoffs, k-points, '
                                   'smearing, spin, convergence thresholds, pseudopotentials), the cell and starting '
                                   'positions, and the output (final and last SCF total energies, number of SCF cycles, '
                                   'relaxation message, final total force, magnetization, Fermi energy, pressure, final '
                                   'positions and, for vc-relax, the final cell).',
        'archive_readme_verbatim': readme,
        'purchase_cost': 1, 'budget': BUDGET,
        'hydrogen_references': {'H2_molecule': 'H-IFE/H2', 'H_atom': 'H-IFE/H'},
        'materials': entries,
        'graded_materials': {'hydrogen': H_SET, 'surface_binding': {m: v for m, v in SBE_SET.items()},
                             'cohesive_energy': COH_SET}}


def lab_data():
    return {k: {'archive_path': v['archive_path'], 'run': v['run'], 'input': v['input'], 'cell_A': v['cell_A'],
                'positions_crystal_initial': v['positions_crystal_initial'], 'output': v['output']} for k, v in P.items()}


def hidden_numbers(ref):
    vals = list(ref['h_formation_energy_eV'].values())
    for v in ref['distinct_site_energies_eV'].values():
        vals += v
    vals += [p['nearest_host_distance_A'] for p in ref['preferred_site'].values()]
    for d in ref['surface_binding_energy_eV'].values():
        vals += [x for x in d.values() if x is not None]
    vals += [x for x in ref['cohesive_energy_eV_per_atom'].values() if x is not None]
    return sorted(set(vals))


def main():
    ref = reference()
    A, info = ablations(ref)
    R = robustness(ref)
    hidden = {'reference': ref, 'thresholds': THRESHOLDS, 'null_also_accepted': NULL_ALSO_ACCEPTED}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog(), indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab_data(), separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(hidden_numbers(ref)) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': R, 'reference': ref},
                                                                 indent=1) + '\n')
    print(json.dumps({'reference': ref, 'ablations': A, 'robustness': R, 'info': info}, indent=1))


if __name__ == '__main__':
    main()
