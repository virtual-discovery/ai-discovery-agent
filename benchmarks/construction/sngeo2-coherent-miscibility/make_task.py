"""Write the catalog, lab data, private targets and ablation values for sngeo2-coherent-miscibility.

    python3 -B make_task.py
Needs pool.json.gz (extract.py). Prints the calibration rows.
"""
import json
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/sngeo2-coherent-miscibility'
P = M.pool()
COST = {'PBE': 1, 'r2SCAN': 2}
BUDGET = 58
X7 = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
PLANES = ('001', '100', '101')
THRESH = {'mixing_enthalpy_window': [100.0, 124.0], 'strain_rel_tol': 0.05, 'incoherent_Tc_window_K': [2000.0, 3000.0],
          'coherent_Tc_rel_tol': {'PBE': 0.12, 'r2SCAN': 0.20}, 'boundary_abs_tol': 0.10}


def r(v, n=6):
    return float(round(v, n))


def scan_ids(fn, plane, member):
    S = M.scans(fn, plane, member)
    keys = sorted(S)
    p1 = sorted({k[0] for k in keys})
    out = {}
    if plane == '001':
        for k in keys:
            out[f'{plane}-{member}-{p1.index(k[0]):02d}'] = k
        return out, p1, None
    cs = sorted({k[1] for k in keys})
    for k in keys:
        out[f'{plane}-{member}-{p1.index(k[0]):02d}-{cs.index(k[1]):02d}'] = k
    return out, p1, cs


def main():
    catalog, data = {}, {'sqs': {}, 'strain': {}, 'hse': {}, 'orderings_r2SCAN': {}}
    # ---- end members (relaxed, from the SQS files' x = 0 / 1 records) ----
    ends = {}
    for fn in ('PBE', 'r2SCAN'):
        ends[fn] = {}
        for rec in P['sqs'][fn]:
            if rec['x_Ge'] in (0.0, 1.0):
                A = np.array(rec['lattice_A'])
                ends[fn]['GeO2' if rec['x_Ge'] == 1.0 else 'SnO2'] = {
                    'a_A': r(np.linalg.norm(A[1]), 4), 'c_A': r(A[0, 2], 4), 'energy_eV_per_6atom_cell': rec['energy_as_released']}
    # ---- orderings ----
    orderings = []
    odir = TASK/'environment/agent/orderings_PBE'
    odir.mkdir(exist_ok=True)
    for f in odir.glob('*.vasp'):
        f.unlink()
    for i, rec in enumerate(P['orderings']['PBE']):
        n = sum(s != 'O' for s in rec['species'])
        oid = f'O{i:02d}'
        order = [e for e in ('Ge', 'Sn', 'O') if e in rec['species']]
        idx = [j for e in order for j, s in enumerate(rec['species']) if s == e]
        lines = [f"{oid} PBE relaxed; total energy of this cell = {rec['energy_eV']:.8f} eV", '1.0']
        lines += ['  ' + '  '.join(f'{v:.8f}' for v in row) for row in rec['lattice_A']]
        lines += ['  '.join(order), '  '.join(str(rec['species'].count(e)) for e in order), 'Direct']
        lines += ['  ' + '  '.join(f'{v:.7f}' for v in rec['frac'][j]) + f"  {rec['species'][j]}" for j in idx]
        (odir/f'{oid}.vasp').write_text('\n'.join(lines) + '\n')
        orderings.append({'id': oid, 'n_cations': n, 'n_Ge': rec['species'].count('Ge'), 'poscar': f'/app/orderings_PBE/{oid}.vasp'})
    for i, rec in enumerate(P['orderings']['r2SCAN']):
        data['orderings_r2SCAN'][f'O{i:02d}'] = {'lattice_A': rec['lattice_A'], 'species': rec['species'], 'frac': rec['frac'],
                                                'energy_eV': rec['energy_eV']}
    # ---- SQS ----
    sqs_av = []
    for fn in ('PBE', 'r2SCAN'):
        for rec in P['sqs'][fn]:
            key = f"{fn}|{rec['x_Ge']:.3f}|{rec['index']}"
            data['sqs'][key] = {'functional': fn, 'x_Ge': rec['x_Ge'], 'index': rec['index'], 'lattice_A': rec['lattice_A'],
                                'species': rec['species'], 'frac': rec['frac'], 'energy_as_released': rec['energy_as_released']}
            sqs_av.append({'functional': fn, 'x_Ge': rec['x_Ge'], 'index': rec['index']})
    # ---- strain scans ----
    strain_cat = {}
    for fn in ('PBE', 'r2SCAN'):
        strain_cat[fn] = {}
        for plane in PLANES:
            for member in ('GeO2', 'SnO2'):
                ids, p1, cs = scan_ids(fn, plane, member)
                raw = {}
                for rec in P['strain'][f'{fn}|{plane}|{member}']:
                    k, oop = M.scan_key(plane, rec['lattice_A'])
                    raw.setdefault(k, []).append((oop, rec))
                for sid, k in ids.items():
                    pts = [rec for oop, rec in sorted(raw[k], key=lambda t: t[0])]
                    data['strain'][f'{fn}|{sid}'] = {'functional': fn, 'id': sid, 'points': [
                        {'lattice_A': p['lattice_A'], 'species': p['species'], 'frac': p['frac'], 'energy_eV': p['energy_eV']} for p in pts]}
            if plane == '001':
                desc = {'fixed': {'|a2| = |a3| (A), a2 perpendicular to a3': [r(v, 4) for v in p1]},
                        'stepped': '|a1| (the rutile c axis)', 'id_format': '001-<member>-<i>'}
            elif plane == '100':
                desc = {'fixed': {'|a2| (A)': [r(v, 4) for v in p1], '|a1| = c (A)': [r(v, 4) for v in cs]},
                        'stepped': '|a3|', 'id_format': '100-<member>-<i>-<j>  (i indexes |a2|, j indexes |a1|)'}
            else:
                desc = {'fixed': {'|a2 - a3| (A)': [r(v, 4) for v in p1], '|a1| = c (A)': [r(v, 4) for v in cs]},
                        'stepped': '|a2 + a3| (|a2| = |a3|)', 'id_format': '101-<member>-<i>-<j>  (i indexes |a2 - a3|, j indexes |a1|)'}
            n_pts = {m: len(next(iter(M.scans(fn, plane, m).values()))) for m in ('GeO2', 'SnO2')}
            strain_cat[fn][plane] = {'archive_file': f'data_{plane}_{fn}.json',
                                     'archive_description': f'data used to compute the coherency strain energy assuming coexistence along the ({plane}) plane',
                                     **desc, 'points_per_scan': n_pts}
    # ---- HSE-setting PBE runs ----
    hse_av = []
    for name, rec in P['hse_setting_pbe'].items():
        data['hse'][name] = {'name': name, **rec}
        hse_av.append({'name': name, 'n_Ge': rec['species'].count('Ge'), 'n_Sn': rec['species'].count('Sn')})
    catalog = {
        'material': 'rutile Sn1-xGexO2 (x = Ge fraction); author-released VASP results of Muller, Kurbay, Zhang, Kioupakis and Natarajan (2026)',
        'note': 'Every purchase retrieves an archived calculation; no new DFT is run.',
        'purchase_costs': {'sqs': COST, 'strain_scan': COST, 'ordering_r2SCAN': 2, 'hse_setting_run': 1},
        'budget': BUDGET,
        'end_members': ends,
        'end_members_note': 'Relaxed rutile end members (6-atom cell, 2 formula units) in the settings of the SQS, ordering and strain data of each functional.',
        'orderings_PBE': orderings,
        'orderings_note': ('The authors\' cluster-expansion training set: the 80 symmetrically distinct Ge/Sn orderings on the rutile cation '
                           'sublattice in cells of up to 6 cations, relaxed with PBE, staged as POSCAR files whose first line gives the total '
                           'energy of that cell (eV). Lattice-row order differs between cells; the rutile c axis is along Cartesian z in every cell.'),
        'orderings_r2SCAN_available': [o['id'] for o in orderings],
        'sqs_available': sqs_av,
        'sqs_note': ('72-atom special quasirandom structures (24 cations, 3 per composition) and the end members, relaxed. A purchase returns the '
                     'relaxed structure and the field energy_as_released copied unchanged from the archive file data_SQS_<functional>.json, '
                     'which the archive README describes as "the associated VASP energy computed with" that functional.'),
        'strain_scans': {'note': ('The archive\'s coherency-strain datasets: rutile GeO2 and SnO2 (6-atom cells, 2 formula units) relaxed with some '
                                  'lattice parameters held fixed. Each scan holds the structures with the listed fixed values and one further lattice '
                                  'parameter stepped; a purchase returns every structure of the scan with its total energy (eV per 6-atom cell). '
                                  'a1, a2, a3 are the POSCAR lattice rows; a1 lies along the rutile c axis. Plane labels are the archive\'s.'),
                         'datasets': strain_cat},
        'hse_setting_runs_available': hse_av,
        'hse_setting_note': ('Separate PBE relaxations of the same SQS (folder Sn_<Sn percent>_<n>) and of the end members, run with the '
                             'hybrid-functional settings (950 eV cutoff, GW PAW datasets). A purchase returns the final energy (sigma -> 0), '
                             'the relaxed structure and the run settings.'),
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')) + '\n')
    # ---- reference and wrong paths ----
    idx0 = {(x, 0) for x in X7}
    ref = M.analysis()
    cases = {
        'sqs_values_read_per_cation': dict(sqs_factor=2.0),
        'strain_per_cell_mixed_with_per_cation': dict(strain_factor=2.0),
        'r2scan_unrelaxed_sqs_kept': dict(sqs_mode={'PBE': 'clean', 'r2SCAN': 'all'}),
        'r2scan_index0_sqs_unscreened': dict(sqs_mode={'PBE': 'clean', 'r2SCAN': 'all'}, sqs_indices={'PBE': None, 'r2SCAN': idx0}),
        'inplane_c_not_constrained': dict(strain_mode='aonly'),
        'no_out_of_plane_relaxation': dict(strain_mode='norelax'),
        'habit_plane_100': dict(coherent_plane='100'),
        'incoherent_model_reported_as_coherent': dict(no_coherency=True),
        'binodal_spinodal_swapped': dict(swap_700K=True),
        'functional_mixing': dict(strain_functional={'PBE': 'r2SCAN', 'r2SCAN': 'PBE'}),
    }
    robust_cases = {
        'lowest_sqs_per_composition': dict(sqs_mode={'PBE': 'cleanmin', 'r2SCAN': 'cleanmin'}),
        'mean_sqs_per_composition': dict(sqs_mode={'PBE': 'cleanmean', 'r2SCAN': 'cleanmean'}, dh05_from='mean'),
        'raw_scan_minima_three_compositions': dict(strain_method='raw', strain_xsel=[0.25, 0.5, 0.75]),
        'one_sqs_per_composition': dict(sqs_indices={'PBE': idx0 | {(0.5, 1), (0.5, 2)}, 'r2SCAN': (idx0 - {(0.75, 0)}) | {(0.75, 2)}},
                                        strain_xsel=[0.25, 0.5, 0.75]),
    }
    abl = {k: M.analysis(**v) for k, v in cases.items()}
    rob = {k: M.analysis(**v) for k, v in robust_cases.items()}
    # Defensible: beyond-mean-field incoherent Tc from a cluster expansion with Monte Carlo (construction estimate,
    # README: Binder-cumulant crossings of SGC Monte Carlo for Lasso CEs, 2350-2550 K)
    rob['cluster_expansion_monte_carlo_Tc'] = {**ref, 'incoherent_Tc_K': 2450.0}
    rob['cluster_expansion_12_clusters_monte_carlo_Tc'] = {**ref, 'incoherent_Tc_K': 2250.0}
    # Wrong path (ce_mc/): over-truncated CE (Lasso lambda 3, 8 clusters, 5-fold CV 20 meV/cation) + the same Monte Carlo
    abl['overtruncated_cluster_expansion_Tc'] = {**ref, 'incoherent_Tc_K': 1650.0}
    xs, ec = {}, {}
    rob['strain_from_subregular_fit'] = {**ref, 'coherency_strain_energy_x050_meV_per_cation':
                                         {pl: float(M.strain_fit('PBE', pl)[0] / 8 * 1000) for pl in PLANES}}
    hidden = {'reference': ref, 'thresholds': THRESH}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    hv = [ref['mixing_enthalpy_x050_meV_per_cation'], *ref['coherency_strain_energy_x050_meV_per_cation'].values(),
          ref['incoherent_Tc_K'], *ref['coherent_spinodal_Tc_K'].values(),
          *ref['coherent_boundaries_700K_PBE']['binodal'], *ref['coherent_boundaries_700K_PBE']['spinodal']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(hv) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob, 'reference': ref}, indent=1) + '\n')
    def short(d):
        return json.dumps({k: (round(v, 3) if isinstance(v, float) else
                               ({kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()} if isinstance(v, dict) and k != 'coherent_boundaries_700K_PBE' else
                                ({kk: (None if vv is None else [round(q, 3) for q in vv]) for kk, vv in v.items()} if isinstance(v, dict) else v)))
                           for k, v in d.items()})
    print('reference', short(ref))
    for k, v in {**abl, **rob}.items():
        print(k, short(v))


if __name__ == '__main__':
    main()
