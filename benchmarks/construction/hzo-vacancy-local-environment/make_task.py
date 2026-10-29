"""Generate the staged, lab, private and ablation data for hzo-vacancy-local-environment from pool.json.

    cd benchmarks/construction/hzo-vacancy-local-environment && python3 -B make_task.py

Writes into benchmarks/tasks/hzo-vacancy-local-environment/:
  environment/agent/catalog.json, environment/agent/structures/*.cif, environment/agent/x050_site_table.csv,
  environment/lab/data.json, environment/lab/private/hidden.json, environment/lab/private/hidden_values.json,
  solution/ablation_values.json
and prints the calibration table (design.json -> calibration).
"""
import csv
import io
import json
import math
from pathlib import Path
import numpy as np
import model as M

TASK = M.HERE.parents[1]/'tasks/hzo-vacancy-local-environment'
P = M.POOL
LABEL = {'0.25': 'x0.25', '0.50': 'x0.50', '0.75': 'x0.75'}
CIFNAME = {'0.25': 'HZO_x0.25_Hf24Zr8O64.cif', '0.50': 'HZO_x0.50_Hf16Zr16O64.cif', '0.75': 'HZO_x0.75_Hf8Zr24O64.cif'}
# Withheld sites (graded). Drawn once (seed 7) by class: x = 0.25 includes its only 4-Zr site and Zr-rich
# sites; x = 0.75 includes both Hf-rich (fourfold, 1 Zr) sites. Frozen here for reproducibility.
HELD = {'0.25': [34, 37, 44, 47, 70, 76, 85, 87, 93, 94], '0.75': [36, 43, 46, 52, 60, 64, 68, 76, 92, 94]}
# The oracle's purchase rule (solution/solve.py implements the same rule): for each composition and
# (coordination, n_Zr) class, the first k purchasable sites in ascending index order.
PLAN = {'0.25': [((4, 0), 1), ((4, 1), 2), ((4, 2), 2), ((4, 3), 1), ((3, 0), 1), ((3, 1), 1), ((3, 2), 1)],
        '0.75': [((4, 4), 3), ((4, 3), 2), ((4, 2), 1), ((3, 1), 1), ((3, 2), 1), ((3, 3), 1)]}
# The oracle's CASTEP-LDA purchases: a fixed spread of oxygen indices over both halves of the site list (the
# environment labels of a CASTEP cell are only known after purchase).
LDA_SITES = [33, 38, 51, 68, 85, 41, 46, 58, 75, 92]
BUDGET = 36
COSTS = {'vacancy': 1, 'pristine': 1, 'castep_vacancy': 1, 'castep_333_vacancy': 2}
THRESHOLDS = {'pred_rmse_max_eV': 0.035, 'pred_bias_max_eV': 0.02, 'distribution_sum_tol': 0.02,
              'mean_nZr_tol': {'x0.25': 0.25, 'x0.75': 0.15}, 'share_4Zr_tol': {'x0.25': 0.03, 'x0.75': 0.08},
              'log10_site_fraction_tol': {'x0.25': 0.20, 'x0.75': 0.15},
              'lda_offset_tol_eV': 0.075, 'lda_slope_tol_eV': 0.05}
TABLE_COLUMNS = {
    'O_index': '0-based atom index of the oxygen in HZO_x0.50_Hf16Zr16O64.cif',
    'E_vac': 'formation energy (eV) of a neutral vacancy at this site, as released in the authors\' table',
    'n_Hf': 'number of Hf neighbours (authors\' count)', 'n_Zr': 'number of Zr neighbours (authors\' count)',
    'mean_MO_distance': 'mean metal-oxygen distance (A) of the site', 'min_MO_distance': 'minimum metal-oxygen distance (A)',
    'mean_angle': 'mean metal-oxygen-metal angle (deg)', 'angle_variance': 'variance of the M-O-M angles',
    'cluster': 'authors\' k-means cluster number', 'angle_group': 'authors\' grouping of the mean angle (109 or 120)'}


def common_slope(items):
    """items: (is_O3, n_Zr, E). Least squares E = E0 + D*is_O3 + b*n_Zr; returns (D, b)."""
    A = np.array([[1.0, o3, z] for o3, z, e in items])
    c = np.linalg.lstsq(A, np.array([e for o3, z, e in items]), rcond=None)[0]
    return float(c[1]), float(c[2])


def lda_items(rows, labels='authors'):
    env50 = M.ENV['0.50']
    if labels == 'authors':
        return [(float(r['vacancy_type'] == 'O3'), r['n_Zr'], r['final_enthalpy_eV']) for r in rows]
    return [(float(env50[r['site']][0] == 3), env50[r['site']][1], r['final_enthalpy_eV']) for r in rows]


def lda_field(D, b):
    return {'offset_O3_eV': D, 'slope_per_Zr_eV': b}


def purchases():
    out = {}
    for x, plan in PLAN.items():
        chosen = []
        for cls, k in plan:
            cand = [s for s in sorted(M.EF[x]) if s not in HELD[x] and M.ENV[x][s] == cls]
            assert len(cand) >= k, (x, cls)
            chosen += cand[:k]
        out[x] = chosen
    return out


def rows_for(x, sites):
    return [(float(x), M.ENV[x][s], M.EF[x][s]) for s in sites]


def preds_from(fn):
    return {f'{LABEL[x]}-O{s}': float(fn(x, M.ENV[x][s])) for x in HELD for s in HELD[x]}


def pops_from_class_energy(fn, freq=M.binom, extra=0.0):
    out = {}
    for x in ('0.25', '0.75'):
        f = M.factors_from_energies(lambda cls: fn(x, cls))
        pop = M.population(float(x), f, freq)
        out[LABEL[x]] = {'distribution': pop['distribution'], 'site_fraction': pop['site_fraction'] * math.exp(extra)}
    return out


def payload_parts(pred_fn, pop=None):
    pop = pop or pops_from_class_energy(pred_fn)
    return {'predicted_formation_energies_eV': preds_from(pred_fn),
            'vacancy_Zr_neighbour_distribution': {k: v['distribution'] for k, v in pop.items()},
            'vacancy_site_fraction': {k: v['site_fraction'] for k, v in pop.items()}}


def evaluate(parts, hidden):
    """Same arithmetic as task_lab.grade, for the calibration printout."""
    t, ref = hidden['thresholds'], hidden['reference']
    out = {}
    for x in ('x0.25', 'x0.75'):
        err = [parts['predicted_formation_energies_eV'][k] - v for k, v in hidden['site_targets'].items() if k.startswith(x)]
        d = parts['vacancy_Zr_neighbour_distribution'][x]
        out[x] = {'rmse': float(np.sqrt(np.mean(np.square(err)))), 'bias': float(np.mean(err)),
                  'mean_nZr': sum(k * p for k, p in enumerate(d)) / sum(d), 'share_4Zr': d[4] / sum(d),
                  'dlog10_sf': math.log10(parts['vacancy_site_fraction'][x] / ref[x]['site_fraction'])}
    lda = parts['lda_common_slope_model']
    out['lda'] = {'dD': lda['offset_O3_eV'] - ref['lda']['offset_O3_eV'], 'db': lda['slope_per_Zr_eV'] - ref['lda']['slope_per_Zr_eV']}
    return out


def main():
    # ---------------- reference (hidden) values
    site_targets = {f'{LABEL[x]}-O{s}': M.EF[x][s] for x in HELD for s in HELD[x]}
    reference = {}
    for x in ('0.25', '0.75'):
        pop, absent = M.truth_population(x)
        reference[LABEL[x]] = {'distribution': pop['distribution'], 'mean_nZr': pop['mean_nZr'],
                               'share_4Zr': pop['distribution'][4], 'site_fraction': pop['site_fraction'],
                               'log10_site_fraction': math.log10(pop['site_fraction']), 'absent_class_weight': absent,
                               'o3_share': pop['o3_share']}
    lda_all = {r['site']: r for r in P['castep']['LDA']}
    D_ref, b_ref = common_slope(lda_items(lda_all.values()))
    reference['lda'] = {'offset_O3_eV': D_ref, 'slope_per_Zr_eV': b_ref}
    hidden = {'site_targets': site_targets, 'reference': reference, 'thresholds': THRESHOLDS,
              'units': 'eV; site fraction per O site; distribution over n_Zr = 0..4'}
    # ---------------- oracle evidence (same rule as solve.py) and models
    bought = purchases()
    table_rows = rows_for('0.50', sorted(M.EF['0.50']))          # free table, converted to O-poor
    train = table_rows + rows_for('0.25', bought['0.25']) + rows_for('0.75', bought['0.75'])
    models = {kind: M.fit(train, kind) for kind in ('local_x', 'local_offsets', 'local', 'pooled')}
    fn = {kind: (lambda kind: lambda x, cls: M.predict(models[kind], float(x), cls, kind))(kind) for kind in models}
    resid = np.array([e - M.predict(models['local_x'], x, cls, 'local_x') for x, cls, e in train])
    sigma = float(np.sqrt(np.mean(resid ** 2)))
    jensen = sigma ** 2 / (2 * M.KT ** 2)          # <exp(-E/kT)> = exp(-<E>/kT) exp(sigma^2/2kT^2) for Gaussian scatter
    oracle = payload_parts(fn['local_x'], pops_from_class_energy(fn['local_x'], extra=jensen))
    lda_bought = [lda_all[s] for s in LDA_SITES]
    oracle['lda_common_slope_model'] = lda_field(*common_slope(lda_items(lda_bought)))
    # ---------------- ablations (realistic wrong paths on the oracle's evidence)
    abl = {}
    abl['dp1_pooled_nZr_regression'] = payload_parts(fn['pooled'])
    abl['dp2_composition_blind_pooled_model'] = payload_parts(fn['local'])
    c50 = M.fit(table_rows, 'local')
    f50 = lambda x, cls: M.predict(c50, float(x), cls, 'local')
    abl['dp3_x050_model_transferred'] = payload_parts(f50)
    # dp1 (second path): geometric descriptors of the staged cells instead of cation identity. In the ideal cells
    # the M-O distances and M-O-M angles are identical within a sublattice, so this is a sublattice-mean model
    # (with per-composition offsets).
    geo_rows = [(x, (cls[0], 0), e) for x, cls, e in train]
    cgeo = M.fit(geo_rows, 'local_offsets')
    fgeo = lambda x, cls: M.predict(cgeo, float(x), (cls[0], 0), 'local_offsets')
    abl['dp1_geometry_descriptor_model'] = payload_parts(fgeo)
    # dp4: statistics of the 96-atom cell's own decoration (model energies for its 64 sites)
    pop_cell, pop_hyper = {}, {}
    for x in ('0.25', '0.75'):
        es = np.array([fn['local_x'](x, M.ENV[x][s]) for s in sorted(M.EF[x])])
        nz = np.array([M.ENV[x][s][1] for s in sorted(M.EF[x])])
        w = np.exp(-es / M.KT)
        pop_cell[LABEL[x]] = {'distribution': [float(w[nz == k].sum() / w.sum()) for k in range(5)],
                              'site_fraction': float(w.mean() * math.exp(jensen))}
    pop_hyper = pops_from_class_energy(fn['local_x'], freq=M.hyper, extra=jensen)
    abl['dp4_supercell_decoration_statistics'] = {**oracle, 'vacancy_Zr_neighbour_distribution': {k: v['distribution'] for k, v in pop_cell.items()},
                                                  'vacancy_site_fraction': {k: v['site_fraction'] for k, v in pop_cell.items()}}
    abl['dp4_fixed_composition_cell_statistics'] = {**oracle, 'vacancy_Zr_neighbour_distribution': {k: v['distribution'] for k, v in pop_hyper.items()},
                                                    'vacancy_site_fraction': {k: v['site_fraction'] for k, v in pop_hyper.items()}}
    # dp5: Boltzmann over environment classes without their multiplicity / lowest-energy class only
    nodeg = {}
    for x in ('0.25', '0.75'):
        f = M.factors_from_energies(lambda cls: fn['local_x'](x, cls))
        z = sum(f.values())
        nodeg[LABEL[x]] = [sum(v for cls, v in f.items() if cls[1] == k) / z for k in range(5)]
    abl['dp5_no_configurational_degeneracy'] = {**oracle, 'vacancy_Zr_neighbour_distribution': nodeg}
    abl['dp5_lowest_energy_class_only'] = {**oracle, 'vacancy_Zr_neighbour_distribution': {k: [0, 0, 0, 0, 1.0] for k in ('x0.25', 'x0.75')}}
    # dp6: class-mean lookup that drops environment classes never purchased at that composition
    drop = {}
    for x in ('0.25', '0.75'):
        seen = {M.ENV[x][s] for s in bought[x]}
        f = {cls: (math.exp(-fn['local_x'](x, cls) / M.KT) if cls in seen else 0.0) for cls in M.CLASSES}
        pop = M.population(float(x), f)
        drop[LABEL[x]] = {'distribution': pop['distribution'], 'site_fraction': pop['site_fraction'] * math.exp(jensen)}
    abl['dp6_unpurchased_classes_dropped'] = {**oracle, 'vacancy_Zr_neighbour_distribution': {k: v['distribution'] for k, v in drop.items()},
                                              'vacancy_site_fraction': {k: v['site_fraction'] for k, v in drop.items()}}
    for k, v in abl.items():
        v.setdefault('lda_common_slope_model', oracle['lda_common_slope_model'])
    # dp7: CASTEP sites classified with the VASP x = 0.50 cell's environments (same index), not the authors' labels
    abl['dp7_castep_sites_given_vasp_cell_environments'] = {**oracle, 'lda_common_slope_model': lda_field(*common_slope(lda_items(lda_bought, 'vasp_cell')))}
    # dp8: functional transfer -- report the PBE description (VASP x = 0.50 table, or CASTEP-PBE) for LDA
    vasp_items = [(float(M.ENV['0.50'][s][0] == 3), M.ENV['0.50'][s][1], e) for s, e in M.EF['0.50'].items()]
    abl['dp8_pbe_description_reported_for_lda'] = {**oracle, 'lda_common_slope_model': lda_field(*common_slope(vasp_items))}
    abl['dp8_castep_pbe_reported_for_lda'] = {**oracle, 'lda_common_slope_model': lda_field(*common_slope(lda_items(P['castep']['PBE'])))}
    # ---------------- robustness (defensible alternatives)
    rob = {}
    rob['per_composition_offsets_shared_slopes'] = payload_parts(fn['local_offsets'], pops_from_class_energy(fn['local_offsets'], extra=jensen))
    off = {x: float(np.mean([M.EF[x][s] - f50(x, M.ENV[x][s]) for s in bought[x]])) for x in ('0.25', '0.75')}
    f2 = lambda x, cls: f50(x, cls) + off[x]
    rob['x050_slopes_plus_composition_offset'] = payload_parts(f2, pops_from_class_energy(f2, extra=jensen))
    rob['no_within_class_scatter_correction'] = payload_parts(fn['local_x'])
    cind = {x: M.fit(rows_for(x, bought[x]), 'local') for x in ('0.25', '0.75')}
    find = lambda x, cls: M.predict(cind[x], float(x), cls, 'local')
    rob['independent_per_composition_fits'] = payload_parts(find, pops_from_class_energy(find, extra=jensen))
    ccs = np.linalg.lstsq(np.array([[1.0, float(cls[0] == 3), cls[1], x] for x, cls, e in train]),
                          np.array([e for x, cls, e in train]), rcond=None)[0]
    fcs = lambda x, cls: float(np.dot(ccs, [1.0, float(cls[0] == 3), cls[1], float(x)]))
    rob['common_Zr_slope_both_sublattices'] = payload_parts(fcs, pops_from_class_energy(fcs, extra=jensen))
    # purchased sites' own Boltzmann averages for the classes bought, model elsewhere
    hyb = {}
    for x in ('0.25', '0.75'):
        f = {}
        for cls in M.CLASSES:
            es = [M.EF[x][s] for s in bought[x] if M.ENV[x][s] == cls]
            f[cls] = float(np.mean(np.exp(-np.array(es) / M.KT))) if es else math.exp(-fn['local_x'](x, cls) / M.KT + jensen)
        pop = M.population(float(x), f)
        hyb[LABEL[x]] = {'distribution': pop['distribution'], 'site_fraction': pop['site_fraction']}
    rob['lda_other_twelve_sites'] = {**oracle, 'lda_common_slope_model': lda_field(*common_slope(lda_items(
        [lda_all[s] for s in [32, 36, 39, 49, 54, 70, 44, 57, 62, 73, 79, 89]])))}
    rob['purchased_site_boltzmann_averages'] = {**oracle, 'vacancy_Zr_neighbour_distribution': {k: v['distribution'] for k, v in hyb.items()},
                                                'vacancy_site_fraction': {k: v['site_fraction'] for k, v in hyb.items()}}
    for k, v in rob.items():
        v.setdefault('lda_common_slope_model', oracle['lda_common_slope_model'])
    # ---------------- staged files
    agent = TASK/'environment/agent'
    (agent/'structures').mkdir(parents=True, exist_ok=True)
    for x in M.XS:
        (agent/'structures'/CIFNAME[x]).write_text(P['cif'][x])
    (agent/'structures/HZO_x0.50_3x3x3_Hf54Zr54O216.cif').write_text(P['cif_333'])
    buf = io.StringIO()
    wr = csv.DictWriter(buf, fieldnames=list(TABLE_COLUMNS), lineterminator='\n')
    wr.writeheader()
    for r in P['table_x050']:
        wr.writerow({k: r[k] for k in TABLE_COLUMNS})
    (agent/'x050_site_table.csv').write_text(buf.getvalue())
    sites_avail = {LABEL[x]: [s for s in sorted(M.EF[x]) if s not in HELD.get(x, [])] for x in M.XS}
    catalog = {
        'material': 'polar orthorhombic (Pca2_1) Hf1-xZrxO2 (x = Zr fraction); neutral oxygen vacancies; released by Taper, Cai, Strkalj, Monserrat and Di Martino (2026)',
        'supercells': {LABEL[x]: {'cif': 'structures/' + CIFNAME[x], 'formula': f, 'x': float(x)}
                       for x, f in zip(M.XS, ['Hf24Zr8O64', 'Hf16Zr16O64', 'Hf8Zr24O64'])},
        'site_numbering': 'Oxygen sites are identified by their 0-based atom index in the CIF (oxygens are atoms 32-95).',
        'vasp_protocol': 'VASP, PBE. Each vacancy energy is the total energy of the 96-atom supercell with one neutral O vacancy, '
                         'ionic positions relaxed at the lattice of the CIF; pristine energies are those of the same supercells without a vacancy.',
        'oxygen_chemical_potential_O_poor_eV': M.MU_POOR,
        'x050_site_table': {'file': 'x050_site_table.csv', 'cell': 'structures/' + CIFNAME['0.50'], 'columns': TABLE_COLUMNS,
                            'note': 'Released by the authors for the x = 0.50 cell (their table FIGS12). The release does not state the oxygen chemical potential of E_vac.'},
        'castep_sets': {'note': 'CASTEP vacancy calculations in 96-atom 2x2x2 supercells for seven exchange-correlation functionals. Only the final '
                                'enthalpy of each vacancy cell was released (no pristine-cell or reference energies). A purchase returns the enthalpy together with the '
                                'authors\' labels of that vacancy (vacancy type and local cation environment) and its fractional coordinates as released.',
                        'sites': {func: [{'site': r['site'], 'frac_coord': r['frac_coord']} for r in sorted(rows, key=lambda r: r['site'])]
                                  for func, rows in sorted(P['castep'].items())}},
        'castep_333': {'cif': 'structures/HZO_x0.50_3x3x3_Hf54Zr54O216.cif',
                       'note': 'CASTEP-PBE neutral vacancies in a 324-atom 3x3x3 supercell at x = 0.50 (final enthalpies only; no pristine-cell energy released). '
                               'Sites are 0-based atom indices of this CIF.',
                       'sites': sorted(int(s) for s in P['castep_333'])},
        'requests': {'vacancy': {'example': {'kind': 'vacancy', 'composition': 'x0.25', 'site': 40}, 'cost': COSTS['vacancy']},
                     'pristine': {'example': {'kind': 'pristine', 'composition': 'x0.25'}, 'cost': COSTS['pristine']},
                     'castep_vacancy': {'example': {'kind': 'castep_vacancy', 'functional': 'LDA', 'site': 57}, 'cost': COSTS['castep_vacancy']},
                     'castep_333_vacancy': {'example': {'kind': 'castep_333_vacancy', 'site': 112}, 'cost': COSTS['castep_333_vacancy']}},
        'budget': BUDGET,
        'vacancy_sites_available': sites_avail,
        'sites_to_predict': [{'id': f'{LABEL[x]}-O{s}', 'composition': LABEL[x], 'site': s} for x in HELD for s in HELD[x]],
    }
    (agent/'catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    lab = {'vacancy': {LABEL[x]: {str(s): P['vasp']['vacancy_eV'][x][str(s)] for s in sites_avail[LABEL[x]]} for x in M.XS},
           'pristine': {LABEL[x]: P['vasp']['bulk_eV'][x] for x in M.XS},
           'castep_vacancy': {func: {str(r['site']): r for r in rows} for func, rows in P['castep'].items()},
           'castep_333_vacancy': P['castep_333']}
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(site_targets.values())
    for x in ('x0.25', 'x0.75'):
        r = reference[x]
        values += [r['mean_nZr'], r['share_4Zr'], r['site_fraction'], r['log10_site_fraction']] + r['distribution']
    values += [reference['lda']['offset_O3_eV'], reference['lda']['slope_per_Zr_eV']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob, 'oracle_expected': oracle},
                                                                 indent=1) + '\n')
    # ---------------- calibration printout
    print('withheld classes', {x: [M.ENV[x][s] for s in HELD[x]] for x in HELD})
    print('oracle purchases', bought, 'sigma %.4f jensen(log10) %.4f' % (sigma, jensen / math.log(10)))
    print('model coefficients', {k: np.round(v, 4).tolist() for k, v in models.items()})
    print('reference', json.dumps({k: {kk: (round(vv, 4) if isinstance(vv, float) else [round(u, 4) for u in vv]) for kk, vv in v.items()} for k, v in reference.items()}))
    print('%-46s %8s %8s %8s %8s | %7s %7s %7s | %7s %7s %7s | %6s %6s' % ('variant', 'b25', 'b75', 'r25', 'r75', 'nZr25', 'p4_25', 'dlf25', 'nZr75', 'p4_75', 'dlf75', 'dD', 'db'))
    for name, parts in [('ORACLE', oracle)] + [('rob:' + k, v) for k, v in rob.items()] + [('abl:' + k, v) for k, v in abl.items()]:
        e = evaluate(parts, hidden)
        print('%-46s %+8.4f %+8.4f %8.4f %8.4f | %7.3f %7.3f %+7.3f | %7.3f %7.3f %+7.3f | %+6.3f %+6.3f' % (
            name, e['x0.25']['bias'], e['x0.75']['bias'], e['x0.25']['rmse'], e['x0.75']['rmse'],
            e['x0.25']['mean_nZr'], e['x0.25']['share_4Zr'], e['x0.25']['dlog10_sf'],
            e['x0.75']['mean_nZr'], e['x0.75']['share_4Zr'], e['x0.75']['dlog10_sf'], e['lda']['dD'], e['lda']['db']))


if __name__ == '__main__':
    main()
