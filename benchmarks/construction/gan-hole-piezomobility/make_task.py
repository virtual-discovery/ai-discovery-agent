"""Generate the staged, lab, private and ablation data for gan-hole-piezomobility from pool.json (+ staged bands.json).

    python3 -B extract.py          # first: pool.json and environment/lab/bands.json
    python3 -B make_task.py        # catalog, lab data, private targets, ablation values, calibration printout
"""
import json
import numpy as np
from pathlib import Path
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gan-hole-piezomobility'
S = M.S
TEMPS = [100, 150, 200, 250, 300, 350, 400, 450, 500]
COSTS = {'bands': 1, 'mobility_serta': 1, 'mobility_ibte': 3}
BUDGET = 50
WITHHELD = ['compress_3%_e3', 'expand_3%_e3', 'expand_3%_e2', 'expand_3%_e4', 'compress_3%_e6']
GRADED = {'compress_3%_e3': ['zz'], 'expand_3%_e3': ['xx', 'yy', 'zz'], 'expand_3%_e2': ['xx', 'yy'],
          'expand_3%_e4': ['xx', 'yy', 'zz', 'yz'], 'compress_3%_e6': ['xx', 'yy', 'zz', 'xy']}
THRESHOLDS = {'diag_rel': 0.10, 'offdiag_rel': 0.15, 'window': {'compress_3%_e3': {'zz': [150.0, 280.0]}},
              'strain_abs': 0.002, 'best_mu_zz_rel': 0.05, 'best_inplane_rel': 0.08}
TWO_PCT = [f'{s}_2%_e{k}' for s in ('compress', 'expand') for k in range(1, 7)]


def r(x, n=6):
    return float(f'{x:.{n}g}')


def mobility_record(label, solver):
    st = S[label]
    out = {'state': label, 'solver': solver, 'units': 'cm^2/Vs', 'hole_density_cm3_nominal': 1e18,
           'tensor_frame': 'Cartesian frame of the catalog cell (x along a1, z along c); rows/columns x, y, z',
           'serta_drift_mobility': {str(T): st['serta'][str(T)]['mu'] for T in TEMPS},
           'fermi_level_eV': {str(T): st['serta'][str(T)]['fermi_eV'] for T in TEMPS}}
    if solver == 'ibte':
        meta = st['ibte_solver']
        out['ibte_drift_mobility'] = {str(T): st['ibte'][str(T)]['mu'] for T in TEMPS}
        out['ibte_hole_density_cm3'] = {str(T): st['ibte'][str(T)]['hole_density_cm3'] for T in TEMPS}
        out['ibte_solver_log'] = {'iterations': meta['iterations'], 'mob_maxiter': 1000,
                                  'last_max_change_cm2_per_Vs': meta['last_max_change_cm2_per_Vs'],
                                  'message': ('The iteration reached the maximum but did not converge.'
                                              if not meta['converged_flag'] else 'converged'),
                                  'crystal_symmetries_used': meta['crystal_symmetries']}
    return out


def main():
    labels = sorted(S, key=lambda s: (s != 'unstrained', s[-1], s))
    states = []
    for lab in labels:
        st = S[lab]
        states.append({'state': lab, 'cell_bohr': st['cell_bohr'], 'atoms_crystal': st['atoms_crystal'],
                       'bands_available': True, 'mobility_available': lab not in WITHHELD})
    catalog = {
        'material': 'wurtzite GaN (4-atom cell), released by Chen, Leveillee, Van de Walle and Giustino (2025)',
        'method': ('Quantum ESPRESSO 7.2 + EPW 5.7: fully relativistic norm-conserving LDA (ecutwfc 120 Ry), '
                   'GW quasiparticle energies read in place of the Kohn-Sham ones (eig_read), Wannier interpolation '
                   '(12 N-p Wannier functions) from 10x10x10 coarse k/q grids to 80x80x80 fine grids, polar long-range '
                   'terms, adaptive smearing, Fermi-surface window 0.2 eV below the valence-band maximum, hole density '
                   '1e18 cm^-3, T = 100-500 K. Each strained state was relaxed and recomputed from scratch.'),
        'labels': 'Archive directory names: unstrained, compress_X%_eK and expand_X%_eK (X = 1, 2, 3; K = 1..6).',
        'structures': 'cell_bohr rows are the lattice vectors a1, a2, a3 in bohr (x along a1 of the unstrained cell); '
                      'atoms_crystal are [species, crystal coordinates].',
        'temperatures_K': TEMPS,
        'bands': {'k_points': 'kgrid_crystal: the 1000 points of the 10x10x10 grid in crystal coordinates of each '
                              "state's own reciprocal lattice (same order as energies_eV rows)",
                  'band_indices': POOL_BANDS, 'note': 'bands 47-52 are the top six valence bands, 53-56 the lowest '
                  'four conduction bands; energies are absolute GW quasiparticle energies in eV as archived; vbm_eV is '
                  'the valence-band maximum the authors used as the Fermi-window reference'},
        'kgrid_crystal': M.K.round(6).tolist(),
        'requests': {
            'bands': {'example': {'kind': 'bands', 'state': 'unstrained'}, 'cost': COSTS['bands'],
                      'returns': 'GW quasiparticle energies on the 10x10x10 grid for one state'},
            'mobility': {'example': {'kind': 'mobility', 'state': 'compress_1%_e1', 'solver': 'ibte'},
                         'cost': {'serta': COSTS['mobility_serta'], 'ibte': COSTS['mobility_ibte']},
                         'returns': ('drift mobility tensors at every temperature in the self-energy relaxation-time '
                                     'approximation (serta), or the iterative Boltzmann solution together with its '
                                     'SERTA starting point and the solver log (ibte)')}},
        'budget': BUDGET,
        'states': states,
        'to_predict': WITHHELD}
    lab_data = {'costs': COSTS, 'withheld': WITHHELD,
                'mobility': {lab: {'serta': mobility_record(lab, 'serta'), 'ibte': mobility_record(lab, 'ibte')}
                             for lab in labels if lab not in WITHHELD},
                'vbm_eV': {lab: S[lab]['vbm_eV'] for lab in labels}}

    # ---- hidden targets -------------------------------------------------------------------------------------
    targets = {lab: M.as_dict(M.mu(lab)) for lab in WITHHELD}
    zz2 = {lab: M.mu(lab)[2, 2] for lab in TWO_PCT}
    best_zz = max(zz2, key=zz2.get)
    inpl2 = {lab: M.inplane_max(M.mu(lab)) for lab in TWO_PCT}
    best_in = max(inpl2, key=inpl2.get)
    assert best_zz == 'compress_2%_e3' and best_in.endswith('_2%_e6')
    e6_states = ['compress_2%_e6', 'expand_2%_e6']
    voigt = {lab: [round(v, 4) + 0.0 for v in S[lab]['strain_voigt_from_cell']] for lab in S}
    hidden = {'targets': targets, 'graded': GRADED, 'thresholds': THRESHOLDS,
              'best_mu_zz_2pct': {'state': best_zz, 'strain_voigt': voigt[best_zz], 'mu_zz': zz2[best_zz]},
              'best_inplane_2pct': {'states': e6_states, 'component_index': 5, 'abs_value': 0.04,
                                    'mu_ref': float(np.mean([inpl2[s] for s in e6_states])),
                                    'mu_values': {s: inpl2[s] for s in e6_states}},
              'two_pct_ranking': {'mu_zz': zz2, 'inplane_max': inpl2}}
    values = []
    for lab, keys in GRADED.items():
        values += [targets[lab][k] for k in keys]
    values += [zz2[best_zz], hidden['best_inplane_2pct']['mu_ref']] + [inpl2[s] for s in e6_states]

    # ---- reference (oracle-equivalent), robustness and ablations ---------------------------------------------
    mu = M.mu
    u = mu('unstrained')
    ref = {
        'compress_3%_e3': {**M.as_dict(M.geometric_next(*M.series(-1, 3))), 'zz': M.bracket_zz()[0]},
        'expand_3%_e3': M.as_dict(M.geometric_next(*M.series(1, 3))),
        'expand_3%_e2': M.as_dict(M.geometric_next(*M.series(1, 2))),
        'expand_3%_e4': M.as_dict(M.mirror(mu('compress_3%_e4'), 4)),
        'compress_3%_e6': M.as_dict(M.mirror(mu('expand_3%_e6'), 6))}

    def preds(**over):
        out = json.loads(json.dumps(ref))
        for lab, d in over.items():
            out[lab] = d
        return out

    def lin_tensor(lab):
        """Linear piezomobility from +-1 % central differences of the same component."""
        sign, pct, comp = (-1 if lab.startswith('compress') else 1), int(lab.split('_')[1][0]), int(lab[-1])
        slope = (mu(M.label(1, 1, comp)) - mu(M.label(-1, 1, comp))) / 2
        return M.as_dict(u + slope * sign * pct)

    def serta_scaled(lab):
        """SERTA tensor scaled by the unstrained IBTE/SERTA ratio (per diagonal component; mean for off-diagonals)."""
        ratio = np.diag(u) / np.diag(mu('unstrained', solver='serta'))
        scale = np.full((3, 3), float(ratio.mean()))
        np.fill_diagonal(scale, ratio)
        return mu(lab, solver='serta') * scale

    def flagged(lab):
        return not S[lab]['ibte_solver']['converged_flag']

    def mu_or_serta(lab):
        return serta_scaled(lab) if flagged(lab) else mu(lab)

    abl = {}
    # dp1: no mirror symmetry -> extrapolate the same-sign series for expand_3%_e4
    abl['W4_same_sign_series'] = preds(**{'expand_3%_e4': M.as_dict(M.linear_next(*M.series(1, 4)))})
    # dp1b: mirror partner copied without changing the sign of the strained off-diagonal element
    abl['shear_offdiag_sign_not_flipped'] = preds(**{'expand_3%_e4': M.as_dict(mu('compress_3%_e4')),
                                                    'compress_3%_e6': M.as_dict(mu('expand_3%_e6'))})
    # dp2: x<->y swap of the e1 partner at the same strain
    abl['W3_swap_from_e1'] = preds(**{'expand_3%_e2': M.as_dict(M.swap_xy(mu('expand_3%_e1')))})
    # dp3: linear piezomobility tensor from +-1 % central differences, applied to every target
    abl['linear_piezomobility_from_1pct'] = {lab: lin_tensor(lab) for lab in WITHHELD}
    # dp4: compress_3%_e3 extrapolated from the pre-crossing (0 -> 1 %) response
    lin01 = u + 3 * (mu('compress_1%_e3') - u)
    abl['W1_precrossing_extrapolation'] = preds(**{'compress_3%_e3': M.as_dict(lin01)})
    # dp5: SERTA survey of the 2 % states picks the state with the largest SERTA mu_zz
    zz2_serta = {lab: mu(lab, solver='serta')[2, 2] for lab in TWO_PCT}
    best_serta = max(zz2_serta, key=zz2_serta.get)
    abl['best_mu_zz_from_serta'] = {'strain_voigt': voigt[best_serta],
                                    'mu_zz_cm2_per_Vs': float(zz2_serta[best_serta] * u[2, 2] / mu('unstrained', solver='serta')[2, 2]),
                                    'state': best_serta}
    # dp6: in-plane maximum taken from the diagonal only
    diag2 = {lab: max(mu(lab)[0, 0], mu(lab)[1, 1]) for lab in TWO_PCT}
    best_diag = max(diag2, key=diag2.get)
    abl['best_inplane_diagonal_only'] = {'strain_voigt': voigt[best_diag], 'mu_inplane_max_cm2_per_Vs': float(diag2[best_diag]),
                                         'state': best_diag}
    # dp7: archive label read as the Voigt shear component
    oracle_in = {'strain_voigt': voigt['compress_2%_e6'], 'mu_inplane_max_cm2_per_Vs': inpl2['compress_2%_e6']}
    wrong = list(voigt['compress_2%_e6'])
    wrong[5] = -0.02
    abl['best_inplane_label_as_voigt'] = {**oracle_in, 'strain_voigt': wrong}
    # dp8: runs flagged 'did not converge' replaced by SERTA x unstrained IBTE/SERTA ratio
    abl['flagged_runs_replaced_by_scaled_serta'] = preds(**{
        'expand_3%_e3': M.as_dict(M.geometric_next(u, mu_or_serta('expand_1%_e3'), mu_or_serta('expand_2%_e3'))),
        'expand_3%_e4': M.as_dict(M.mirror(mu_or_serta('compress_3%_e4'), 4))})

    rob = {
        'W3_linear_series': preds(**{'expand_3%_e2': M.as_dict(M.linear_next(*M.series(1, 2)))}),
        'W5_linear_compress_series': preds(**{'compress_3%_e6': M.as_dict(M.linear_next(*M.series(-1, 6)))}),
        'W1_saturated_at_2pct': preds(**{'compress_3%_e3': {**ref['compress_3%_e3'], 'zz': M.bracket_zz()[1]}}),
        'W1_linear_1_to_2pct': preds(**{'compress_3%_e3': {**ref['compress_3%_e3'], 'zz': M.bracket_zz()[2]}}),
        'W2_flat_from_1pct': preds(**{'expand_3%_e3': M.as_dict(mu('expand_1%_e3'))}),
        'best_inplane_from_expand_e6': {'strain_voigt': voigt['expand_2%_e6'],
                                        'mu_inplane_max_cm2_per_Vs': inpl2['expand_2%_e6']}}
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio_e2e1 = [mu(M.label(1, p, 2)) / M.swap_xy(mu(M.label(1, p, 1))) for p in (1, 2)]
    r3 = ratio_e2e1[1] + (ratio_e2e1[1] - ratio_e2e1[0])
    rob['W3_anisotropy_ratio'] = preds(**{'expand_3%_e2': M.as_dict(np.nan_to_num(r3 * M.swap_xy(mu('expand_3%_e1'))))})

    # ---- write -----------------------------------------------------------------------------------------------
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab_data) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'reference': ref, 'ablations': abl, 'robustness': rob},
                                                                 indent=1) + '\n')
    for path in (TASK/'environment/lab/public.json', TASK/'tests/public.json'):
        pub = json.loads(path.read_text())
        pub['solver_budget'] = BUDGET
        if 'lab' in str(path):
            pub['requests'] = {'bands': {'kind': 'bands', 'state': '<state>', 'cost': 1},
                               'mobility': {'kind': 'mobility', 'state': '<state>', 'solver': 'serta | ibte',
                                            'cost': {'serta': 1, 'ibte': 3}}}
            pub['catalog'] = '/app/catalog.json'
        path.write_text(json.dumps(pub, indent=1) + '\n')

    # ---- calibration printout --------------------------------------------------------------------------------
    def errs(pred):
        out = {}
        for lab, keys in GRADED.items():
            for k in keys:
                out[f'{lab}:{k}'] = round(pred[lab][k] / targets[lab][k] - 1, 4)
        return out
    report = {'targets': {lab: {k: round(targets[lab][k], 2) for k in GRADED[lab]} for lab in GRADED},
              'reference_errors': errs(ref),
              'robustness_errors': {k: errs(v) for k, v in rob.items() if 'expand_3%_e2' in v},
              'robustness_inplane': rob['best_inplane_from_expand_e6'],
              'ablation_errors': {k: errs(v) for k, v in abl.items() if isinstance(v, dict) and 'expand_3%_e2' in v},
              'best_zz': hidden['best_mu_zz_2pct'], 'zz2_ranking': {k: round(v, 2) for k, v in sorted(zz2.items(), key=lambda x: -x[1])[:5]},
              'serta_zz2_ranking': {k: round(v, 2) for k, v in sorted(zz2_serta.items(), key=lambda x: -x[1])[:5]},
              'inplane2_ranking': {k: round(v, 2) for k, v in sorted(inpl2.items(), key=lambda x: -x[1])[:4]},
              'abl_best_serta': abl['best_mu_zz_from_serta'], 'abl_best_diag': abl['best_inplane_diagonal_only'],
              'bracket_zz': M.bracket_zz()}
    print(json.dumps(report, indent=1, default=float))


POOL_BANDS = list(range(M.POOL['band_indices'][0], M.POOL['band_indices'][1] + 1))

if __name__ == '__main__':
    main()
