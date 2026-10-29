"""Generate the staged, lab, private and ablation data for mg2x-midir-forbidden-transition.

    python3 -B benchmarks/construction/mg2x-midir-forbidden-transition/make_task.py

Reads the archive members extracted by fetch.py (model.ARCHIVE), asserts the archive facts the task rests on,
writes environment/agent/{catalog,experiment}.json, environment/lab/data.json, private/hidden*.json,
solution/ablation_values.json and the budget, and prints the calibration numbers (also saved to calibration.json).
"""
import json
import re
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mg2x-midir-forbidden-transition'
A = M.ARCHIVE
BUDGET = 32
COSTS = {'bands': 1, 'x_scan': 1, 'qsgw_box': 5, 'momentum_matrix': 1, 'tote': 1, 'questaal': 1,
         'pseudopotential_g0w0': 1, 'electron_phonon': 1}
MAX_PMAT_POINTS = 40
TOL = {'x3_minus_x1_eV': 0.010, 'vertex_meV': 8.0, 'ratio_rel': 0.04, 'predicted_eV': 0.008}
TOTE_COLUMNS = ('Columns: q (three, units of 2 pi / a), band, q index, the eigenvalue of the Hamiltonian that '
                'entered the calculation, the quasiparticle energy eQP (Z = 1), eQP without Z, and Z; energies in eV. '
                'The one-shot G0W0 values of the paper are eQP; in a QSGW run the file of iteration N holds, in the '
                'first energy column, the eigenvalues of the static QSGW Hamiltonian built from the self-energy of '
                'iteration N - 1.')
R5 = lambda x: round(float(x), 5)
R6 = lambda x: round(float(x), 6)


def sig(x, n=6):
    return float(f'{x:.{n}g}')


# ----------------------------------------------------------------------------- lab data
def bands_entry(c, level):
    p = M.COMPOUNDS[c]
    keep = list(range(p['nv'] - 2, p['nv'] + 5))
    names = [('L', 'Gamma'), ('Gamma', 'X'), ('X', 'W'), ('W', 'Gamma')]
    segs = []
    for (a, b), seg in zip(names, M.path_bands(c, level)):
        ref = seg[keep[0]]
        segs.append({'from': a, 'to': b, 'k_2pi_over_a': [[round(v, 6) for v in k] for _, k, _ in ref],
                     'E_eV': {f'band_{n}': [R5(e) for _, _, e in seg[n]] for n in keep}})
    return {'compound': c, 'level': level, 'valence_bands': p['nv'], 'bands_listed': keep,
            'energy_zero': 'valence-band maximum of this calculation', 'segments': segs}


def scan_entry(c, level, direction):
    p = M.COMPOUNDS[c]
    s = M.x_scan(c, level, direction)
    keep = list(range(p['nv'], p['nv'] + 4))
    q = [x for x, _ in s[keep[0]]]
    kvec = '(1 - q, 0, 0)' if direction == 'X-Gamma' else '(1, q, 0)'
    return {'compound': c, 'level': level, 'direction': direction, 'k_2pi_over_a': kvec,
            'q_from_X_2pi_over_a': q, 'valence_bands': p['nv'], 'energy_zero': 'valence-band maximum of this calculation',
            'E_eV': {f'band_{n}': [R6(e) for _, e in s[n]] for n in keep}}


def box_entry(c, name):
    p = M.COMPOUNDS[c]
    K, E = M.load_box(c, name, [p['lo'], p['hi']])
    d, h, r = M.BOXES[(c, name)]
    return {'compound': c, 'box': name, 'step_2pi_over_a': h, 'q_max_2pi_over_a': {'q_par': r[0], 'q_y': r[1], 'q_z': r[2]},
            'points': [[round(v, 4) for v in k] for k in K],
            'E_eV': {f"band_{b}": [R6(e) for e in E[b]] for b in (p['lo'], p['hi'])},
            'energy_zero': 'valence-band maximum of the QSGW calculation'}


TOTE_RUNS = {
    'g0w0_lda/Mg2Si': ('ecalj-supplementary/mg2si-g0w0-oneshot/TOTE.UP', 'Mg2Si',
                       'All-electron one-shot G0W0 from the LDA (ecalj), production settings: SCF 8^3, GW k mesh 6^3, room-temperature lattice.'),
    'g0w0_lda/Mg2Ge': ('ecalj-supplementary/mg2ge-g0w0-oneshot/TOTE.UP', 'Mg2Ge', 'As g0w0_lda/Mg2Si, for Mg2Ge.'),
    'g0w0_lda/Mg2Sn': ('ecalj-supplementary/mg2sn-g0w0-oneshot/TOTE.UP', 'Mg2Sn', 'As g0w0_lda/Mg2Si, for Mg2Sn.'),
    'g0w0_pbe/Mg2Si': ('ecalj-supplementary/mg2si-g0w0-oneshot-pbe/TOTE.UP', 'Mg2Si',
                       'All-electron one-shot G0W0 from the PBE (ecalj), otherwise as g0w0_lda.'),
    'g0w0_pbe/Mg2Ge': ('ecalj-supplementary/mg2ge-g0w0-oneshot-pbe/TOTE.UP', 'Mg2Ge', 'As g0w0_pbe/Mg2Si, for Mg2Ge.'),
    'qsgw_final_iteration/Mg2Si': ('ecalj-supplementary/qsgw-iterations/mg2si-666-a6338-scf8-cont/TOTE.UP.it10', 'Mg2Si',
                                   'Production all-electron QSGW run (ecalj, SCF 8^3, GW k mesh 6^3, room-temperature lattice a = 6.3380 A): TOTE file of iteration 10, the last.'),
    'qsgw_final_iteration/Mg2Ge': ('ecalj-supplementary/qsgw-iterations/mg2ge-666-a6393-scf8-cont/TOTE.UP.it10', 'Mg2Ge',
                                   'As qsgw_final_iteration/Mg2Si, for Mg2Ge (a = 6.3933 A).'),
    'qsgw80/Mg2Si': ('ecalj-supplementary/mg2si-666-a6338-scf8-cont-ss080/TOTE.UP', 'Mg2Si',
                     'QSGW80: QSGW with the self-energy scaled by 0.8, continued from the production run; TOTE file of the last iteration.'),
    'qsgw80/Mg2Ge': ('ecalj-supplementary/mg2ge-666-a6393-scf8-cont-ss080/TOTE.UP', 'Mg2Ge', 'As qsgw80/Mg2Si, for Mg2Ge.'),
    'qsgw80/Mg2Sn': ('ecalj-supplementary/mg2sn-666-a6760-scf8-cont-ss080/TOTE.UP', 'Mg2Sn', 'As qsgw80/Mg2Si, for Mg2Sn.'),
    'qsgw_low_temperature_lattice/Mg2Si': ('ecalj-supplementary/mg2si-666-a6322-scf8-cont/TOTE.UP', 'Mg2Si',
                                           'QSGW of Mg2Si as the production run but at a = 6.3218 A, the lattice constant at 100 K (the lowest temperature at which it was measured; Kasai et al. 2017 thermal expansion applied to the room-temperature 6.3380 A); TOTE file of the last iteration.'),
}
QUESTAAL_RUNS = {
    'prod/lda': None,
    'prod/qsgw': ('questaal/prod-k4-es/eigen-qsgw-final.eval.txt', 'Questaal QSGW of Mg2Si, production basis (LMTO with an empty sphere at the 4b site, no APW), GW k mesh 4^3, final self-energy.'),
    'prod/qsgw_ladder': ('questaal/prod-k4-es/eigen-bsw-final.eval.txt', 'As prod/qsgw, with ladder diagrams in the screened interaction W (lmgw --bsw, window of 4 valence and 4 conduction bands), final self-energy.'),
    'prod/qsgw_ladder_nc8': ('questaal/bsw-nc8/eigen-bsw-nc8-final.eval.txt', 'As prod/qsgw_ladder with 8 conduction bands in the ladder window.'),
    'prod/qsgw_ladder_nc12': ('questaal/bsw-nc12/eigen-bsw-nc12-final.eval.txt', 'As prod/qsgw_ladder with 12 conduction bands in the ladder window.'),
    'prod/lda_only': ('questaal/k4-es-lda/eigen-lda.eval.txt', 'LDA of Mg2Si in the production Questaal basis.'),
    'basis2/lda': ('questaal/basis2-k4-es/eigen-lda.eval.txt', 'LDA of Mg2Si in a larger Questaal LMTO basis (with the empty sphere).'),
    'basis2/qsgw': ('questaal/basis2-k4-es/eigen-qsgw-final.eval.txt', 'QSGW of Mg2Si in the larger basis, final self-energy.'),
    'basis2/qsgw_ladder': ('questaal/basis2-k4-es/eigen-bsw-final.eval.txt', 'QSGW with ladder diagrams in W in the larger basis (window 4/4), final self-energy.'),
}
QUESTAAL_RUNS.pop('prod/lda')
BGW_RUNS = {
    'Mg2Si/gpp_all_electron_density': ('mg2si', 'plasmon-pole model (Hybertsen-Louie GPP) with the charge density of all electrons'),
    'Mg2Si/gpp_valence_density': ('mg2si-valrho', 'GPP with the charge density of the valence bands only'),
    'Mg2Si/full_frequency': ('mg2si-ff', 'full frequency dependence (contour deformation), no plasmon-pole model'),
    'Mg2Ge/gpp_all_electron_density': ('mg2ge', 'GPP with the charge density of all electrons'),
    'Mg2Ge/gpp_valence_density': ('mg2ge-valrho', 'GPP with the charge density of the valence bands only'),
    'Mg2Ge/gpp_semicore_d_valence_density': ('mg2ge-dval', 'GPP with the charge density of the Ge 3d and valence bands'),
    'Mg2Ge/full_frequency': ('mg2ge-ff', 'full frequency dependence (contour deformation), no plasmon-pole model'),
}


def tote_entry(run):
    path, c, desc = TOTE_RUNS[run]
    p = M.COMPOUNDS[c]
    lines = (A/path).read_text().splitlines()
    rows = []
    for ln in lines[1:]:
        f = ln.split()
        if len(f) < 9:
            continue
        q = [float(x.replace('D', 'E')) for x in f[:3]]
        if (np.allclose(q, 0, atol=1e-6) or np.allclose(q, (0, 1, 0), atol=1e-6)) and p['nv'] - 1 <= int(f[3]) <= p['nv'] + 3:
            rows.append(ln.rstrip())
    return {'run': run, 'compound': c, 'description': desc, 'file': path.split('/')[-1], 'valence_bands': p['nv'],
            'first_line': lines[0].rstrip(), 'rows_at_Gamma_and_X': rows, 'column_note_from_archive_README': TOTE_COLUMNS}


def questaal_entry(run):
    path, desc = QUESTAAL_RUNS[run]
    L = M.questaal_levels(A/path)
    return {'run': run, 'compound': 'Mg2Si', 'description': desc, 'file': path.split('/')[-1], 'valence_bands': 4,
            'units': 'Ry (absolute scale of the calculation)', 'k_points_2pi_over_alat': {'Gamma': [0, 0, 0], 'X': [0, 1, 0]},
            'eigenvalues_Ry': {'Gamma': [round(e, 10) for e in L['Gamma'][:12]], 'X': [round(e, 10) for e in L['X'][:12]]}}


def bgw_entry(run):
    d, desc = BGW_RUNS[run]
    rows = M.bgw_rows(A/f'data/berkeleygw/{d}/eqp1.dat')
    c = run.split('/')[0]
    first = {'Mg2Si': 13, 'Mg2Ge': 22}[c]
    return {'run': run, 'compound': c, 'description': 'Pseudopotential G0W0 (BerkeleyGW 4.0) from the PBE (PseudoDojo pseudopotentials, Quantum ESPRESSO), 6^3 k mesh, about 900 bands, dielectric cutoff 30 Ry; ' + desc + '. Quasiparticle energies are Eqp1 (linearized Dyson equation) of eqp1.dat.',
            'first_conduction_band': first, 'k_points_crystal': {'Gamma': [0, 0, 0], 'X': [0, 0.5, 0.5]},
            'columns': ['band', 'E_DFT_eV', 'E_qp_eV'],
            'rows': {'Gamma': [list(r) for r in rows[(0.0, 0.0, 0.0)]], 'X': [list(r) for r in rows[(0.0, 0.5, 0.5)]]}}


def eph_entries():
    J = M.eph_json()
    ahc, zg = {}, {}
    names = {'Mg2Si': ('Mg2Si', 'Ge 3d not applicable'), 'Mg2Ge': ('Mg2Ge', 'core'), 'Mg2Ge-3d-valence': ('Mg2Ge', 'valence')}
    for key, (c, d3) in names.items():
        for grid, variants in J['ahc'][key].items():
            for T in (0, 300):
                states = {}
                for v, vals in variants.items():
                    if not v.startswith(f'T{T}_'):
                        continue
                    states[v.split('_', 1)[1]] = {s: [round(x, 4) for x in vals[s]] for s in ('G15v', 'X1', 'X3')}
                k = f"{c}|AHC|{grid}|{T}|{'valence' if d3 == 'valence' else 'core'}"
                ahc[k] = {'compound': c, 'method': 'AHC', 'q_grid': int(grid), 'temperature_K': T,
                          'ge_3d': None if c == 'Mg2Si' else d3,
                          'description': 'Allen-Heine-Cardona renormalization (Quantum ESPRESSO 7.5 ph.x + postahc.x), LDA, norm-conserving pseudopotentials; perturbed potential on a 6^3 grid Fourier-interpolated (with the dipole term) to the stated q grid. Keys eta*: non-adiabatic, broadening eta in eV; eta0.020_adiabatic: adiabatic. States: G15v = valence-band maximum at Gamma; X1, X3 = the two lowest conduction states at X (bands 5 and 6 of the LDA calculation).',
                          'columns_meV': ['Re[Total]', 'Upper_DW', 'Lower_DW', 'Upper_Fan', 'Re[Lower_Fan]', 'Im[Total]'],
                          'renormalization_meV': states}
    for key, val in J['zg'].items():
        m = re.match(r'(Mg2\w\w) (\d)x\dx\d (\d+) K', key)
        c, sc, T = m.group(1), int(m.group(2)), int(m.group(3))
        zg[f'{c}|ZG|{sc}|{T}'] = {'compound': c, 'method': 'ZG', 'supercell': sc, 'temperature_K': T,
                                  'description': 'Special-displacement (Zacharias-Giustino) supercell, Quantum ESPRESSO 7.3.1, LDA; band energies unfolded to the primitive cell; average of the special and the opposite configuration. Misses the long-range (Froehlich) part beyond the supercell.',
                                  'change_of_E_X3_minus_E_X1_meV': round(val['dE'], 4), 'change_of_gap_X1_minus_G15v_meV': round(val['gap'], 4)}
    return ahc, zg


# ----------------------------------------------------------------------------- staged experiment
def experiment():
    D = A/'data/digitized'
    k = json.loads((D/'koenig1961-band-curves-th140.json').read_text())
    lott = json.loads((D/'20260914-DATA-lott1966-mg2ge-curves.json').read_text())
    edge = json.loads((D/'20260916-DATA-koenig1961-edge.json').read_text())
    areas = json.loads((D/'20260914-DATA-koenig1961-areas.json').read_text())
    meV = lambda xs: [round(1e3 * x, 1) for x in xs]
    s4 = lambda xs: [sig(x, 4) for x in xs]
    return {
        'note': 'Points read from published figures by the archive author (data/digitized of the archive, CC0-1.0). Photon energies in meV. Only the curves are given; no peak readings.',
        'koenig1961_Mg2Si_band': {
            'source': 'Koenig, Lynch and Danielson, J. Phys. Chem. Solids 20, 122 (1961), Fig. 1, n-type Mg2Si sample 27B01; background A*E^-2.5+B subtracted by the archive author; vertical axis in relative units (pixels of the figure).',
            'carrier_concentration_cm3_by_T': {kk: sig(v, 4) for kk, v in areas['n_of_T'].items()},
            'curves': {T: {'E_meV': meV(c['E_eV']), 'absorbance_minus_background': s4(c['absorbance_minus_background'])}
                       for T, c in k['curves'].items()}},
        'lott1966_Mg2Ge_band': {
            'source': 'Lott and Lynch, Phys. Rev. 141, 681 (1966), n-type Mg2Ge; vertical axis as read from the figure above its background (relative units).',
            'curves': {T: {'E_meV': meV(c[0]), 'absorption': s4(c[1])} for T, c in lott['curves'].items()}},
        'koenig1961_interband_edges': {
            'source': 'Koenig, Lynch and Danielson (1961), Fig. 3: absorption coefficient (cm^-1) near the interband edge. The key MgSi_93K is the archive label of the Mg2Si curve at 93 K.',
            'curves': {T: {'E_meV': meV(c['E']), 'alpha_per_cm': s4(c['alpha'])} for T, c in edge['curves'].items()}},
    }


# ----------------------------------------------------------------------------- references and wrong paths
def references():
    ref, chk = {}, {}
    # X3 - X1 at X: LDA and QSGW from the fine scans (identical at X to the path and box files); G0W0 = column 7.
    split = {}
    for c, p in M.COMPOUNDS.items():
        split[c] = {}
        for lev in ('LDA', 'QSGW'):
            s = M.x_scan(c, lev, 'X-Gamma')
            split[c][lev] = s[p['x3']][0][1] - s[p['x1']][0][1]
            seg = M.path_bands(c, lev)[1]
            assert abs((seg[p['x3']][-1][2] - seg[p['x1']][-1][2]) - split[c][lev]) < 2e-6
        t = A/f"ecalj-supplementary/mg2{p['tag']}-g0w0-oneshot/TOTE.UP"
        split[c]['G0W0'] = M.tote_split(t, p['x1'], p['x3'], 6)
        split[c]['G0W0_col8_Z1'] = M.tote_split(t, p['x1'], p['x3'], 7)
        split[c]['G0W0_col6_input'] = M.tote_split(t, p['x1'], p['x3'], 5)
        # masses identify X1 / X3: X3 has the heavier longitudinal and the lighter transverse mass
        ml = {b: M.parabolic_mass(*zip(*M.x_scan(c, 'QSGW', 'X-Gamma')[b]), p['a']) for b in (p['x1'], p['x3'])}
        assert ml[p['x3']] > 1.7 * ml[p['x1']]
    q = json.loads((A/'data/quasiparticle/20260921-DATA-quasiparticle-path.json').read_text())['materials']
    tl = json.loads((A/'ecalj-supplementary/tote_levels.json').read_text())
    for c, p in M.COMPOUNDS.items():
        assert abs(1e3 * split[c]['LDA'] - q[c]['lda']['X1_X3_meV']) < 0.06 and abs(1e3 * split[c]['QSGW'] - q[c]['qsgw']['X1_X3_meV']) < 0.06
        assert abs(1e3 * split[c]['G0W0'] - tl[f"mg2{p['tag']}-g0w0-oneshot/TOTE.UP"]['X3mX1']) < 1e-6
        tote = A/f"ecalj-supplementary/mg2{p['tag']}-g0w0-oneshot/TOTE.UP"
        _, rows = M.tote_rows(tote)
        for r in rows:   # column 7 is the Z-renormalized linearized solution: (c7 - c6) = Z (c8 - c6)
            assert abs((r[6] - r[5]) - r[8] * (r[7] - r[5])) < 1e-6
    ref['x3_minus_x1_eV'] = {c: {lev: R6(split[c][lev]) for lev in ('LDA', 'G0W0', 'QSGW')} for c in M.COMPOUNDS}
    # QSGW variants (column 6 of the TOTE files of the last iteration)
    var = {}
    for run in ('qsgw80/Mg2Si', 'qsgw80/Mg2Ge', 'qsgw80/Mg2Sn', 'qsgw_low_temperature_lattice/Mg2Si',
                'qsgw_final_iteration/Mg2Si', 'qsgw_final_iteration/Mg2Ge', 'g0w0_pbe/Mg2Si', 'g0w0_pbe/Mg2Ge'):
        path, c, _ = TOTE_RUNS[run]
        p = M.COMPOUNDS[c]
        col = 6 if run.startswith('g0w0') else 5
        var[run] = M.tote_split(A/path, p['x1'], p['x3'], col)
        var[run + ':col7'] = M.tote_split(A/path, p['x1'], p['x3'], 6)
    assert abs(var['qsgw_final_iteration/Mg2Si'] - split['Mg2Si']['QSGW']) < 2e-4
    assert abs(var['qsgw_final_iteration/Mg2Ge'] - split['Mg2Ge']['QSGW']) < 2e-4
    # Questaal: vertex (ladder) correction within one code and basis
    ql = {r: M.questaal_levels(A/QUESTAAL_RUNS[r][0]) for r in QUESTAAL_RUNS}
    qs = {r: (L['X'][5] - L['X'][4]) * M.RY_EV * 1e3 for r, L in ql.items()}
    qgap = {r: (L['X'][4] - L['Gamma'][3]) * M.RY_EV * 1e3 for r, L in ql.items()}
    qj = json.loads((A/'questaal/questaal_levels.json').read_text())
    assert abs(qs['prod/qsgw'] - qj['prod-k4-es/eigen-qsgw-final.eval.txt']['X3mX1']) < 1e-6
    ref['vertex_correction_Mg2Si_meV'] = R6(qs['prod/qsgw_ladder'] - qs['prod/qsgw'])
    chk['vertex'] = {'nc8': qs['prod/qsgw_ladder_nc8'] - qs['prod/qsgw'], 'nc12': qs['prod/qsgw_ladder_nc12'] - qs['prod/qsgw'],
                     'basis2': qs['basis2/qsgw_ladder'] - qs['basis2/qsgw'],
                     'cross_code_questaal_ladder_minus_ecalj_qsgw': qs['prod/qsgw_ladder'] - 1e3 * split['Mg2Si']['QSGW'],
                     'qsgw80_minus_qsgw': 1e3 * (var['qsgw80/Mg2Si'] - split['Mg2Si']['QSGW']),
                     'gap_change_with_ladder': qgap['prod/qsgw_ladder'] - qgap['prod/qsgw'],
                     'questaal_qsgw_minus_ecalj_qsgw': qs['prod/qsgw'] - 1e3 * split['Mg2Si']['QSGW']}
    # IR / Raman ratio
    si_hot = M.nested('Mg2Si', ['coarse'], 5)
    si_cold = M.nested('Mg2Si', ['coarse', 'fine_0.01', 'fine_0.005'], 5)
    sn = M.nested('Mg2Sn', ['coarse'], 10)
    sn_x1 = M.nested('Mg2Sn', ['coarse'], 11)
    aSi, aSn = M.COMPOUNDS['Mg2Si']['a'], M.COMPOUNDS['Mg2Sn']['a']
    R = lambda hot, cold, n, a, **kw: M.fd_mean(hot, n, 300, a, **kw) / M.fd_mean(cold, n, 77, a, **kw)
    ir = json.loads((A/'analysis/ir_raman_ratio.json').read_text())
    archive_box = {p['n_cm3']: p['ratio'] for p in ir['points_box']}
    r16, r318 = R(si_hot, si_cold, 1e16, aSi), R(si_hot, si_cold, 3e18, aSi)
    assert abs(r16 - archive_box[1e16]) < 1e-9 and abs(r318 - archive_box[3e18]) < 1e-9
    assert abs(M.boltzmann_mean(si_hot, 300) / M.boltzmann_mean(si_cold, 77) - ir['boltzmann_ratio_box']) < 1e-9
    rsn = R(sn, sn, 1e16, aSn)
    ref['ir_raman_ratio_300K_over_77K'] = {'Mg2Si_1e16': R6(archive_box[1e16]), 'Mg2Si_3e18': R6(archive_box[3e18]),
                                           'Mg2Sn_1e16': R6(rsn)}
    # parabolic and line-scan wrong paths
    def masses(c, band):
        a = M.COMPOUNDS[c]['a']
        ml = M.parabolic_mass(*zip(*M.x_scan(c, 'QSGW', 'X-Gamma')[band]), a)
        mt = M.parabolic_mass(*zip(*M.x_scan(c, 'QSGW', 'X-W')[band]), a)
        return ml, mt, (ml * mt * mt) ** (1 / 3)
    _, _, mdos_si = masses('Mg2Si', 5)
    _, _, mdos_sn = masses('Mg2Sn', 10)
    em = json.loads((A/'data/mass-scan/effective-masses.json').read_text())['materials']
    assert abs(masses('Mg2Si', 5)[1] - em['Mg2Si']['qsgw']['X1']['parabolic']['m_t']) < 1e-4
    line_si = M.line_isotropic_ratio(*zip(*M.x_scan('Mg2Si', 'QSGW', 'X-W')[5]), aSi)
    assert abs(1 / line_si - em['Mg2Si']['qsgw']['X1']['thermal_average']['strength_ratio_77_over_300']) < 1e-3
    line_sn = M.line_isotropic_ratio(*zip(*M.x_scan('Mg2Sn', 'QSGW', 'X-W')[10]), aSn)
    q2all = lambda K: (K ** 2).sum(1)
    unit = lambda K: np.ones(len(K))

    def nohalf(box):
        K, e, w = box
        h = np.cbrt(np.abs(w).max())
        return K, e, np.where(w > 0, h ** 3, 0.0)
    abl_ratio = {
        'parabolic': {'Mg2Si_1e16': M.parabolic_ratio(1e16, mdos_si), 'Mg2Si_3e18': M.parabolic_ratio(3e18, mdos_si),
                      'Mg2Sn_1e16': M.parabolic_ratio(1e16, mdos_sn)},
        'xw_line_isotropic': {'Mg2Si_1e16': line_si, 'Mg2Si_3e18': line_si * M.parabolic_ratio(3e18, mdos_si) / M.parabolic_ratio(1e16, mdos_si),
                              'Mg2Sn_1e16': line_sn},
        'constant_matrix_element': {'Mg2Si_1e16': R(si_hot, si_cold, 1e16, aSi, weight=unit), 'Mg2Si_3e18': R(si_hot, si_cold, 3e18, aSi, weight=unit),
                                    'Mg2Sn_1e16': R(sn, sn, 1e16, aSn, weight=unit)},
        'isotropic_q2': {'Mg2Si_1e16': R(si_hot, si_cold, 1e16, aSi, weight=q2all), 'Mg2Si_3e18': R(si_hot, si_cold, 3e18, aSi, weight=q2all),
                         'Mg2Sn_1e16': R(sn, sn, 1e16, aSn, weight=q2all)},
        'boltzmann': {'Mg2Si_1e16': M.boltzmann_mean(si_hot, 300) / M.boltzmann_mean(si_cold, 77),
                      'Mg2Si_3e18': M.boltzmann_mean(si_hot, 300) / M.boltzmann_mean(si_cold, 77),
                      'Mg2Sn_1e16': M.boltzmann_mean(sn, 300) / M.boltzmann_mean(sn, 77)},
        'spin_not_counted': {'Mg2Si_1e16': R(si_hot, si_cold, 1e16, aSi, degeneracy=3), 'Mg2Si_3e18': R(si_hot, si_cold, 3e18, aSi, degeneracy=3),
                             'Mg2Sn_1e16': R(sn, sn, 1e16, aSn, degeneracy=3)},
        'one_octant_only': {'Mg2Si_1e16': R(si_hot, si_cold, 1e16, aSi, octants=1), 'Mg2Si_3e18': R(si_hot, si_cold, 3e18, aSi, octants=1),
                            'Mg2Sn_1e16': R(sn, sn, 1e16, aSn, octants=1)},
        'no_symmetry_plane_weights': {'Mg2Si_1e16': R(nohalf(si_hot), nohalf(si_hot), 1e16, aSi),
                                      'Mg2Si_3e18': R(nohalf(si_hot), nohalf(si_hot), 3e18, aSi),
                                      'Mg2Sn_1e16': R(nohalf(sn), nohalf(sn), 1e16, aSn)},
        'mg2si_value_for_mg2sn': {'Mg2Sn_1e16': archive_box[1e16]},
        'x1_occupied_in_mg2sn': {'Mg2Sn_1e16': R(sn_x1, sn_x1, 1e16, aSn)},
    }
    pm = {r['q']: r for r in M.pmat_rows()}
    K, e, w = si_hot
    P = np.array([pm[tuple(k)]['P2'] for k in K])
    bm = lambda T, W: float(((w * np.exp(-e / (M.KB_EV * T))) * W).sum() / (w * np.exp(-e / (M.KB_EV * T))).sum())
    abl_ratio['full_pbe_matrix_element'] = {'Mg2Si_1e16': bm(300, P) / bm(77, P)}
    def unfolded(box):     # octant unfolded to the full neighbourhood of X and summed plainly (outer faces at full weight)
        K, e, w = box
        h = np.cbrt(np.abs(w).max())
        w2 = np.full(len(K), h ** 3)
        for j in range(3):
            w2 *= np.where(np.abs(K[:, j]) < 1e-9, 0.5, 1.0)
        return K, e, w2
    rob_ratio = {'unfolded_plain_sum': {'Mg2Si_1e16': R(unfolded(si_hot), unfolded(si_hot), 1e16, aSi),
                                        'Mg2Si_3e18': R(unfolded(si_hot), unfolded(si_hot), 3e18, aSi),
                                        'Mg2Sn_1e16': R(unfolded(sn), unfolded(sn), 1e16, aSn)},
                 'coarse_box_both_T': {'Mg2Si_1e16': R(si_hot, si_hot, 1e16, aSi), 'Mg2Si_3e18': R(si_hot, si_hot, 3e18, aSi),
                                       'Mg2Sn_1e16': rsn},
                 'nested_boxes_both_T': {'Mg2Si_1e16': R(si_cold, si_cold, 1e16, aSi), 'Mg2Si_3e18': R(si_cold, si_cold, 3e18, aSi),
                                         'Mg2Sn_1e16': M.boltzmann_mean(sn, 300) / M.boltzmann_mean(sn, 77)}}
    # predicted splittings: QSGW + temperature corrections (lattice, electron-phonon)
    J = M.eph_json()
    dE = lambda key, grid, v: J['ahc'][key][grid][v]['dE'][0]
    zpr_si, th_si = dE('Mg2Si', '24', 'T0_eta0.010'), dE('Mg2Si', '24', 'T300_eta0.010')
    th_ge, zpr_ge = dE('Mg2Ge', '24', 'T300_eta0.010'), dE('Mg2Ge', '24', 'T0_eta0.010')
    for key in ('Mg2Si', 'Mg2Ge'):     # archived dE = X3 - X1 of the per-state values
        v = J['ahc'][key]['24']['T0_eta0.010']
        assert abs(v['dE'][0] - (v['X3'][0] - v['X1'][0])) < 1e-6
    lowT = var['qsgw_low_temperature_lattice/Mg2Si']
    qsi, qge = split['Mg2Si']['QSGW'], split['Mg2Ge']['QSGW']
    ref['predicted_splitting_eV'] = {'Mg2Si_T0': R6(lowT + 1e-3 * zpr_si), 'Mg2Si_300K': R6(qsi + 1e-3 * th_si),
                                     'Mg2Ge_300K': R6(qge + 1e-3 * th_ge)}
    zg = {k: v['dE'] for k, v in J['zg'].items()}
    rob_pred = {'zg_4x4x4': {'Mg2Si_T0': lowT + 1e-3 * zg['Mg2Si 4x4x4 0 K'], 'Mg2Si_300K': qsi + 1e-3 * zg['Mg2Si 4x4x4 300 K'],
                             'Mg2Ge_300K': qge + 1e-3 * zg['Mg2Ge 4x4x4 300 K']},
                'ahc_adiabatic': {'Mg2Si_T0': lowT + 1e-3 * dE('Mg2Si', '24', 'T0_eta0.020_adiabatic'),
                                  'Mg2Si_300K': qsi + 1e-3 * dE('Mg2Si', '24', 'T300_eta0.020_adiabatic'),
                                  'Mg2Ge_300K': qge + 1e-3 * dE('Mg2Ge', '24', 'T300_eta0.020_adiabatic')},
                'ahc_q16_or_ge3d': {'Mg2Si_T0': lowT + 1e-3 * dE('Mg2Si', '16', 'T0_eta0.005'),
                                    'Mg2Si_300K': qsi + 1e-3 * dE('Mg2Si', '16', 'T300_eta0.005'),
                                    'Mg2Ge_300K': qge + 1e-3 * dE('Mg2Ge-3d-valence', '24', 'T300_eta0.040')}}
    abl_pred = {'no_electron_phonon': {'Mg2Si_T0': lowT, 'Mg2Si_300K': qsi, 'Mg2Ge_300K': qge},
                'room_temperature_lattice_at_T0': {'Mg2Si_T0': qsi + 1e-3 * zpr_si},
                'zero_point_only_at_300K': {'Mg2Si_300K': qsi + 1e-3 * zpr_si, 'Mg2Ge_300K': qge + 1e-3 * zpr_ge},
                'electron_phonon_sign_reversed': {'Mg2Si_T0': lowT - 1e-3 * zpr_si, 'Mg2Si_300K': qsi - 1e-3 * th_si,
                                                  'Mg2Ge_300K': qge - 1e-3 * th_ge},
                'low_temperature_lattice_at_300K': {'Mg2Si_300K': lowT + 1e-3 * th_si},
                'qsgw_column7_low_temperature': {'Mg2Si_T0': var['qsgw_low_temperature_lattice/Mg2Si:col7'] + 1e-3 * zpr_si}}
    abl_split = {'energy_order_labels': {c: {lev: abs(split[c][lev]) for lev in ('LDA', 'G0W0', 'QSGW')} for c in M.COMPOUNDS},
                 'g0w0_z_equal_1_column': {c: {'G0W0': split[c]['G0W0_col8_Z1']} for c in M.COMPOUNDS},
                 'g0w0_input_eigenvalue_column': {c: {'G0W0': split[c]['G0W0_col6_input']} for c in M.COMPOUNDS},
                 'g0w0_from_pbe': {'Mg2Si': {'G0W0': var['g0w0_pbe/Mg2Si']}, 'Mg2Ge': {'G0W0': var['g0w0_pbe/Mg2Ge']}},
                 'qsgw80_as_qsgw': {c: {'QSGW': var[f'qsgw80/{c}']} for c in M.COMPOUNDS}}
    abl_vertex = {'cross_code': chk['vertex']['cross_code_questaal_ladder_minus_ecalj_qsgw'],
                  'qsgw80_scaling': chk['vertex']['qsgw80_minus_qsgw'], 'gap_change': chk['vertex']['gap_change_with_ladder']}
    extra = {'eph_dE_meV': {'Mg2Si_T0': zpr_si, 'Mg2Si_300K': th_si, 'Mg2Ge_300K': th_ge, 'Mg2Ge_T0': zpr_ge},
             'lattice_shift_meV': 1e3 * (lowT - qsi), 'questaal': qs, 'questaal_gap': qgap, 'variants': var, 'checks': chk,
             'masses': {'Mg2Si_X1': masses('Mg2Si', 5), 'Mg2Si_X3': masses('Mg2Si', 6), 'Mg2Sn_lower': masses('Mg2Sn', 10),
                        'Mg2Sn_upper': masses('Mg2Sn', 11)}}
    return ref, {'ratio': abl_ratio, 'pred': abl_pred, 'split': abl_split, 'vertex': abl_vertex}, \
        {'ratio': rob_ratio, 'pred': rob_pred, 'vertex': {k: chk['vertex'][k] for k in ('nc8', 'nc12', 'basis2')}}, extra


def flat(ref):
    out = []
    for c in M.COMPOUNDS:
        out += [ref['x3_minus_x1_eV'][c][lev] for lev in ('LDA', 'G0W0', 'QSGW')]
    out.append(ref['vertex_correction_Mg2Si_meV'])
    out += list(ref['ir_raman_ratio_300K_over_77K'].values()) + list(ref['predicted_splitting_eV'].values())
    return out


def main():
    ref, abl, rob, extra = references()
    catalog_items = {
        'bands': [{'compound': c, 'level': lev} for c in M.COMPOUNDS for lev in ('LDA', 'QSGW')],
        'x_scan': [{'compound': c, 'level': lev, 'direction': d} for c in M.COMPOUNDS for lev in ('LDA', 'QSGW') for d in ('X-Gamma', 'X-W')],
        'qsgw_box': [{'compound': c, 'box': b} for (c, b) in M.BOXES],
        'tote': [{'run': r, 'description': TOTE_RUNS[r][2]} for r in TOTE_RUNS],
        'questaal': [{'run': r, 'description': QUESTAAL_RUNS[r][1]} for r in QUESTAAL_RUNS],
        'pseudopotential_g0w0': [{'run': r, 'description': BGW_RUNS[r][1]} for r in BGW_RUNS],
    }
    ahc, zg = eph_entries()
    data = {
        'costs': COSTS, 'max_momentum_matrix_points': MAX_PMAT_POINTS,
        'bands': {f'{c}|{lev}': bands_entry(c, lev) for c in M.COMPOUNDS for lev in ('LDA', 'QSGW')},
        'x_scan': {f'{c}|{lev}|{d}': scan_entry(c, lev, d) for c in M.COMPOUNDS for lev in ('LDA', 'QSGW') for d in ('X-Gamma', 'X-W')},
        'qsgw_box': {f'{c}|{b}': box_entry(c, b) for (c, b) in M.BOXES},
        'momentum_matrix': {'|'.join(f'{v:.2f}' for v in r['q']): [r['E13'], r['E14'], r['P2']] for r in M.pmat_rows()},
        'tote': {r: tote_entry(r) for r in TOTE_RUNS},
        'questaal': {r: questaal_entry(r) for r in QUESTAAL_RUNS},
        'pseudopotential_g0w0': {r: bgw_entry(r) for r in BGW_RUNS},
        'electron_phonon': {**ahc, **zg},
    }
    eph_items = sorted({(v['compound'], v['method'], v.get('q_grid') or v.get('supercell'), v['temperature_K'], v.get('ge_3d')) for v in data['electron_phonon'].values()},
                       key=lambda t: (t[0], t[1], t[2], t[3], str(t[4])))
    catalog_items['electron_phonon'] = [
        ({'compound': c, 'method': m, 'q_grid': g, 'temperature_K': T} if m == 'AHC' else {'compound': c, 'method': m, 'supercell': g, 'temperature_K': T})
        | ({'ge_3d': d} if (m == 'AHC' and c == 'Mg2Ge') else {}) for c, m, g, T, d in eph_items]
    catalog = {
        'title': 'n-type Mg2Si, Mg2Ge, Mg2Sn: archived calculations of the conduction states at X (Sugawara 2026, Zenodo 23319151)',
        'codes': 'ecalj (all-electron LDA, one-shot G0W0, QSGW; PMT basis), Questaal 7.14.1 (QSGW and QSGW with ladder diagrams in W, Mg2Si), BerkeleyGW 4.0 + Quantum ESPRESSO 7.3.1 (pseudopotential PBE and G0W0; PBE momentum matrix elements), Quantum ESPRESSO 7.3.1/7.5 (electron-phonon, LDA). No spin-orbit coupling.',
        'structure': 'antifluorite (Fm-3m); X = (1, 0, 0) in units of 2 pi / a (three equivalent X valleys); W = (1, 1/2, 0).',
        'lattice_constants_A': {'Mg2Si': 6.3380, 'Mg2Ge': 6.3933, 'Mg2Sn': 6.7601,
                                'note': 'room-temperature values, used by the production calculations'},
        'valence_bands_ecalj': {c: p['nv'] for c, p in M.COMPOUNDS.items()},
        'requests': {
            'bands': {'cost': COSTS['bands'], 'format': {'kind': 'bands', 'compound': 'Mg2Si|Mg2Ge|Mg2Sn', 'level': 'LDA|QSGW'},
                      'returns': 'ecalj band energies along L-Gamma-X-W-Gamma (bands N_v-2 ... N_v+4), production calculations; QSGW at the final self-energy.'},
            'x_scan': {'cost': COSTS['x_scan'], 'format': {'kind': 'x_scan', 'compound': '...', 'level': 'LDA|QSGW', 'direction': 'X-Gamma|X-W'},
                       'returns': 'ecalj fine scan from X, 21 points with step 0.004 (2 pi / a), bands N_v ... N_v+3.'},
            'qsgw_box': {'cost': COSTS['qsgw_box'], 'format': {'kind': 'qsgw_box', 'compound': '...', 'box': 'coarse|fine_0.01|fine_0.005'},
                         'returns': 'QSGW (final self-energy) energies of bands N_v+1 and N_v+2 on a box of wave vectors on one side of X: q = k - X with 0 <= q_par <= q_max, 0 <= q_y, q_z <= q_max (one octant of the neighbourhood of the valley along x), on a uniform grid. coarse: step 0.02, q_par <= 0.24, q_y, q_z <= 0.16 (all three compounds); fine_0.01: step 0.01, 0.20 / 0.14; fine_0.005: step 0.005, 0.10 / 0.07 (Mg2Si only).'},
            'momentum_matrix': {'cost': COSTS['momentum_matrix'], 'format': {'kind': 'momentum_matrix', 'points': [[0.0, 0.0, 0.02]]},
                                'returns': f'Mg2Si only, pseudopotential PBE (Quantum ESPRESSO): energies of bands 13 and 14 (the two lowest conduction states at X; bands 1-12 occupied) and |<13|p|14>|^2 summed over the three Cartesian components (atomic units), at up to {MAX_PMAT_POINTS} points (q_par, q_y, q_z) of the coarse-box grid per request.'},
            'tote': {'cost': COSTS['tote'], 'format': {'kind': 'tote', 'run': '<run>'},
                     'returns': 'Rows of an ecalj TOTE file at Gamma and X (bands N_v-1 ... N_v+3), as written by the code, with its first line.'},
            'questaal': {'cost': COSTS['questaal'], 'format': {'kind': 'questaal', 'run': '<run>'},
                         'returns': 'Questaal eigenvalues (Ry) of bands 1-12 at Gamma and X for Mg2Si (band 4 is the valence-band maximum).'},
            'pseudopotential_g0w0': {'cost': COSTS['pseudopotential_g0w0'], 'format': {'kind': 'pseudopotential_g0w0', 'run': '<run>'},
                                     'returns': 'BerkeleyGW eqp1.dat rows at Gamma and X: band, E_DFT (PBE) and E_qp (eV).'},
            'electron_phonon': {'cost': COSTS['electron_phonon'],
                                'format': [{'kind': 'electron_phonon', 'compound': 'Mg2Si|Mg2Ge', 'method': 'AHC', 'q_grid': 24, 'temperature_K': 0, 'ge_3d': 'core|valence (Mg2Ge only)'},
                                           {'kind': 'electron_phonon', 'compound': 'Mg2Si|Mg2Ge', 'method': 'ZG', 'supercell': 4, 'temperature_K': 300}],
                                'returns': 'Electron-phonon renormalization of the band edges (meV), LDA with pseudopotentials: AHC per state (all broadenings), ZG as changes of X3 - X1 and of the gap.'},
        },
        'items': catalog_items,
        'tote_note_from_archive_README': TOTE_COLUMNS,
        'budget': BUDGET,
        'experiment': '/app/experiment.json',
    }
    hidden = {'targets': ref, 'tolerances': TOL}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/agent/experiment.json').write_text(json.dumps(experiment(), separators=(',', ':')) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(flat(ref)) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1, default=float) + '\n')
    for f in ('environment/lab/public.json', 'tests/public.json'):
        pub = json.loads((TASK/f).read_text())
        pub['solver_budget'] = BUDGET
        (TASK/f).write_text(json.dumps(pub, indent=1) + '\n')
    calib = {'reference': ref, 'ablations': abl, 'robustness': rob, 'extra': extra, 'tolerances': TOL}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    print(json.dumps(calib, indent=1, default=lambda x: round(float(x), 5)))


if __name__ == '__main__':
    main()
