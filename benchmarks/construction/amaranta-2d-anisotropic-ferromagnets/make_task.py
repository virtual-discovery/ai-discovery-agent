"""Write the catalog, lab model, private targets and ablation values of amaranta-2d-anisotropic-ferromagnets.

    cd benchmarks/construction/amaranta-2d-anisotropic-ferromagnets && python3 -B make_task.py

(run inside ai-theorist-benchmark-runtime:20261007; needs pool.json from extract.py and the Monte
Carlo results mc_CrI3_full.json / mc_CrI3_nn_only.json from mc_reference.py). Prints the
calibration table and writes calibration.json next to this script.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
import model as M
import reference as R
import tc_analysis as TA

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/amaranta-2d-anisotropic-ferromagnets'
NAMES = ['CrCl3', 'CrBr3', 'CrI3', 'CrGeTe3', 'CrSiSe3', 'CrSBr']
PREDICT = {'CrI3': (6, 6), 'CrGeTe3': (6, 6), 'CrSBr': (7, 5)}
BUDGET = 245
NOISE_MEV = 0.02
THRESHOLDS = {'polar_tol_deg': 10.0, 'anisotropy_abs_meV': 0.005, 'anisotropy_rel': 0.10,
              'bond_J_abs_meV': 0.10, 'bond_K_rel': 0.20, 'bond_Gamma_abs_meV': 0.06,
              'Tc_rel': 0.15, 'dm_abs_meV': 0.05, 'gap_rel': 0.15, 'prediction_rmse_meV_per_site': 0.03}
K_B = 0.08617333  # meV/K


def r4(x):
    return [round(float(v), 4) + 0.0 for v in x]


def seed(text):
    return int(hashlib.sha256(text.encode()).hexdigest()[:8], 16)


def angles(v):
    v = v / np.linalg.norm(v)
    th = math.degrees(math.acos(max(-1.0, min(1.0, v[2]))))
    ph = math.degrees(math.atan2(v[1], v[0])) % 360.0
    return [round(th, 2), round(ph, 2) % 360.0]


def from_angles(a):
    th, ph = np.radians(a[0]), np.radians(a[1])
    return np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])


def lab_noise(compound, n1, n2, spins):
    e = np.asarray(spins, float)
    e = e / np.linalg.norm(e, axis=1, keepdims=True)
    key = json.dumps([compound, n1, n2, np.round(e, 6).tolist()]).encode()
    return float(np.random.default_rng(int(hashlib.sha256(key).hexdigest()[:16], 16)).normal() * NOISE_MEV)


# ------------------------------------------------------------------ prediction configurations
def prediction_configs(c, name):
    n1, n2 = PREDICT[name]
    nb = len(c.basis)
    pos = np.array([b + i * c.cell[0] + j * c.cell[1] for b in c.basis for i in range(n1) for j in range(n2)])
    sub = np.repeat(np.arange(nb), n1 * n2)
    rng = np.random.default_rng(seed('predict-' + name))
    out = {}
    v = rng.normal(size=(len(pos), 3))
    out['P1'] = v / np.linalg.norm(v, axis=1, keepdims=True)
    if name == 'CrSBr':
        q = 2 * np.pi / (n1 * np.linalg.norm(c.cell[0]))
        ph = q * pos[:, 0] + 0.4
        out['P2'] = np.c_[np.sin(ph), np.zeros(len(pos)), np.cos(ph)]            # cycloid, xz plane
        out['P3'] = np.c_[-np.sin(ph), np.zeros(len(pos)), np.cos(ph)]           # opposite rotation sense
        s = np.where(sub == 0, 1.0, -1.0)
        out['P4'] = np.c_[s * np.sin(ph), np.zeros(len(pos)), np.cos(ph)]        # the two chains counter-rotating
    else:
        v = rng.normal(size=(len(pos), 3))
        out['P2'] = v / np.linalg.norm(v, axis=1, keepdims=True)
        B = R.bvec(c)
        qv = np.array([1 / n1, 0.0]) @ B
        ph = pos[:, :2] @ qv + 0.3
        a1 = c.cell[0] / np.linalg.norm(c.cell[0])
        out['P3'] = np.outer(np.sin(ph), a1) + np.outer(np.cos(ph), [0, 0, 1.0])  # cycloid in the (a1, z) plane
    res = {}
    for tag, conf in out.items():
        ang = [angles(s) for s in conf]
        res[f'{name}-{tag}'] = {'compound': name, 'supercell': [n1, n2], 'spin_angles_deg': ang,
                                'spins': np.array([from_angles(a) for a in ang])}
    return res


def per_site(c, table, K, n1, n2, spins):
    nb = len(c.basis)
    N = nb * n1 * n2
    fm = np.tile([0, 0, 1.0], (N, 1))
    return (M.supercell_energy(table, K, nb, n1, n2, spins) - M.supercell_energy(table, K, nb, n1, n2, fm)) / N


# ------------------------------------------------------------------ model variants (wrong paths)
def is_nn(c, k0, k, m1, m2):
    d = c.basis[k] + m1 * c.cell[0] + m2 * c.cell[1] - c.basis[k0]
    return abs(np.linalg.norm(d) / c.L[0] - 1) < M.CLASS_REL_TOL


def table_with(c, fn):
    """Bond table where each 1NN tensor J (bond k0 -> k, offset) is replaced by fn(k0, k, m1, m2, J)."""
    return [(k0, k, m1, m2, fn(k0, k, m1, m2, J) if is_nn(c, k0, k, m1, m2) else J) for k0, k, m1, m2, J in c.bond_table()]


def table_remapped(c, T0):
    """Teacher-model table with the archived-bond tensor replaced by T0 (symmetry-mapped as in the lab)."""
    c.symmetry_ops()
    saved = dict(c._cache)
    c._cache = {'ops': saved['ops'], 'Tsym': np.array(T0, float)}
    tab = c.bond_table()
    c._cache = saved
    return tab


def iso_table(c):
    return [(k0, k, m1, m2, np.eye(3) * np.trace(J) / 3) for k0, k, m1, m2, J in c.bond_table()]


def aniso_matrix_table(c, table, K, n1, n2, pattern):
    N = len(pattern)
    E = lambda n: M.supercell_energy(table, K, len(c.basis), n1, n2, np.outer(pattern, n)) / N
    ex = np.eye(3)
    Mm = np.zeros((3, 3))
    for a in range(3):
        Mm[a, a] = E(ex[a])
    for a in range(3):
        for b in range(a + 1, 3):
            Mm[a, b] = Mm[b, a] = (E((ex[a] + ex[b]) / np.sqrt(2)) - E((ex[a] - ex[b]) / np.sqrt(2))) / 2
    return Mm


def polar_deg(Mm):
    w, V = np.linalg.eigh(Mm)
    n = V[:, 0]
    return float(np.degrees(np.arccos(min(1.0, abs(n[2])))))


def load_mc(fname, window=None):
    f = HERE/fname
    if f.exists():
        return TA.analyse(json.loads(f.read_text()), window)
    if os.environ.get('PROVISIONAL'):
        print('WARNING: provisional Tc placeholder for', fname)
        v = 37.0 if 'full' in fname else 25.0
        return {'Tc_binder_largest_pair': v, 'Tc_binder_smallest_pair': v, 'Tc_chi_peak_largest_L': v, 'provisional': True}
    raise FileNotFoundError(f)


def main():
    pool = M.load_pool()
    C = {n: M.Compound(n, pool[n]) for n in NAMES}
    lab_model, catalog_compounds, predictions, to_predict, pconf = {}, {}, {}, {}, {}
    gs_type, polar, aniso, patterns = {}, {}, {}, {}
    calib = {'model_validation_vs_archived_DFT': {n: v for n, v in M.validate(pool).items() if n in NAMES},
             'model_validation_K_equals_A': {n: v for n, v in M.validate(pool, 1.0).items() if n in NAMES}}
    for n, c in C.items():
        rec = pool[n]
        tab = c.bond_table()
        first = rec['four_state_sets'][0]
        nmag0 = sum(s['species'] in M.MAGNETIC for s in first['supercell']['sites'])
        lab_model[n] = {'cell': c.cell.tolist(), 'basis': c.basis.tolist(), 'K': c.K, 'noise_meV': NOISE_MEV,
                        'e0_eV': round(first['energies_eV'][0] / nmag0, 6),
                        'bonds': [[k0, k, m1, m2, J.tolist()] for k0, k, m1, m2, J in tab]}
        st = rec['structure']
        mag_raw = [s['position'] for s in st['sites'] if s['species'] in M.MAGNETIC]
        catalog_compounds[n] = {
            'magnetic_species': 'Cr', 'cell_A': [r4(v) for v in st['cell']],
            'sites': [{'species': s['species'], 'position_A': r4(s['position'])} for s in st['sites']],
            'magnetic_sites': [r4(b) for b in mag_raw],
            'archive_spin_label': 'spin (magnetic moment): %.1f' % rec['spin_label'], 'archive_pk': rec['pk']}
        if n == 'CrSiTe3':
            catalog_compounds[n]['note'] = 'structure from the AiiDA provenance'
        # ground state: best collinear pattern of the isotropic part, confirmed by Luttinger-Tisza
        q, e_lt, vec, _ = R.luttinger_tisza(c)
        e_col, n1, n2, pat = R.best_collinear(c)
        assert abs(e_col - e_lt) < 1e-6 and n1 == n2 == 1, (n, e_col, e_lt)
        gs_type[n] = 'ferromagnet' if np.all(pat > 0) else 'Neel antiferromagnet'
        patterns[n] = pat
        Mm = R.anisotropy_matrix(c, 1, 1, pat)
        Ms = R.anisotropy_matrix(c, 1, 1, pat, 'sia')
        polar[n] = polar_deg(Mm)
        aniso[n] = {'total': float(Mm[0, 0] - Mm[2, 2]), 'single_ion': float(Ms[0, 0] - Ms[2, 2])}
        if n in PREDICT:
            for pid, p in prediction_configs(c, n).items():
                n1p, n2p = p['supercell']
                predictions[pid] = per_site(c, tab, c.K, n1p, n2p, p['spins'])
                to_predict[pid] = {k: v for k, v in p.items() if k != 'spins'}
                pconf[pid] = p
    # CrI3 bond exchange in the ideal cubic frame
    ci3 = C['CrI3']
    a, b = ci3.couples[0]
    e1, e2, e3, nrm = R.kitaev_frame(ci3, a, b)
    T_i3 = ci3.tensor(a, b - a)
    (bJ, bK, bG, bGp), res = R.jkgg(T_i3, e1, e2, e3)
    bond = {'J': bJ, 'K': bK, 'Gamma': bG, 'Gamma_prime': bGp}
    # CrSBr: DM vector of the bond site 0 -> site 0 + a1, and the zone-centre gap
    csb = C['CrSBr']
    T_sb = [J for k0, k, m1, m2, J in csb.bond_table() if (k0, k, m1, m2) == (0, 0, 1, 0)][0]
    D_ref = R.dm_vector(T_sb)
    Msb = R.anisotropy_matrix(csb, 1, 1, patterns['CrSBr'])
    gap_ref, wsb, Vsb = R.gap_uniform(Msb)
    gap_lin = R.gap_linearized(csb, csb.bond_table(), csb.K, [1, 1], Vsb[:, 0])
    assert abs(gap_lin[0] - gap_ref) < 1e-4
    # CrI3 ordering temperature from the reference Monte Carlo
    tc_full = load_mc('mc_CrI3_full.json', (33.5, 38.5))   # cubic fits of U(T) on 34-38 K around the crossing
    Tc_ref = tc_full['Tc_binder_largest_pair']

    # ------------------------------------------------------------------ wrong paths
    pids = {n: [p for p in predictions if p.startswith(n + '-')] for n in PREDICT}

    def preds_from(n, table, K):
        return {pid: per_site(C[n], table, K, *pconf[pid]['supercell'], pconf[pid]['spins']) for pid in pids[n]}

    abl, rob = {}, {}
    # dp1: CrSiSe3 assumed ferromagnetic like the others
    cse = C['CrSiSe3']
    Mfm = R.anisotropy_matrix(cse, 1, 1, np.array([1.0, 1.0]))
    abl['CrSiSe3_assumed_ferromagnet'] = {'type': 'ferromagnet', 'polar': polar_deg(Mfm),
                                          'anisotropy': {'total': float(Mfm[0, 0] - Mfm[2, 2]), 'single_ion': aniso['CrSiSe3']['single_ion']}}
    # dp2: whole ordered-state anisotropy attributed to the single-ion term (ferromagnet-only probes)
    abl['mae_as_single_ion'] = {'anisotropy': {n: {'total': aniso[n]['total'], 'single_ion': aniso[n]['total']} for n in NAMES},
                                'predictions': {}}
    for n in PREDICT:
        abl['mae_as_single_ion']['predictions'].update(preds_from(n, iso_table(C[n]), -aniso[n]['total']))
    # dp3: exchange anisotropy neglected (single-ion term decides the axis and the anisotropy)
    abl['single_ion_only'] = {'anisotropy': {n: {'total': aniso[n]['single_ion'], 'single_ion': aniso[n]['single_ion']} for n in NAMES},
                              'polar': {n: (0.0 if aniso[n]['single_ion'] > 0 else 90.0) for n in NAMES}, 'predictions': {}}
    for n in PREDICT:
        abl['single_ion_only']['predictions'].update(preds_from(n, iso_table(C[n]), C[n].K))
    # dp11: single-site four-state rotations (site 0 along +-z, +-x; the rest along y) in a 3x3 cell, at the noise level
    ss = {}
    for n in NAMES:
        c = C[n]
        n1, n2 = 3, 3
        N = len(c.basis) * n1 * n2
        Es = []
        for v in ([0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0]):
            sp = np.tile([0, 1.0, 0], (N, 1))
            sp[0] = v
            Es.append(M.supercell_energy(c.bond_table(), c.K, len(c.basis), n1, n2, sp) + lab_noise(n, n1, n2, sp))
        Kest = (Es[0] + Es[1] - Es[2] - Es[3]) / 2
        ss[n] = {'total': aniso[n]['total'], 'single_ion': -Kest}
    abl['single_site_four_state'] = {'anisotropy': ss}
    # dp4: CrI3 tensor of the measured bond copied to every bond (transpose for reversed bonds)
    Ti3 = ci3.symmetrized_T()
    tcopy = table_with(ci3, lambda k0, k, m1, m2, J: Ti3 if k0 == 0 else Ti3.T)
    Mc = aniso_matrix_table(ci3, tcopy, ci3.K, 1, 1, np.array([1.0, 1.0]))
    abl['CrI3_tensor_copied'] = {'polar': polar_deg(Mc), 'anisotropy': {'total': float(Mc[0, 0] - Mc[2, 2]), 'single_ion': aniso['CrI3']['single_ion']},
                                 'predictions': preds_from('CrI3', tcopy, ci3.K)}

    def rz(v, th):
        Rm = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        return Rm @ v
    # dp5: mirror-image octahedral frame (cubic axes on the wrong side of the bond)
    (mJ, mK, mG, mGp), _ = R.jkgg(T_i3, rz(e1, np.pi), rz(e2, np.pi), rz(e3, np.pi))
    abl['CrI3_mirror_frame'] = {'J': mJ, 'K': mK, 'Gamma': mG, 'Gamma_prime': mGp}
    # dp6: collinear-only mapping (off-diagonal elements never measured)
    Td = np.diag(np.diag(Ti3))
    (dJ, dK, dG, dGp), _ = R.jkgg(np.diag(np.diag(T_i3)), e1, e2, e3)
    abl['CrI3_diagonal_only'] = {'bond': {'J': dJ, 'K': dK, 'Gamma': dG, 'Gamma_prime': dGp},
                                 'predictions': preds_from('CrI3', table_remapped(ci3, Td), ci3.K)}
    # dp7: CrSBr off-diagonal xz read as symmetric (only J_xz measured) and applied to every 1NN bond
    Tsb = csb.symmetrized_T()
    Tsym = Tsb.copy()
    Tsym[2, 0] = Tsym[0, 2]
    tsym = table_with(csb, lambda k0, k, m1, m2, J: Tsym)
    Ms2 = aniso_matrix_table(csb, tsym, csb.K, 1, 1, np.array([1.0, 1.0]))
    abl['CrSBr_symmetric_offdiagonal'] = {'D': [0.0, 0.0, 0.0], 'polar': polar_deg(Ms2),
                                          'anisotropy': {'total': float(Ms2[0, 0] - Ms2[2, 2]), 'single_ion': aniso['CrSBr']['single_ion']},
                                          'gap': float(R.gap_uniform(Ms2)[0]), 'predictions': preds_from('CrSBr', tsym, csb.K)}
    # dp8: CrSBr DM vector of chain 0 copied to chain 1 (uniform instead of staggered)
    T0 = T_sb
    tuni = table_with(csb, lambda k0, k, m1, m2, J: T0 if m1 > 0 else T0.T)
    abl['CrSBr_uniform_dm'] = {'predictions': preds_from('CrSBr', tuni, csb.K)}
    # dp9: four-state mapping in a cell two lattice constants along a1: J(+a1) + J(-a1) is measured, D cancels
    Tal = (T0 + T0.T) / 2
    tal = table_with(csb, lambda k0, k, m1, m2, J: Tal)
    abl['CrSBr_aliased_2a'] = {'D': [0.0, 0.0, 0.0], 'predictions': preds_from('CrSBr', tal, csb.K)}
    # dp10: uniaxial gap from one in-plane axis (x = a)
    abl['CrSBr_uniaxial_gap'] = 2.0 / R.S_SPIN * float(Msb[0, 0] - Msb[2, 2])
    abl['CrSBr_uniaxial_gap_b'] = 2.0 / R.S_SPIN * float(Msb[1, 1] - Msb[2, 2])
    # dp11: ordering temperature from mean-field theory, or from Monte Carlo of the nearest-neighbour model only
    J0 = np.zeros((3, 3))
    for k0, k, m1, m2, J in ci3.bond_table():
        if k0 == 0:
            J0 += J
    abl['Tc_mean_field'] = TA.mean_field_tc(J0, ci3.K)   # classical single-site mean field with the full J0 and K
    tc_nn = load_mc('mc_CrI3_nn_only.json')
    abl['Tc_nn_only'] = tc_nn['Tc_binder_largest_pair']

    # robustness
    (aJ, aK, aG, aGp), _ = R.jkgg(T_i3, *actual_normal_frame(nrm))
    rob['CrI3_bond_actual_plane_normal'] = {'J': aJ, 'K': aK, 'Gamma': aG, 'Gamma_prime': aGp}
    rob['Tc_susceptibility_peak'] = tc_full['Tc_chi_peak_largest_L']
    rob['Tc_binder_small_pair'] = tc_full['Tc_binder_smallest_pair']
    # anisotropy read from noisy lab energies of collective states in 6x3 (36-site) cells
    rob_an = {}
    for n in NAMES:
        c = C[n]
        tab = c.bond_table()
        n1, n2 = 6, 3
        pat_sc = np.repeat(patterns[n], n1 * n2)
        def Elab(nv):
            sp = np.outer(pat_sc, nv)
            return M.supercell_energy(tab, c.K, len(c.basis), n1, n2, sp) + lab_noise(n, n1, n2, sp)
        # partner pattern flips the 1NN products (Neel <-> ferro for the honeycomb; chains AF along a for CrSBr)
        if n == 'CrSBr':
            alt = np.array([(-1) ** i for k in range(2) for i in range(n1) for j in range(n2)], float)
        else:
            alt = np.repeat(patterns[n] * np.array([1.0, -1.0]), n1 * n2)
        def Ealt(nv):
            sp = np.outer(alt, nv)
            return M.supercell_energy(tab, c.K, len(c.basis), n1, n2, sp) + lab_noise(n, n1, n2, sp)
        N = len(pat_sc)
        tot = (Elab([1, 0, 0]) - Elab([0, 0, 1])) / N
        other = (Ealt([1, 0, 0]) - Ealt([0, 0, 1])) / N
        rob_an[n] = {'total': tot, 'single_ion': (tot + other) / 2}
    rob['anisotropy_from_noisy_collective_states'] = rob_an
    # predictions from the archived (unsymmetrized) tensor instead of the stabilizer-averaged one
    rob['predictions_unsymmetrized_tensor'] = {}
    for n in PREDICT:
        rob['predictions_unsymmetrized_tensor'].update(preds_from(n, table_remapped_raw(C[n]), C[n].K))

    # ------------------------------------------------------------------ write files
    hidden = {'ground_states': gs_type, 'easy_axis_polar_deg': polar, 'anisotropy_meV_per_Cr': aniso,
              'CrI3_bond_exchange_meV': bond, 'CrI3_Tc_K': Tc_ref, 'CrSBr_dm_vector_meV': D_ref.tolist(),
              'CrSBr_spin_wave_gap_meV': gap_ref, 'predictions': predictions, 'thresholds': THRESHOLDS,
              'reference_notes': {'Tc_analysis': tc_full, 'CrI3_frame_e3': e3.tolist(), 'jkgg_residual': res}}
    catalog = {
        'source': 'A. Putatunda, G. Cuono, F. Orlando, S. Picozzi, "AMaRaNTA: automated first-principles exchange parameters in 2D magnets [Results]", '
                  'Materials Cloud Archive (2026), doi:10.24435/materialscloud:he-3k (CC BY 4.0); method paper arXiv:2509.20226.',
        'energies': 'Total energies (eV) of a supercell for classical spin configurations (each spin a unit vector). They are generated by the '
                    'authors\' AMaRaNTA spin model (bilinear exchange of limited range plus single-ion anisotropy), fitted by the authors to '
                    'their noncollinear LSDA+U calculations with spin-orbit coupling (VASP, PAW, U = 2.8 eV and J_H = 0.8 eV on Cr, '
                    'constrained moments). The model parameters are withheld. Deterministic noise of DFT-like size (standard deviation '
                    '%.2f meV per supercell energy) is added: an identical request returns identical energies. Energies contain an arbitrary '
                    'constant per magnetic site, the same for all supercells of a compound. Magnetic dipole-dipole interactions are not '
                    'included (as in the DFT).' % NOISE_MEV,
        'frame': 'Cartesian frame of each catalog cell (positions in A): x along a1, z normal to the layer.',
        'structures': 'Released monolayer structures as archived (cell rows a1, a2, c; all atoms, Cartesian A).',
        'supercells': 'Diagonal n1 x n2 supercells of the catalog cell with at most 36 magnetic sites. Site order: magnetic site k of the cell '
                      '(order of magnetic_sites), then i = 0..n1-1, then j = 0..n2-1; site position r = m_k + i a1 + j a2 (a1, a2 = first two '
                      'rows of cell_A). Spins are given in that order.',
        'request': {'kind': 'energies', 'compound': '<name>', 'supercell': [3, 3],
                    'spins': '<list of 1-6 configurations; each is a list of [sx, sy, sz], one per magnetic site>'},
        'cost': 'per configuration ceil(n_sites / 12) units (1-12 sites: 1; 13-24: 2; 25-36: 3)',
        'budget': BUDGET,
        'compounds': catalog_compounds,
        'to_predict': {pid: dict(p, definition='(E(configuration) - E(all spins along +z)) per magnetic site, in meV, same supercell')
                       for pid, p in to_predict.items()},
        'prediction_angles': 'spin_angles_deg = [theta, phi] per site in lab site order: theta from +z, phi from +x in the xy plane.'}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/private/model.json').write_text(json.dumps(lab_model) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1, default=float) + '\n')
    values = []
    for n in NAMES:
        values += [aniso[n]['total'], aniso[n]['single_ion']]
    values += [bJ, bK, bG, bGp, Tc_ref, gap_ref] + [x for x in D_ref.tolist() if abs(x) > 1e-3] + list(predictions.values())
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps([float(v) for v in values]) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1, default=float) + '\n')
    for f in ('environment/lab/public.json', 'tests/public.json'):
        p = json.loads((TASK/f).read_text())
        p['solver_budget'] = BUDGET
        if f.startswith('environment'):
            p['requests'] = {'energies': {'kind': 'energies', 'compound': '<name>', 'supercell': [3, 3], 'spins': '<1-6 configurations>',
                                          'cost': 'ceil(n_sites/12) per configuration'}}
            p['catalog'] = '/app/catalog.json'
        (TASK/f).write_text(json.dumps(p, indent=1) + '\n')

    # ------------------------------------------------------------------ calibration
    def rmse(n, pr):
        return math.sqrt(np.mean([(pr[p] - predictions[p]) ** 2 for p in pids[n]]))
    calib['targets'] = {'ground_states': gs_type, 'polar': polar, 'anisotropy': aniso, 'CrI3_bond': bond, 'CrI3_Tc_K': Tc_ref,
                        'CrSBr_D': D_ref.tolist(), 'CrSBr_gap': gap_ref, 'CrSBr_gap_eigs_rel': (wsb - wsb[0]).tolist(),
                        'predictions': predictions}
    calib['dipolar_mae_3muB'] = {n: R.dipolar_mae(C[n], 3.0, 150.0, patterns[n]) for n in NAMES}
    calib['ablations'] = {
        'CrSiSe3_assumed_ferromagnet': abl['CrSiSe3_assumed_ferromagnet'],
        'mae_as_single_ion': {'single_ion': {n: abl['mae_as_single_ion']['anisotropy'][n]['single_ion'] for n in NAMES},
                              'pred_rmse': {n: rmse(n, abl['mae_as_single_ion']['predictions']) for n in PREDICT}},
        'single_ion_only': {'total': {n: aniso[n]['single_ion'] for n in NAMES}, 'polar': abl['single_ion_only']['polar'],
                            'pred_rmse': {n: rmse(n, abl['single_ion_only']['predictions']) for n in PREDICT}},
        'single_site_four_state': {n: ss[n]['single_ion'] - aniso[n]['single_ion'] for n in NAMES},
        'CrI3_tensor_copied': {'polar': abl['CrI3_tensor_copied']['polar'], 'total': abl['CrI3_tensor_copied']['anisotropy']['total'],
                               'pred_rmse': rmse('CrI3', abl['CrI3_tensor_copied']['predictions'])},
        'CrI3_mirror_frame': abl['CrI3_mirror_frame'],
        'CrI3_diagonal_only': {'bond': abl['CrI3_diagonal_only']['bond'], 'pred_rmse': rmse('CrI3', abl['CrI3_diagonal_only']['predictions'])},
        'CrSBr_symmetric_offdiagonal': {k: v for k, v in abl['CrSBr_symmetric_offdiagonal'].items() if k != 'predictions'} |
                                       {'pred_rmse': rmse('CrSBr', abl['CrSBr_symmetric_offdiagonal']['predictions'])},
        'CrSBr_uniform_dm': {'pred_rmse': rmse('CrSBr', abl['CrSBr_uniform_dm']['predictions'])},
        'CrSBr_aliased_2a': {'pred_rmse': rmse('CrSBr', abl['CrSBr_aliased_2a']['predictions'])},
        'CrSBr_uniaxial_gap': [abl['CrSBr_uniaxial_gap'], abl['CrSBr_uniaxial_gap_b']],
        'Tc_mean_field': abl['Tc_mean_field'], 'Tc_nn_only': abl['Tc_nn_only'], 'Tc_nn_only_analysis': tc_nn}
    calib['robustness'] = {'CrI3_bond_actual_plane_normal': rob['CrI3_bond_actual_plane_normal'],
                           'Tc_susceptibility_peak': rob['Tc_susceptibility_peak'], 'Tc_binder_small_pair': rob['Tc_binder_small_pair'],
                           'anisotropy_from_noisy_collective_states': rob_an,
                           'predictions_unsymmetrized_tensor_rmse': {n: rmse(n, rob['predictions_unsymmetrized_tensor']) for n in PREDICT}}
    calib['thresholds'] = THRESHOLDS
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    print(json.dumps({k: calib[k] for k in ('targets', 'dipolar_mae_3muB', 'ablations', 'robustness')}, indent=1, default=float))


def actual_normal_frame(nrm):
    """Orthonormal frame with e3 = the actual Cr-L-Cr plane normal and e1 + e2 in the plane of e3 and z."""
    z = np.array([0, 0, 1.0])
    p = z - (z @ nrm) * nrm
    p /= np.linalg.norm(p)
    q = np.cross(nrm, p)
    return (p + q) / np.sqrt(2), (p - q) / np.sqrt(2), nrm


def table_remapped_raw(c):
    """Teacher-model table built from the archived tensor without stabilizer averaging (first matching operation)."""
    return table_remapped(c, c.T.copy())


if __name__ == '__main__':
    main()
