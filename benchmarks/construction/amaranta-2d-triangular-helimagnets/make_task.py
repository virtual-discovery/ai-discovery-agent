"""Write the catalog, lab model, private targets and ablation values of amaranta-2d-triangular-helimagnets.

    cd benchmarks/construction/amaranta-2d-triangular-helimagnets && python3 -B make_task.py

(run inside ai-theorist-benchmark-runtime:20261007; needs pool.json from extract.py).
Prints the calibration table and writes calibration.json next to this script.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import model as M
import reference as R

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/amaranta-2d-triangular-helimagnets'
BUDGET = 240
NOISE_MEV = 0.02
THRESHOLDS = {'q_tol_inv_A': 0.02, 'prediction_rmse_meV_per_site': 0.3, 'axis_tol_deg': 7.0}
DISTORTED = ('VBr2', 'VI2', 'MnCl2')
FORMULA_NOTE = {'NiCl2': 'structure from the AiiDA provenance (input of AmarantaWorkChain 11815); the readable zip lacks this POSCAR'}


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


def prediction_configs(c, name):
    nb = len(c.basis)
    n1, n2 = (7, 7) if nb == 1 else (4, 7)
    pos = c.sites(n1, n2)
    rng = np.random.default_rng(seed('predict-' + name))
    q = (np.array([2 / 7, 1 / 7]) if nb == 1 else np.array([1 / 4, 2 / 7])) @ c.bvec()
    ph = pos[:, :2] @ q + 0.3
    spiral = np.c_[np.sin(ph), np.zeros(len(pos)), np.cos(ph)]
    v = rng.normal(size=(len(pos), 3))
    random_vec = v / np.linalg.norm(v, axis=1, keepdims=True)
    ising = np.outer(rng.choice([-1.0, 1.0], len(pos)), [0, 0, 1.0])
    out = {}
    for tag, conf in (('P1', spiral), ('P2', random_vec), ('P3', ising)):
        ang = [angles(s) for s in conf]
        out[f'{name}-{tag}'] = {'compound': name, 'supercell': [n1, n2], 'spin_angles_deg': ang,
                                'spins': np.array([from_angles(a) for a in ang])}
    return out


def per_site(table, K, c, n1, n2, spins):
    nb = len(c.basis)
    N = nb * n1 * n2
    fm = np.tile([0, 0, 1.0], (N, 1))
    return (R.energy(table, K, nb, n1, n2, spins) - R.energy(table, K, nb, n1, n2, fm)) / N


def lab_noise(compound, n1, n2, spins):
    e = np.asarray(spins, float)
    e = e / np.linalg.norm(e, axis=1, keepdims=True)
    key = json.dumps([compound, n1, n2, np.round(e, 6).tolist()]).encode()
    return float(np.random.default_rng(int(hashlib.sha256(key).hexdigest()[:16], 16)).normal() * NOISE_MEV)


def dist(c, d):
    return float(np.linalg.norm(d[:2]))


def modified_iso(c, rule):
    """Isotropic table where rule(d_vector_2d, length) returns J or None."""
    out = []
    for k0, b0 in enumerate(c.basis):
        for k, b in enumerate(c.basis):
            for m1 in range(-4, 5):
                for m2 in range(-4, 5):
                    d = (b + m1 * c.cell[0] + m2 * c.cell[1] - b0)
                    L = np.linalg.norm(d[:2])
                    if 0.1 < L < 8.6:
                        J = rule(d[:2], L)
                        if J is not None and J != 0:
                            out.append((k0, k, m1, m2, np.eye(3) * J))
    return out


def state(c, table):
    q, e = R.lt(table, c, grid=101)
    return {'type': R.classify(c, q), 'q_rlu': R.canonical_q(c, q)}, e


def commensurate_only(c, table):
    """Lowest of FM, 120-degree and collinear-stripe states (all symmetry equivalents)."""
    B = c.bvec()
    G = R.primitive_G(c)
    cands = []
    for h in np.arange(-1, 1.0001, 1 / 6):
        for k in np.arange(-1, 1.0001, 1 / 6):
            q = np.array([h, k])
            typ = R.classify(c, q)
            if typ != 'incommensurate spiral':
                cands.append((R.Jq(table, c, q @ B), typ, q))
    e, typ, q = min(cands, key=lambda x: x[0])
    return {'type': typ, 'q_rlu': R.canonical_q(c, q)}


def aniso_lt_q(c):
    """Generalised Luttinger-Tisza with the full exchange tensors (SIA omitted)."""
    tab = R.bond_table(c)
    B = c.bvec()
    def lam(x):
        q = np.array(x) @ B
        Mq = np.zeros((3, 3), complex)
        for k0, k, m1, m2, J in tab:
            if k0 == 0:
                d = (c.basis[k] + m1 * c.cell[0] + m2 * c.cell[1] - c.basis[0])[:2]
                Mq += 0.5 * J * np.exp(1j * q @ d)
        return float(np.linalg.eigvalsh(0.5 * (Mq + Mq.conj().T))[0])
    q0, _ = R.lt(R.bond_table(c, iso=True), c, grid=61)
    r = minimize(lam, q0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13})
    return r.x


def main():
    pool = M.load_pool()
    C = {n: M.Compound(n, pool[n]) for n in M.COMPOUNDS}
    lab_model, catalog_compounds, predictions, to_predict = {}, {}, {}, {}
    targets, equivalents, reciprocal = {}, {}, {}
    tables, iso_tables = {}, {}
    calib = {'model_validation_vs_archived_DFT': M.validate(pool)}
    for n, c in C.items():
        rec = pool[n]
        tab, tiso = R.bond_table(c), R.bond_table(c, iso=True)
        tables[n], iso_tables[n] = tab, tiso
        first = rec['four_state_sets'][0]
        nmag0 = sum(s['species'] in M.MAGNETIC for s in first['supercell']['sites'])
        lab_model[n] = {'cell': c.cell.tolist(), 'basis': c.basis.tolist(), 'K': c.K, 'noise_meV': NOISE_MEV,
                        'e0_eV': round(first['energies_eV'][0] / nmag0, 6),
                        'bonds': [[k0, k, m1, m2, J.tolist()] for k0, k, m1, m2, J in tab]}
        st = rec['structure']
        txt_spin = rec['spin_label']
        catalog_compounds[n] = {
            'magnetic_species': c.species, 'cell_A': [r4(v) for v in st['cell']],
            'sites': [{'species': s['species'], 'position_A': r4(s['position'])} for s in st['sites']],
            'magnetic_sites': [r4(b) for b in c.basis],
            'archive_spin_label': 'spin (magnetic moment): %.1f' % txt_spin,
            'archive_pk': rec['pk']}
        if n in FORMULA_NOTE:
            catalog_compounds[n]['note'] = FORMULA_NOTE[n]
        gs, e = state(c, tiso)
        targets[n] = {'type': gs['type'], 'q': gs['q_rlu'], 'E_gs_meV_per_site': e}
        equivalents[n] = {'q': np.round(R.equivalents(c, gs['q_rlu']), 6).tolist()}
        reciprocal[n] = c.bvec().tolist()
        for pid, p in prediction_configs(c, n).items():
            n1, n2 = p['supercell']
            predictions[pid] = per_site(tab, c.K, c, n1, n2, p['spins'])
            to_predict[pid] = {k: v for k, v in p.items() if k != 'spins'}
    vi2 = C['VI2']
    axis = R.easy_axis(tables['VI2'], vi2.K, 2, 1, 1, np.array([1.0, -1.0]))

    # ------------------------------------------------------------------ wrong paths
    pids = {n: [p for p in predictions if p.startswith(n + '-')] for n in C}
    def preds_from(n, table, K=0.0, factor=1.0):
        out = {}
        for pid in pids[n]:
            p = prediction_configs(C[n], n)[pid]
            out[pid] = factor * per_site(table, K, C[n], *p['supercell'], p['spins'])
        return out
    abl, rob = {}, {}
    # dp1: six equivalent nearest neighbours (ideal triangular shells, one representative bond per shell)
    gs1, pr1 = {}, {}
    for n in DISTORTED:
        c = C[n]
        a = np.linalg.norm(c.cell[0]); b = np.linalg.norm(c.cell[1])
        Jdiag = np.trace(c.T) / 3
        J_a = c.J3   # representative of the sqrt(3)-type shell: the bond along a
        nn = c.L[0]
        def rule(d, L, nn=nn, Jdiag=Jdiag, J_a=J_a):
            if L < nn * 1.05:
                return Jdiag
            if L < nn * 1.8:
                return J_a
            return None
        t = modified_iso(c, rule)
        gs1[n], _ = state(c, t)
        pr1.update(preds_from(n, t))
    abl['ideal_triangular'] = {'ground_states': gs1, 'predictions': pr1}
    # dp2 / dp3: third-neighbour coupling doubled by images, or omitted
    for key, f in (('J3_doubled', 2.0), ('no_J3', 0.0)):
        g, p = {}, {}
        for n, c in C.items():
            t = [(k0, k, m1, m2, J * (f if abs(np.linalg.norm((c.basis[k] + m1 * c.cell[0] + m2 * c.cell[1] - c.basis[k0])[:2]) - c.L[2]) < 0.01 else 1.0))
                 for k0, k, m1, m2, J in iso_tables[n]]
            t = [x for x in t if np.any(x[4])]
            g[n], _ = state(c, t)
            p.update(preds_from(n, t))
        abl[key] = {'ground_states': g, 'predictions': p}
    # dp4: compare only FM, 120-degree and stripe states
    abl['commensurate_only'] = {'ground_states': {n: commensurate_only(c, iso_tables[n]) for n, c in C.items()}}
    # dp5: easy axis from the three Cartesian axes (or from the single-ion term alone)
    E_axes = {ax: R.energy(tables['VI2'], vi2.K, 2, 1, 1, np.array([v, -v])) for ax, v in
              (('x', np.array([1.0, 0, 0])), ('y', np.array([0, 1.0, 0])), ('z', np.array([0, 0, 1.0])))}
    cart = {'x': [1.0, 0, 0], 'y': [0, 1.0, 0], 'z': [0, 0, 1.0]}[min(E_axes, key=E_axes.get)]
    abl['VI2_axis_cartesian_or_SIA_only'] = cart
    # dp6: the archived tensor copied to every nearest-neighbour bond (no site-symmetry transformation)
    def copied(d):
        return vi2.T if d[0] * vi2.d[0][0] + d[1] * vi2.d[0][1] > 0 else vi2.T.T
    tcopy = []
    for k0, k, m1, m2, J in tables['VI2']:
        d = vi2.basis[k] + m1 * vi2.cell[0] + m2 * vi2.cell[1] - vi2.basis[k0]
        tcopy.append((k0, k, m1, m2, copied(d) if abs(np.linalg.norm(d) - vi2.L[0]) < 0.01 else J))
    abl['VI2_axis_tensor_copied'] = R.easy_axis(tcopy, vi2.K, 2, 1, 1, np.array([1.0, -1.0])).tolist()
    # dp7: energies evaluated with every ordered pair counted (J from the pair-once mapping)
    abl['pairs_double_counted'] = {pid: v for n in C for pid, v in preds_from(n, iso_tables[n], factor=2.0).items()}
    # robustness
    rob['isotropic_predictions'] = {pid: v for n in C for pid, v in preds_from(n, iso_tables[n]).items()}
    rob['anisotropic_LT_q'] = {n: {'type': targets[n]['type'], 'q_rlu': [float(x) for x in aniso_lt_q(c)]} for n, c in C.items()}
    zz = {}
    for n, c in C.items():
        t = [(k0, k, m1, m2, (np.eye(3) * c.T[2, 2] if abs(np.linalg.norm((c.basis[k] + m1 * c.cell[0] + m2 * c.cell[1] - c.basis[k0])[:2]) - c.L[0]) < 0.01 else J))
             for k0, k, m1, m2, J in iso_tables[n]]
        zz[n], _ = state(c, t)
    rob['collinear_zz_mapping_q'] = zz
    # VI2 axis from six noisy probe directions in a 3x3 cell, through the lab's own noise
    probes = [np.array(v, float) / np.linalg.norm(v) for v in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 0], [0, 1, 1])]
    pat = np.array([1.0] * 9 + [-1.0] * 9)
    rows, ys = [], []
    for v in probes:
        spins = np.outer(pat, v)
        y = R.energy(tables['VI2'], vi2.K, 2, 3, 3, spins) + lab_noise('VI2', 3, 3, spins)
        rows.append([1, v[0] ** 2, v[1] ** 2, 2 * v[0] * v[1], 2 * v[0] * v[2], 2 * v[1] * v[2]])
        ys.append(y)
    cf = np.linalg.lstsq(np.array(rows), np.array(ys), rcond=None)[0]
    Mm = np.array([[cf[1], cf[3], cf[4]], [cf[3], cf[2], cf[5]], [cf[4], cf[5], 0.0]])
    w, V = np.linalg.eigh(Mm)
    rob['VI2_axis_six_probes'] = V[:, 0].tolist()

    # ------------------------------------------------------------------ write files
    hidden = {'ground_states': {n: {'type': t['type'], 'q': t['q']} for n, t in targets.items()},
              'equivalents': equivalents, 'reciprocal': reciprocal, 'predictions': predictions,
              'VI2_axis': axis.tolist(), 'thresholds': THRESHOLDS,
              'reference_notes': {'E_gs_meV_per_site': {n: t['E_gs_meV_per_site'] for n, t in targets.items()}}}
    catalog = {
        'source': 'A. Putatunda, G. Cuono, F. Orlando, S. Picozzi, "AMaRaNTA: automated first-principles exchange parameters in 2D magnets [Results]", '
                  'Materials Cloud Archive (2026), doi:10.24435/materialscloud:he-3k (CC BY 4.0); method paper arXiv:2509.20226.',
        'energies': 'Total energies (eV) of a supercell for classical spin configurations (each spin a unit vector). They are generated by the '
                    'authors\' AMaRaNTA spin model (bilinear exchange of limited range plus single-ion anisotropy), fitted by the authors to '
                    'their noncollinear LSDA+U calculations with spin-orbit coupling (VASP 6.4, PAW, U = 2.8 eV and J_H = 0.8 eV in the '
                    'Liechtenstein scheme, constrained moments). The model parameters are withheld. Deterministic noise of DFT-like size '
                    '(standard deviation %.2f meV per supercell energy) is added: an identical request returns identical energies. '
                    'Energies contain an arbitrary constant per magnetic site, the same for all supercells of a compound.' % NOISE_MEV,
        'frame': 'Cartesian frame of each catalog cell (positions in A): x along a1, z normal to the layer.',
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
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = []
    for n, t in targets.items():
        values += [x for x in t['q'] if abs(x) > 1e-9 and abs(abs(x) - 1) > 1e-9]
    values += list(predictions.values()) + [x for x in axis.tolist() if abs(x) > 1e-6]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    for f in ('environment/lab/public.json', 'tests/public.json'):
        p = json.loads((TASK/f).read_text())
        p['solver_budget'] = BUDGET
        if f.startswith('environment'):
            p['requests'] = {'energies': {'kind': 'energies', 'compound': '<name>', 'supercell': [3, 3], 'spins': '<1-6 configurations>',
                                          'cost': 'ceil(n_sites/12) per configuration'}}
            p['catalog'] = '/app/catalog.json'
        (TASK/f).write_text(json.dumps(p, indent=1) + '\n')

    # ------------------------------------------------------------------ calibration
    def qerr(n, q):
        return R.q_distance(C[n], q, targets[n]['q'])
    def rmse(n, pr):
        return math.sqrt(np.mean([(pr[p] - predictions[p]) ** 2 for p in pids[n]]))
    calib['targets'] = {n: {'type': t['type'], 'q': t['q']} for n, t in targets.items()}
    calib['VI2_axis'] = axis.tolist()
    calib['VI2_axis_tilt_from_c_deg'] = R.axis_angle_deg(axis, [0, 0, 1])
    calib['VI2_stripe_energy_meV_per_site'] = {k: v / 2 for k, v in E_axes.items()}
    calib['ablations'] = {
        'ideal_triangular': {n: {'gs': abl['ideal_triangular']['ground_states'][n], 'q_err': qerr(n, abl['ideal_triangular']['ground_states'][n]['q_rlu']),
                                 'pred_rmse': rmse(n, abl['ideal_triangular']['predictions'])} for n in DISTORTED},
        **{key: {n: {'gs': abl[key]['ground_states'][n], 'q_err': qerr(n, abl[key]['ground_states'][n]['q_rlu']),
                     'pred_rmse': rmse(n, abl[key]['predictions'])} for n in C} for key in ('J3_doubled', 'no_J3')},
        'commensurate_only': {n: {'gs': g, 'q_err': qerr(n, g['q_rlu'])} for n, g in abl['commensurate_only']['ground_states'].items()},
        'VI2_axis_cartesian_or_SIA_only': R.axis_angle_deg(cart, axis),
        'VI2_axis_tensor_copied': R.axis_angle_deg(abl['VI2_axis_tensor_copied'], axis),
        'pairs_double_counted': {n: rmse(n, abl['pairs_double_counted']) for n in C}}
    calib['robustness'] = {
        'isotropic_predictions_rmse': {n: rmse(n, rob['isotropic_predictions']) for n in C},
        'anisotropic_LT_q_err': {n: qerr(n, rob['anisotropic_LT_q'][n]['q_rlu']) for n in C},
        'collinear_zz_q_err': {n: qerr(n, zz[n]['q_rlu']) for n in C},
        'VI2_axis_six_probes_deg': R.axis_angle_deg(rob['VI2_axis_six_probes'], axis)}
    calib['thresholds'] = THRESHOLDS
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    print(json.dumps(calib, indent=1, default=float))


if __name__ == '__main__':
    main()
