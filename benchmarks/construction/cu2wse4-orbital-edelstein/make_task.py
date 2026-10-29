"""Write the staged, lab, private and ablation files of cu2wse4-orbital-edelstein.

    bash run_in_image.sh make_task.py

Needs /data/x (fetch.py), /data/model/*.npz (make_model.py) and /data/cache/{wosoc,soc}_N.npz (make_cache.py,
N = 16 24 32 48 64 96 128, CWS_DEGEN = 5e-3). Writes into benchmarks/tasks/cu2wse4-orbital-edelstein and
calibration.json here. All references use the lab's own evaluator (cws_tb.py) on the lab's own cache data.
"""
import gzip
import itertools
import json
import shutil
from pathlib import Path
import numpy as np
import cws_tb as T

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/cu2wse4-orbital-edelstein'
X, MODEL, CACHE = Path('/data/x'), Path('/data/model'), Path('/data/cache')
LAB, AGENT = TASK/'environment/lab', TASK/'environment/agent'
MESHES = [16, 24, 32, 48, 64, 96, 128]
NAMES = {'soc': 'soc', 'nosoc': 'wosoc'}
NVAL = {'soc': 88, 'nosoc': 44}
SPIN = {'soc': 1, 'nosoc': 2}
KT, GAMMA, EX, P_TARGET = 0.030, 0.030, 1.0e4, 5.0e20
MUB = 9.2740100783e-24
PREF = -EX / (2 * GAMMA) * MUB * 1e20        # M (A/m) per K (mu_B / Angstrom^2)
WINDOW = (-2.4, 2.4)                          # eV around each model's header Fermi level
BUDGET = 30
REF_MESH, ROB_MESH, COARSE_MESH = 128, 96, 32
THRESHOLDS = {'ratio_rel': 0.25, 'magnetization_rel': 0.15, 'soc_ratio_rel': 0.10, 'symmetry_abs': 0.05,
              'nlh_ratio_abs': 0.15}


def fermi(name):
    return float(np.load(MODEL/f'{NAMES[name]}.npz')['fermi_Ha']) * T.HARTREE


def full_cache(name, n):
    return dict(np.load(CACHE/f'{NAMES[name]}_{n}.npz'))


def compact(name, n):
    """Lab copy: window around the Fermi level, float32 histograms, states below the window in nlow."""
    c = full_cache(name, n)
    ef = fermi(name)
    lo = int(np.floor((ef + WINDOW[0] - float(c['emin'])) / T.EBIN))
    hi = int(np.ceil((ef + WINDOW[1] - float(c['emin'])) / T.EBIN))
    out = {'nk': np.array(float(c['nk'])), 'e0': np.array(float(c['emin']) + lo * T.EBIN),
           'nlow': np.array(float(c['nlow']) + float(c['count'][:lo].sum())), 'count': c['count'][lo:hi].astype(np.float32)}
    for q in T.QUANTS:
        if q == 'spin' and name == 'nosoc':
            continue
        out[q] = c[q][lo:hi].astype(np.float32)
    return out


def lab_cache(name, n):
    c = {k: np.asarray(v, float) for k, v in compact(name, n).items()}
    c['centres'] = c['e0'] + (np.arange(len(c['count'])) + 0.5) * T.EBIN
    return c


def K(c, name, q, mu):
    """Per-volume Fermi-surface tensor for the crystal (spin factor included), mu relative to the header."""
    if '+' in q:      # sum of moments, e.g. 'orb+atomic'
        return sum(K(c, name, part, mu) for part in q.split('+'))
    return SPIN[name] * T.surface_integral(c, q, fermi(name) + np.atleast_1d(mu), KT)


def holes(c, name, mu, spin=None):
    s = SPIN[name] if spin is None else spin
    x = np.clip((c['centres'][None, :] - fermi(name) - np.atleast_1d(mu)[:, None]) / KT, -300, 300)
    ne = (c['nlow'] + (1 / (np.exp(x) + 1)) @ c['count']) / c['nk']
    return s * (NVAL[name] - ne) / T.VOL * 1e24


def mu_at_p(c, name, vbm, p, spin=None):
    grid = np.arange(vbm - 0.3, vbm + 0.05, 0.0005)
    h = holes(c, name, grid, spin)
    return float(np.interp(np.log(p), np.log(h[::-1]), grid[::-1]))


def band_edges(name):
    from scipy.optimize import minimize
    m = T.Model(MODEL/f'{NAMES[name]}.npz')
    ef, nv = fermi(name), NVAL[name]
    band = lambda kf, n: m.bands(np.atleast_2d(kf) @ T.B)[:, n] - ef
    g = np.array([[i / 16, j / 16, l / 8] for i in range(-8, 8) for j in range(-8, 8) for l in range(-4, 4)])
    out = []
    for n, sgn in ((nv - 1, -1), (nv, 1)):
        E = band(g, n)
        seeds = list(g[np.argsort(sgn * E)[:8]]) + [np.array([0.06, 0.0, 0.0]), np.array([0.49, 0.49, 0.0])]
        best = min((minimize(lambda k: sgn * band(k, n)[0], s0, method='Nelder-Mead',
                             options={'xatol': 1e-7, 'fatol': 1e-9, 'maxiter': 4000}) for s0 in seeds), key=lambda r: r.fun)
        out.append((float(sgn * best.fun), best.x.round(5).tolist()))
    return out


def point_ops(atoms, rotate45=False):
    ops = []
    for perm in itertools.permutations(range(2)):
        for s in itertools.product((1, -1), repeat=3):
            g = np.zeros((3, 3)); g[0, perm[0]], g[1, perm[1]], g[2, 2] = s
            x0 = atoms[0][1]
            for sp, x in atoms:
                if sp != atoms[0][0]:
                    continue
                t = x - g @ x0
                if all(any(sp2 == sp1 and np.allclose(((g @ x1 + t - x2 + 0.5) % 1) - 0.5, 0, atol=1e-4) for sp2, x2 in atoms)
                       for sp1, x1 in atoms):
                    ops.append(g); break
    if rotate45:      # the same group in axes rotated by 45 degrees about z (twofold axes on the diagonals)
        r = np.array([[1, -1, 0], [1, 1, 0], [0, 0, np.sqrt(2)]]) / np.sqrt(2)
        ops = [r @ g @ r.T for g in ops]
    return ops


def allowed(ops):
    rng = np.random.default_rng(1)
    A = rng.normal(size=(3, 3)); B = rng.normal(size=(3, 3, 3)); B = B + B.transpose(0, 2, 1)
    S2 = sum(np.linalg.det(g) * g @ A @ g.T for g in ops) / len(ops)
    S3 = sum(np.einsum('ia,jb,lc,abc->ijl', g, g, g, B) for g in ops) / len(ops)
    ax = 'xyz'
    e = [ax[i] + ax[j] for i in range(3) for j in range(3) if abs(S2[i, j]) > 1e-8]
    n = [ax[i] + ax[j] + ax[l] for i in range(3) for j in range(3) for l in range(3) if abs(S3[i, j, l]) > 1e-8]
    return e, n, S2


def read_cif():
    atoms = []
    for line in (X/'Fig1a.cif').read_text().splitlines():
        w = line.split()
        if len(w) == 7 and w[0][:2] in ('Cu', 'W6', 'Se'):
            atoms.append((w[0].rstrip('+-0123456789'), np.array([float(v) for v in w[3:6]])))
    return atoms


def values(n_soc, n_nosoc, edges, moment='orb', alpha_spin_nosoc=2, count_spin_nosoc=2, nosoc_mu=None):
    s, ns = lab_cache('soc', n_soc), lab_cache('nosoc', n_nosoc)
    vbm, cbm = edges['soc'][0][0], edges['soc'][1][0]
    ko, ks = K(s, 'soc', moment, [vbm - 0.10, cbm + 0.10]), K(s, 'soc', 'spin', [vbm - 0.10, cbm + 0.10])
    ratio = {'VBM-0.10': float(ko[0, 0, 0] / ks[0, 0, 0]), 'CBM+0.10': float(ko[1, 0, 0] / ks[1, 0, 0])}
    mag = float(PREF * ko[0, 0, 0])
    mu_s = mu_at_p(s, 'soc', vbm, P_TARGET)
    mu_n = nosoc_mu(mu_s) if nosoc_mu else mu_at_p(ns, 'nosoc', edges['nosoc'][0][0], P_TARGET, spin=count_spin_nosoc)
    a_s = K(s, 'soc', moment, mu_s)[0, 0, 0]
    a_n = alpha_spin_nosoc / 2 * K(ns, 'nosoc', moment, mu_n)[0, 0, 0]
    return {'ratio': ratio, 'mag': mag, 'soc_nosoc': float(a_s / a_n), 'mu_soc_p': mu_s, 'mu_nosoc_p': mu_n}


def main():
    for d in (AGENT/'archive', LAB/'data', LAB/'private'):
        d.mkdir(parents=True, exist_ok=True)
    # ---- staged archive files (unchanged; band files gzip-compressed; HWR headers only)
    staged = ['Fig1a.cif'] + sorted(p.name for p in X.glob('Fig1c_*.d'))
    for f in staged:
        shutil.copy(X/f, AGENT/'archive'/f)
    for f in ('Fig1b_soc.d', 'Fig1b_wosoc.d'):
        with open(X/f, 'rb') as src, gzip.GzipFile(AGENT/'archive'/(f + '.gz'), 'wb', mtime=0) as dst:
            shutil.copyfileobj(src, dst)
    for f in ('Cu2WSe4_soc.HWR', 'Cu2WSe4_wosoc.HWR'):
        (AGENT/'archive'/(f + '.header')).write_text('\n'.join((X/f).read_text().splitlines()[:9]) + '\n')
    # ---- lab data: compact models and caches, evaluator
    for name in ('soc', 'nosoc'):
        z = np.load(MODEL/f'{NAMES[name]}.npz')
        np.savez_compressed(LAB/'data'/f'model_{name}.npz', R=z['R'].astype(np.int8), deg=z['deg'].astype(np.int8),
                            H_re=np.round(z['H'].real * 1e9).astype(np.int32), H_im=np.round(z['H'].imag * 1e9).astype(np.int32),
                            fermi_Ha=z['fermi_Ha'])
        for n in MESHES:
            np.savez_compressed(LAB/'data'/f'mesh_{name}_{n}.npz', **compact(name, n))
    shutil.copy(HERE/'cws_tb.py', LAB/'cws_tb.py')
    # ---- reference analysis
    edges = {name: band_edges(name) for name in ('soc', 'nosoc')}
    atoms = read_cif()
    e_ok, n_ok, S2 = allowed(point_ops(atoms))
    e_rot, n_rot, _ = allowed(point_ops(atoms, rotate45=True))
    ref = values(REF_MESH, REF_MESH, edges)
    rob = values(ROB_MESH, ROB_MESH, edges)
    # BCD ratio from the reference mesh (raw, kT = 10 meV, at CBM + 0.10 eV and VBM - 0.10 eV)
    s128 = lab_cache('soc', REF_MESH)
    bcd = {}
    for lab_, mu in (('VBM-0.10', edges['soc'][0][0] - 0.1), ('CBM+0.10', edges['soc'][1][0] + 0.1)):
        D = T.surface_integral(s128, 'bcd', fermi('soc') + np.array([mu]), 0.010)[0]
        bcd[lab_] = {'diag': np.diag(D).tolist(), 'raw_zxy_over_xyz': float((D[0, 0] - D[1, 1]) / (D[1, 1] - D[2, 2]))}
    # spin convergence: 1/N extrapolation of the valence-edge ratio (defensible alternative)
    sp = {}
    for n in (64, 96, 128):
        c = lab_cache('soc', n)
        mu = edges['soc'][0][0] - 0.10
        sp[n] = (K(c, 'soc', 'orb', mu)[0, 0, 0], K(c, 'soc', 'spin', mu)[0, 0, 0])
    slope = (sp[96][1] - sp[128][1]) / (1 / 96 - 1 / 128)
    spin_inf = sp[128][1] - slope / 128
    ratio_extrap = float(sp[128][0] / spin_inf)
    targets = {'orbital_to_spin_ratio': ref['ratio'], 'orbital_magnetization_A_per_m': ref['mag'],
               'soc_to_nosoc_orbital_ratio': ref['soc_nosoc'], 'edelstein_allowed_components': e_ok,
               'alpha_yy_over_alpha_xx': -1.0, 'chirality_indicator': 0.0,
               'nonlinear_hall_allowed_components': n_ok, 'sigma_zxy_over_sigma_xyz': -2.0}
    assert e_ok == ['xx', 'yy'] and abs(S2[1, 1] / S2[0, 0] + 1) < 1e-9 and abs(S2[2, 2]) < 1e-9
    assert sorted(n_ok) == sorted(['xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx'])
    # ---- ablations (realistic wrong paths, computed from the lab data)
    atomic = values(REF_MESH, REF_MESH, edges, moment='atomic')
    coarse = values(COARSE_MESH, COARSE_MESH, edges)
    summed = values(REF_MESH, REF_MESH, edges, moment='orb+atomic')
    abl = {
        'rotated_axes_setting': {'edelstein_allowed_components': e_rot, 'nonlinear_hall_allowed_components': n_rot},
        'chiral_d2_assumption': {'edelstein_allowed_components': ['xx', 'yy', 'zz']},
        'intra_atomic_orbital_moment': {'orbital_to_spin_ratio': atomic['ratio'], 'orbital_magnetization_A_per_m': atomic['mag'],
                                        'soc_to_nosoc_orbital_ratio': atomic['soc_nosoc']},
        'berry_plus_intra_atomic_moment': {'orbital_to_spin_ratio': summed['ratio'], 'orbital_magnetization_A_per_m': summed['mag'],
                                           'soc_to_nosoc_orbital_ratio': summed['soc_nosoc']},
        'no_spin_factor_in_nosoc_alpha': {'soc_to_nosoc_orbital_ratio': values(REF_MESH, REF_MESH, edges, alpha_spin_nosoc=1)['soc_nosoc']},
        'no_spin_factor_in_nosoc_holes': {'soc_to_nosoc_orbital_ratio': values(REF_MESH, REF_MESH, edges, count_spin_nosoc=1)['soc_nosoc']},
        'no_spin_factor_anywhere_nosoc': {'soc_to_nosoc_orbital_ratio': values(REF_MESH, REF_MESH, edges, alpha_spin_nosoc=1,
                                                                               count_spin_nosoc=1)['soc_nosoc']},
        'common_header_energy_axis': {'soc_to_nosoc_orbital_ratio': values(REF_MESH, REF_MESH, edges, nosoc_mu=lambda m: m)['soc_nosoc']},
        'coarse_mesh_32': {'orbital_to_spin_ratio': coarse['ratio'], 'orbital_magnetization_A_per_m': coarse['mag'],
                           'soc_to_nosoc_orbital_ratio': coarse['soc_nosoc']},
        'tau_hbar_over_gamma': {'orbital_magnetization_A_per_m': 2 * ref['mag']},
        'positive_charge_sign': {'orbital_magnetization_A_per_m': -ref['mag']},
        'conduction_ratio_from_valence': {'orbital_to_spin_ratio': {'VBM-0.10': ref['ratio']['VBM-0.10'],
                                                                    'CBM+0.10': ref['ratio']['VBM-0.10']}},
        'nlh_full_permutation_symmetry': {'sigma_zxy_over_sigma_xyz': 1.0},
        'chirality_from_chiral_class': {'chirality_indicator': 2.0 + float(K(s128, 'soc', 'orb', edges['soc'][0][0] - 0.1)[0, 2, 2]
                                                                         / K(s128, 'soc', 'orb', edges['soc'][0][0] - 0.1)[0, 0, 0]),
                                        'alpha_yy_over_alpha_xx': 1.0},
    }
    robust = {'mesh_96': {'orbital_to_spin_ratio': rob['ratio'], 'orbital_magnetization_A_per_m': rob['mag'],
                          'soc_to_nosoc_orbital_ratio': rob['soc_nosoc']},
              'valence_spin_extrapolated_1_over_N': {'orbital_to_spin_ratio': {'VBM-0.10': ratio_extrap,
                                                                                'CBM+0.10': ref['ratio']['CBM+0.10']}},
              'raw_bcd_ratio_mesh128': {'sigma_zxy_over_sigma_xyz': bcd['CBM+0.10']['raw_zxy_over_xyz']}}
    # ---- private and solution files
    hidden = {'targets': targets, 'thresholds': THRESHOLDS,
              'reference': {'edges_eV': edges, 'mu_soc_at_p': ref['mu_soc_p'], 'mu_nosoc_at_p': ref['mu_nosoc_p'],
                            'mesh': REF_MESH, 'kT_eV': KT, 'gamma_eV': GAMMA, 'p_cm3': P_TARGET}}
    (LAB/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    hv = list(ref['ratio'].values()) + [ref['mag'], ref['soc_nosoc']]
    (LAB/'private/hidden_values.json').write_text(json.dumps(hv) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': robust}, indent=1) + '\n')
    # ---- per-mesh table and calibration
    per_mesh = {}
    for n in MESHES:
        v = values(n, n, edges)
        per_mesh[n] = {'ratio': v['ratio'], 'mag': v['mag'], 'soc_nosoc': v['soc_nosoc']}
    cal = {'targets': targets, 'robustness': robust, 'ablations': abl, 'per_mesh': per_mesh, 'edges': edges,
           'bcd_mesh128': bcd, 'thresholds': THRESHOLDS, 'mu_at_p': {'soc': ref['mu_soc_p'], 'nosoc': ref['mu_nosoc_p']}}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1) + '\n')
    # ---- catalog and budgets
    catalog = {
        'material': 'Cu2WSe4, space group I-42m (No. 121), conventional cell Cu4W2Se8; OpenMX GGA (PBE), fully relativistic '
                    'with spin-orbit coupling or scalar-relativistic without',
        'source': 'Zenodo record 20131063 (DOI 10.5281/zenodo.20131063, CC BY 4.0), CWS_PRM.zip, data for K. Nakazawa, '
                  'T. Yamaguchi, A. Yamakage, Phys. Rev. Materials 10, 064603 (2026)',
        'staged_files': {'root': '/app/archive', 'files': sorted(p.name for p in (AGENT/'archive').iterdir()),
                         'note': 'archive files unchanged (the two Fig1b band files gzip-compressed); *.HWR.header are the '
                                 'first 9 lines of each Wannier Hamiltonian file; the archive gives no column descriptions. '
                                 'Not staged: the Hamiltonians (evaluated by the service) and the response data Fig2-Fig4 (withheld).'},
        'wannier_models': {
            'soc': {'file': 'Cu2WSe4_soc.HWR', 'num_wann': 108, 'basis': 'Cu d (4 atoms), W d (2 atoms), Se p (8 atoms), atom order of '
                    'Fig1a.cif, OpenMX real harmonics d = (z2, x2-y2, xy, xz, yz), p = (x, y, z); functions 1-54 spin up, '
                    '55-108 spin down (quantised along z)'},
            'nosoc': {'file': 'Cu2WSe4_wosoc.HWR', 'num_wann': 54, 'basis': 'the same 54 orbitals, without spin'}},
        'requests': {
            'bands': {'format': {'kind': 'bands', 'model': 'soc | nosoc', 'k_frac or k_cart_invA': '[[x, y, z], ...] (<= 1000)'}, 'cost': 1},
            'band_properties': {'format': {'kind': 'band_properties', 'model': 'soc | nosoc',
                                           'k_frac or k_cart_invA': '[[x, y, z], ...] (<= 200)'}, 'cost': 1},
            'fermi_surface': {'format': {'kind': 'fermi_surface', 'model': 'soc | nosoc', 'mesh': 'N, one of 16, 24, 32, 48, 64, 96, 128 '
                                         '(N x N x N/2 mesh)', 'mu_eV': '[...] (1-400, -1.9 to 1.9 eV)', 'kT_eV': '0.005-0.05'},
                              'cost': {'soc': {str(k): v for k, v in {16: 1, 24: 1, 32: 1, 48: 2, 64: 3, 96: 6, 128: 10}.items()},
                                       'nosoc': {str(k): v for k, v in {16: 1, 24: 1, 32: 1, 48: 1, 64: 2, 96: 3, 128: 5}.items()}}}},
        'budget': BUDGET,
        'request_size_limit': 'a request body may not exceed 65536 bytes',
        'conventions': {
            'cell_A': T.CELL.tolist(), 'volume_A3': T.VOL, 'reciprocal_basis_B_invA': T.B.round(8).tolist(),
            'k_frac': 'coordinates in the reciprocal basis of the conventional cell (k_cart = k_frac @ B)',
            'energies': 'eV relative to the Fermi level in the header of the model used (Hartree in the header)',
            'velocity': 'dE/dk in eV Angstrom (divide by hbar for the group velocity)',
            'berry_curvature': 'Omega = curl_k i<u|grad_k u>, Angstrom^2',
            'moments': 'mu_B; spin -mu_B <sigma>; Berry-phase orbital (|e|/2hbar) Im <grad u| x (H - E)|grad u>; intra-atomic '
                       'orbital -mu_B <L> (L about the own atom of each d or p function)',
            'tight_binding': 'H(k) with Wannier centres at the atomic positions; band pairs within 5 meV treated as degenerate '
                             '(their mutual terms omitted)',
            'fermi_surface': 'tensors_per_cell[q][i][j] = (1/N_k) sum_nk (-df/dE) X_ij, X = m_i dE/dk_j (mu_B Angstrom) or '
                             '(dE/dk_i) Omega_j (Angstrom^3); sums over the bands of the model per cell of the model; '
                             'electrons_per_cell counts one state per band and k-point'}}
    (AGENT/'catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    for p in (LAB/'public.json', TASK/'tests/public.json'):
        d = json.loads(p.read_text()); d['solver_budget'] = BUDGET
        p.write_text(json.dumps(d, indent=1) + '\n')
    print(json.dumps({'targets': targets, 'robust': robust, 'edges': edges, 'bcd': bcd,
                      'ablations': {k: v for k, v in abl.items()}}, indent=1, default=float))


if __name__ == '__main__':
    main()
