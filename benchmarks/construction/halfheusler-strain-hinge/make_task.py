"""Generate the staged, lab, private and ablation data of halfheusler-strain-hinge.

    bash benchmarks/construction/halfheusler-strain-hinge/run_in_image.sh make_task.py [--reuse]

Runs in the runtime image with the archive (fetch.py DEST) mounted at /data.  It
  1. converts the four FPLO +hamdata files to environment/lab/data/<model>.npz (hamdata.py; values unchanged) and copies
     tb.py into the lab;
  2. stages the archive's code file, the header block of each +hamdata file and the record description under
     environment/agent/archive/, and writes environment/agent/catalog.json;
  3. computes every reference value with the lab's own code path (task_lab.run on the lab's npz files, or tb.Model for
     dense bulk meshes beyond the lab's limits), asserts the archive facts the task rests on, and writes
     private/hidden.json, private/hidden_values.json, solution/ablation_values.json, the budget in both public.json
     files, and refs.json + calibration.json here.
--reuse skips step 3's calculations and rebuilds the outputs from refs.json.
"""
import json
import shutil
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
import hamdata
import tb

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/halfheusler-strain-hinge'
LAB, AGENT = TASK/'environment/lab', TASK/'environment/agent'
SRC = Path('/data/x')
FILES = {'LiSbZn_orig': '+hamdata_lisbzn_eq_orig_setup.txt', 'LiSbZn_strained_setup': '+hamdata_lisbzn_eq_strained_setup.txt',
         'LiBiZn_orig': '+hamdata_libizn_eq_orig_setup.txt', 'LiBiZn_strained_setup': '+hamdata_libizn_eq_strained_setup.txt'}
BUDGET = 30
I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
# The orig prism (axis a1, sides a2, a3 of the orig file) expressed in the strained-setup file's lattice basis:
# strained-setup frame = orig frame rotated by -45 deg about z; a1' = (a/2)(0,1,1), a2' = (a/2)(0,-1,1), a3' = (a/2)(1,0,-1).
STR_BASIS = [[1, 0, 0], [1, 1, 1], [1, 0, 1]]
STR_BASIS_C2 = [[0, 1, 0], [0, 0, -1], [-1, 0, -1]]   # the same prism rotated by C2z (frame rotation 135 deg)
THRESHOLDS = {'strain_percent_abs': 0.25,
              'gap_abs_meV': {'LiSbZn_orig': 3.5, 'LiSbZn_strained_setup': 3.5, 'LiBiZn_orig': 15.0},
              'fermi_level_abs_meV': 12.0, 'hinge_energy_abs_meV': 6.0, 'band_bottom_abs_meV': 6.0,
              'band_top_abs_meV': 15.0, 'zn_weight_abs': 0.08, 'removal_potential_abs_eV': 0.035}


def convert():
    (LAB/'data').mkdir(parents=True, exist_ok=True)
    for name, f in FILES.items():
        d = hamdata.read(SRC/'wannier'/f)
        np.savez_compressed(LAB/'data'/f'{name}.npz', lattice_A=d['lattice_bohr'] * hamdata.BOHR,
                            centers_A=d['centers_bohr'] * hamdata.BOHR, R=d['R'].astype(np.int16), H=d['H'].astype(np.complex128))
    shutil.copy2(HERE/'tb.py', LAB/'tb.py')


def header(path):
    lines = Path(path).read_text().splitlines()
    return '\n'.join(lines[:lines.index('Tij, Hij:')]) + '\n[... hopping blocks "Tij, Hij:" not staged ...]\n'


def stage():
    arch = AGENT/'archive'
    if arch.exists():
        shutil.rmtree(arch)
    (arch/'wannier').mkdir(parents=True)
    (arch/'code').mkdir()
    shutil.copy2(SRC/'code/strained_hH.py', arch/'code/strained_hH.py')
    models = {}
    for name, f in FILES.items():
        (arch/'wannier'/(f + '.header.txt')).write_text(header(SRC/'wannier'/f))
        d = hamdata.read(SRC/'wannier'/f)
        models[name] = {'archive_file': 'wannier/' + f, 'header': 'archive/wannier/' + f + '.header.txt',
                        'compound': name.split('_')[0],
                        'lattice_vectors_bohr': np.round(d['lattice_bohr'], 10).tolist(),
                        'lattice_vectors_A': np.round(d['lattice_bohr'] * hamdata.BOHR, 8).tolist(),
                        'wannier_functions': d['names'],
                        'wannier_centres_bohr': np.round(d['centers_bohr'], 10).tolist(),
                        'wannier_centres_A': np.round(d['centers_bohr'] * hamdata.BOHR, 8).tolist(),
                        'hopping_vectors': int(len(d['R']))}
    (arch/'README.txt').write_text(
        'Zenodo record 16927281, "Hinge electronic structure of strained half-Heuslers" (CC BY 4.0), file strained_hH_zenodo.zip.\n'
        'Record description (verbatim): This repository contains data and code necessary to reproduce the numerical results in:\n'
        'Sanjib Kumar Das, Ion Cosma Fulga, Rakshanda Dhawan, Hem C. Kandpal, Jeroen van den Brink, and Jorge I. Facio,\n'
        '"Hinge electronic structure of strained half-Heuslers" arXiv:2504.20190.\n'
        'The "wannier" folder contains Wannier Hamiltonians as obtained using density-functional theory. The "kwant_ham" folder\n'
        'contains processed data used to construct the tight-binding models. In the "code" folder we include a python file which\n'
        'reproduces some of our main results.\n\n'
        'Staged here: code/strained_hH.py (unchanged) and the header block of each wannier/+hamdata file (everything before the\n'
        'first "Tij, Hij:" hopping block).  The hopping blocks and the kwant_ham pickles (the same hoppings, processed by\n'
        'process_hamdata) are evaluated by the calculation service only.\n')
    return models


def catalog(models):
    req = {
        'bulk_bands': {'cost': 1, 'format': {'kind': 'bulk_bands', 'model': '<model>', 'k_frac or k_cart_invA': '[[x, y, z], ...] (<= 300)'},
                       'returns': 'band energies (eV) and the Zn-4s weight of each band at each k'},
        'band_edges': {'cost': '1 (mesh <= 24) or 2 (mesh 25-48)',
                       'format': {'kind': 'band_edges', 'model': '<model>', 'band': 'n (1-7)', 'mesh': 'm (4-48)'},
                       'returns': 'the 6 highest values of band n, the 6 lowest of band n+1 and the smallest direct gap E_{n+1}-E_n '
                                  'on the Gamma-centred m^3 grid of k_frac (no refinement)'},
        'fermi_level': {'cost': '1 (mesh <= 32) or 2 (mesh 33-64)',
                        'format': {'kind': 'fermi_level', 'model': '<model>', 'electrons_per_cell': 'x (0-8)', 'mesh': 'm (4-64)',
                                   'offset': '0 (Gamma-centred, default) or 0.5'},
                        'returns': 'the energy that holds x electrons per cell on the m^3 grid (each band holds one electron per '
                                   'cell) and the filling of each band'},
        'slab': {'cost': 1,
                 'format': {'kind': 'slab', 'model': '<model>', 'basis': '[[p1], [p2], [n]] integer rows, det +-1',
                            'layers': 'L (2-40)', 'k_phases': '[[phi1, phi2], ...] (<= 30)', 'window_eV': '[emin, emax] (<= 1.2 eV wide)',
                            'region_depth': '1-4 (default 2)'},
                 'returns': 'slab states in the window with their weight near the low and high end of row 3 and their Zn-4s weight'},
        'prism': {'cost': '1 (size <= 12), 2 (13-20), 3 (21-28)',
                  'format': {'kind': 'prism', 'model': '<model>', 'basis': '[[axis], [side1], [side2]] integer rows, det +-1',
                             'size': 'N (4-28)', 'k_phases': '[k, ...] (<= 8 for N <= 12, <= 4 for N 13-20, <= 2 for N 21-28)',
                             'window_eV': '[emin, emax] (<= 1.0 eV wide)', 'region_depth': '1-4 (default 2)',
                             'potential': '[{"region": "c00"|"c0N"|"cN0"|"cNN"|"f1lo"|"f1hi"|"f2lo"|"f2hi", "depth": 1-6, '
                                          '"V_eV": -1.5..1.5, "species": "all"|"Zn"|"X"}, ...] (optional, <= 4)'},
                  'returns': 'prism states in the window (at most 150 per k) with their weight in each corner, facet and bulk '
                             'region and their Zn-4s weight'},
    }
    conv = {
        'energies': 'eV on the energy scale of the archived file; its zero is the zero of the FPLO calculation that produced it',
        'bands': 'numbered from 1 (lowest) among the 8 Wannier bands; with spin-orbit coupling and no inversion centre each band '
                 'is non-degenerate at a generic k and holds one electron per cell',
        'k_bulk': 'k_cart_invA: Cartesian, 1/Angstrom with 2 pi, in the Cartesian frame of that model\'s file; k_frac: coordinates '
                  'in the reciprocal basis of that model\'s lattice vectors (k_cart = k_frac @ B, B = 2 pi inv(A)^T, rows of A = '
                  'lattice vectors)',
        'bulk_hamiltonian': 'H_ij(k) = sum_R H_ij(R) exp(i k.(R + w_j - w_i)), with H_ij(R) the archived matrix element between '
                            'orbital i in the home cell and orbital j displaced by the archived hopping vector T = R + w_j - w_i',
        'supercell_basis': 'slabs and prisms take an integer, unimodular basis matrix whose rows are the new cell vectors in units '
                           'of that model\'s own lattice vectors',
        'site_selection': 'a slab (prism) contains every Zn and pnictogen (X) site whose reduced coordinates along the finite '
                          'basis row(s), measured from an X site, lie in [0, L) ([0, N)); cell t holds the X site at reduced '
                          'coordinates t and the Zn site at t + zn_site_reduced (reported in each result\'s geometry block)',
        'prism_regions': 'cells (t1, t2) along sides 1 and 2, t = 0..N-1; corner cAB = cells within region_depth of the low (0) '
                         'or high (N) end of both sides, e.g. c0N: t1 < depth and t2 >= N - depth; facet f1lo = t1 < depth '
                         '(corners excluded), f1hi, f2lo, f2hi likewise; bulk = the rest.  In a potential entry the facet '
                         'regions include the corners and depth may differ from region_depth',
        'k_prism': 'k_phases are Bloch phases (radians) per period of the periodic row(s) of the basis',
        'onsite_potential': 'V_eV is added to the onsite energy of the selected orbitals (all, the two Zn-4s or the six X-p) of '
                            'every cell in the region',
    }
    data = {'material': 'half-Heusler LiSbZn and LiBiZn (space group 216, Zn on the 4c site), FPLO 18.57 fully relativistic '
                        'DFT, projective Wannier functions (Zn-4s and pnictogen-p spinors, 8 per cell); Li carries no Wannier '
                        'function',
            'source': 'Zenodo record 16927281 (CC BY 4.0), file strained_hH_zenodo.zip; see /app/archive/README.txt',
            'models': models, 'budget': BUDGET, 'requests': req, 'conventions': conv,
            'archive_staged': ['archive/README.txt', 'archive/code/strained_hH.py'] +
                              [m['header'] for m in models.values()]}
    (AGENT/'catalog.json').write_text(json.dumps(data, indent=1) + '\n')
    for p in (LAB/'public.json', TASK/'tests/public.json'):
        pub = json.loads(p.read_text())
        pub['solver_budget'] = BUDGET
        if p.parent.name == 'lab':
            pub['requests'] = {k: v['cost'] for k, v in req.items()}
            pub['catalog'] = '/app/catalog.json'
        p.write_text(json.dumps(pub, indent=1) + '\n')


# ---------------------------------------------------------------- reference calculations
def lab():
    hidden = LAB/'private/hidden.json'
    if not hidden.exists() or hidden.read_text().strip() in ('', '{}'):
        hidden.write_text(json.dumps({'targets': {}, 'thresholds': {}}) + '\n')
    sys.path.insert(0, str(LAB))
    import task_lab
    return task_lab


def refine_edges(m):
    """Global max of band 6 and min of band 7: 40^3 grid then Nelder-Mead from the grid extremum."""
    f = np.arange(40) / 40
    kf = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T
    e = m.bands(kf @ m.B)
    out = {}
    for name, band, sgn in (('vbm', 5, -1), ('cbm', 6, 1)):
        j = int(np.argmax(-sgn * e[:, band]))
        r = minimize(lambda x: sgn * m.bands(x)[0, band], kf[j] @ m.B, method='Nelder-Mead',
                     options={'xatol': 1e-7, 'fatol': 1e-9, 'maxiter': 6000})
        out[name] = {'E_eV': float(sgn * r.fun), 'k_frac': (r.x @ np.linalg.inv(m.B)).tolist(), 'grid_E_eV': float(e[j, band]),
                     'grid_k_frac': kf[j].tolist()}
    return out


def fermi(m, mesh, offset):
    f = (np.arange(mesh) + offset) / mesh
    kf = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T
    flat = np.sort(m.bands(kf @ m.B).ravel())
    n = 6 * mesh ** 3
    return float(0.5 * (flat[n - 1] + flat[n]))


def prism_states(L, model, basis, n, k, window, potential=None, depth=2):
    req = {'kind': 'prism', 'model': model, 'basis': basis, 'size': n, 'k_phases': [k], 'window_eV': window, 'region_depth': depth}
    if potential:
        req['potential'] = potential
    assert L.validate(req) >= 1
    return L.run(req, None)


def hinge_level(L, model, basis, n, k, vbm, cbm, potential=None, corner='cNN'):
    """Lowest state in [vbm - 0.03, cbm] at k whose weight in the 2x2 corner exceeds 0.3."""
    res = prism_states(L, model, basis, n, k, [vbm - 0.03, cbm], potential)['results'][0]
    for s in res['states']:
        if s['weights'][corner] > 0.3:
            return s
    return None


PARTIAL = {}


def compute(L):
    global PARTIAL
    if (HERE/'refs_partial.json').exists():
        PARTIAL = json.loads((HERE/'refs_partial.json').read_text())
    refs = {}
    M = {name: L.model(name) for name in FILES}
    # --- geometry: strain of each model relative to the orig cubic lattice parameter
    strain = {}
    for name in FILES:
        A = M[name].A
        orig = M[name.split('_')[0] + '_orig'].A
        a0 = 2 * abs(orig[0, 1])                       # orig: rows (0, a/2, a/2), ...
        if 'strained' in name:                          # bct rows (-at/2, at/2, c/2), ...
            c, at = 2 * abs(A[0, 2]), 2 * abs(A[0, 0])
            assert np.allclose(np.abs(A), np.abs(A[0])[None], atol=1e-8)
        else:
            c, at = a0, a0 / np.sqrt(2)
        strain[name] = {'c_A': c, 'a_inplane_cubic_equiv_A': at * np.sqrt(2), 'delta_percent': 100 * (a0 - c) / a0,
                        'volume_A3': abs(np.linalg.det(A))}
        assert abs(strain[name]['delta_percent']) < 1e-3 and abs(at * np.sqrt(2) - a0) < 1e-4, (name, strain[name])
    refs['strain'] = strain
    # --- Gamma levels, inversion, band edges
    gamma, edges = {}, {}
    for name in FILES:
        r = L.run({'kind': 'bulk_bands', 'model': name, 'k_frac': [[0, 0, 0]]}, None)
        gamma[name] = {'E_eV': r['energies_eV'][0], 'zn_s_weight': r['zn_s_weight'][0]}
        edges[name] = refine_edges(M[name])
    refs['gamma'], refs['edges'] = gamma, edges
    for name in FILES:
        e, zn = np.array(gamma[name]['E_eV']), np.array(gamma[name]['zn_s_weight'])
        assert np.ptp(e[2:6] if 'Sb' in name else e[4:8]) < 2e-4, (name, e)   # fourfold Gamma8, intact: no strain
    # LiSbZn: Gamma6 (Zn-s doublet) above the Gamma8 quartet; LiBiZn: below it
    assert gamma['LiSbZn_orig']['zn_s_weight'][6] > 0.99 and gamma['LiBiZn_orig']['zn_s_weight'][2] > 0.99
    # --- LiBiZn Fermi level (6 electrons) on dense meshes
    m = M['LiBiZn_orig']
    refs['libizn_fermi'] = {f'{n}_{o}': fermi(m, n, o) for n, o in ((8, 0.0), (12, 0.0), (16, 0.0), (24, 0.0), (32, 0.0),
                                                                    (48, 0.0), (64, 0.0), (64, 0.5), (96, 0.5))}
    lab_ef = L.run({'kind': 'fermi_level', 'model': 'LiBiZn_orig', 'electrons_per_cell': 6, 'mesh': 64, 'offset': 0.5}, None)
    assert abs(lab_ef['fermi_level_eV'] - refs['libizn_fermi']['64_0.5']) < 1e-6
    # --- hinge state, LiSbZn orig and strained setup
    hinge = PARTIAL.get('hinge_sizes', {})
    for name, basis in (('LiSbZn_orig', I3), ('LiSbZn_strained_setup', STR_BASIS)):
        if name in hinge:
            continue
        vbm, cbm = edges[name]['vbm']['E_eV'], edges[name]['cbm']['E_eV']
        rows = {}
        for n in (6, 8, 10, 12, 16, 20, 24, 28):
            s = hinge_level(L, name, basis, n, 0.0, vbm, cbm)
            assert s is not None and max(s['weights'], key=s['weights'].get) == 'cNN', (name, n, s)
            rows[n] = {'E_eV': s['E_eV'], 'zn_s_weight': s['zn_s_weight'], 'cNN': s['weights']['cNN']}
            print(name, 'N', n, rows[n], flush=True)
        hinge[name] = {'by_size_k0': {str(k): v for k, v in rows.items()}, 'vbm_eV': vbm, 'cbm_eV': cbm}
        PARTIAL['hinge_sizes'] = hinge
        (HERE/'refs_partial.json').write_text(json.dumps(PARTIAL, indent=1, default=float) + '\n')
    for name in hinge:
        hinge[name].update(vbm_eV=edges[name]['vbm']['E_eV'], cbm_eV=edges[name]['cbm']['E_eV'])
    refs['hinge'] = hinge
    # geometry of the orig prism: the hinge carrying the state is the high-high corner; which species is outermost there?
    geo = L.geometry(M['LiSbZn_orig'], I3)
    zr = np.array(geo['zn_site_reduced'])
    assert np.allclose(zr, 0.25)       # Zn at +1/4 along both sides: X outermost at t = 0 facets, Zn outermost at t = N-1
    refs['hinge_geometry'] = geo
    # strained-setup prism in its own basis reproduces the orig cell content
    assert np.allclose(L.geometry(M['LiSbZn_strained_setup'], STR_BASIS)['zn_site_reduced'], 0.25)
    # hinge band over the 1D BZ (orig): both spin-split branches; coarse scan at N = 12 and N = 20, refined top
    vbm, cbm = edges['LiSbZn_orig']['vbm']['E_eV'], edges['LiSbZn_orig']['cbm']['E_eV']
    m = M['LiSbZn_orig']

    def corner_states(n, k, window, onsite=None, corner='cNN', wmin=0.5, model=m, basis=I3):
        """The lab's prism evaluator (tb.finite_matrix/finite_states, task_lab.region_masks) without the 40-state cap."""
        Mx, cells = model.finite_matrix(basis, [n, n], [float(k)], onsite=onsite)
        e, v, _ = model.finite_states(Mx, window[0], window[1], nmax=2000)
        w = (np.abs(v) ** 2).reshape(len(cells), tb.NW, -1)
        mask = L.region_masks(cells, n, 2)[corner]
        cw = w[mask].sum((0, 1))
        return [float(x) for x, c in zip(e, cw) if c > wmin]

    def stage(key, fn):
        if key not in PARTIAL:
            PARTIAL[key] = fn()
            (HERE/'refs_partial.json').write_text(json.dumps(PARTIAL, indent=1, default=float) + '\n')
        print(key, json.dumps(PARTIAL[key], default=float)[:1500], flush=True)
        return PARTIAL[key]

    win = [vbm + 0.02, 0.95]
    scan = stage('scan', lambda: {str(n): [[float(k)] + corner_states(n, k, win)[:2] for k in np.linspace(0, np.pi, npts)]
                                  for n, npts in ((12, 25), (20, 13))})

    def top_of(n):
        r = minimize_scalar(lambda k: -corner_states(n, k, win)[1], bounds=(0.85, 1.25), method='bounded', options={'xatol': 1e-3})
        return {'k': float(r.x), 'E_eV': float(-r.fun)}
    tops = stage('top', lambda: {str(n): top_of(n) for n in (12, 20)})
    bottom = stage('bottom_k_pi', lambda: {str(n): corner_states(n, np.pi, [vbm + 0.02, cbm])[0] for n in (12, 20, 28)})
    grid5 = stage('grid5', lambda: [[float(k)] + corner_states(12, k, win)[:2] for k in np.linspace(0, np.pi, 5)])
    grid9 = stage('grid9', lambda: [[float(k)] + corner_states(12, k, win)[:2] for k in np.linspace(0, np.pi, 9)])
    refs['hinge_band'] = {'scan': scan, 'top': tops, 'bottom_k_pi': bottom, 'grid5': grid5, 'grid9': grid9}

    # removal potential: V on all orbitals of the 2x2 cells at the cNN corner; the k = 0 level reaches the bulk VBM
    def level(V, n, corner='cNN', wmin=0.3):
        t1, t2 = [x.ravel() for x in np.meshgrid(np.arange(n), np.arange(n), indexing='ij')]
        sel = {'cNN': (t1 >= n - 2) & (t2 >= n - 2), 'c00': (t1 < 2) & (t2 < 2)}[corner]
        ons = np.zeros((n * n, tb.NW))
        ons[sel] = V
        return corner_states(n, 0.0, [vbm - 0.08, cbm], ons, 'cNN', wmin)[0]

    def removal(n):
        pts = {V: level(V, n) for V in (-0.2, -0.3, -0.35, -0.4, -0.45)}
        return {'levels': {str(k): v for k, v in pts.items()},
                'brentq': brentq(lambda V: level(V, n) - vbm, -0.45, -0.3, xtol=1e-4),
                'linear_0.2_0.4': -0.2 - 0.2 * (pts[-0.2] - vbm) / (pts[-0.2] - pts[-0.4]),
                'quadratic_0.30_0.38_0.46': None}
    rem = stage('removal', lambda: {str(n): removal(n) for n in (12, 16, 20)})
    # the oracle's protocol: three potentials, quadratic interpolation
    def quad(n):
        vs = [-0.30, -0.38, -0.46]
        lv = [level(V, n) - vbm for V in vs]
        roots = [r.real for r in np.roots(np.polyfit(vs, lv, 2)) if abs(r.imag) < 1e-9 and -0.6 < r.real < -0.2]
        return roots[0]
    rem_q = stage('removal_quadratic', lambda: {str(n): quad(n) for n in (12, 16)})
    # the wrong corner (X/X, c00): the cNN level does not move
    wrong = stage('removal_wrong_corner', lambda: {str(V): level(V, 16, 'c00', 0.15) for V in (-1.5, -1.0, 0.0, 1.0, 1.5)})
    refs['removal'] = {**rem, 'quadratic': rem_q, 'wrong_corner_levels_eV': wrong}

    # --- ablation input: the orig integer basis applied to the strained-setup file (a (100)/(111) prism)
    def identity_basis():
        name = 'LiSbZn_strained_setup'
        vbm_s, cbm_s = edges[name]['vbm']['E_eV'], edges[name]['cbm']['E_eV']
        res = prism_states(L, name, I3, 16, 0.0, [vbm_s, cbm_s])['results'][0]
        best = max(res['states'], key=lambda s: max(s['weights'][c] for c in ('c00', 'c0N', 'cN0', 'cNN')))
        return {'n_in_gap_states': len(res['states']), 'max_corner_state': best,
                'facet_dominated': sum(max(s['weights'], key=s['weights'].get).startswith('f') for s in res['states'])}
    refs['str_identity_basis'] = stage('str_identity_basis', identity_basis)
    # the C2z-equivalent basis of the strained-setup file gives the same level
    refs['str_c2z_basis_N20'] = stage('str_c2z_basis_N20', lambda: hinge_level(
        L, 'LiSbZn_strained_setup', STR_BASIS_C2, 20, 0.0,
        edges['LiSbZn_strained_setup']['vbm']['E_eV'], edges['LiSbZn_strained_setup']['cbm']['E_eV']))
    return refs


def removal_to_zero(levels):
    """V at which the k = 0 hinge level crosses the FPLO zero (mistaken for the VBM / E_F), linear between tabulated V."""
    pts = sorted((float(k), v) for k, v in levels.items())
    for (v1, e1), (v2, e2) in zip(pts, pts[1:]):
        if e1 * e2 <= 0:
            return v1 + (v2 - v1) * (0 - e1) / (e2 - e1)
    raise ValueError('no crossing')


def write(refs):
    E = refs['edges']
    gap = {k: 1000 * (E[k]['cbm']['E_eV'] - E[k]['vbm']['E_eV']) for k in ('LiSbZn_orig', 'LiSbZn_strained_setup', 'LiBiZn_orig')}
    g8 = refs['gamma']['LiBiZn_orig']['E_eV'][4]
    ef = 1000 * (refs['libizn_fermi']['96_0.5'] - g8)
    H = refs['hinge']
    for k in H:   # geometric (exponential) extrapolation in N from N = 20, 24, 28
        e20, e24, e28 = (H[k]['by_size_k0'][str(n)]['E_eV'] for n in (20, 24, 28))
        r = (e28 - e24) / (e24 - e20)
        assert 0 < r < 1
        H[k]['extrapolated_E_eV'] = e28 + (e28 - e24) * r / (1 - r)
    ek = {k: 1000 * (H[k]['extrapolated_E_eV'] - H[k]['vbm_eV']) for k in H}
    vbm = H['LiSbZn_orig']['vbm_eV']
    hb = refs['hinge_band']
    window = [1000 * (hb['bottom_k_pi']['28'] - vbm), 1000 * (hb['top']['20']['E_eV'] - vbm)]
    zn = H['LiSbZn_orig']['by_size_k0']['28']['zn_s_weight']
    rp = refs['removal']['20']['brentq']
    targets = {'strain_percent': {k: 0.0 for k in FILES},
               'band_inverted': {'LiSbZn': False, 'LiBiZn': True},
               'fundamental_gap_meV': gap, 'libizn_fermi_level_meV': ef,
               'hinge_termination': 'Zn/Zn', 'hinge_energy_k0_meV': ek, 'hinge_band_window_meV': window,
               'hinge_zn_s_weight_k0': zn, 'hinge_removal_potential_eV': rp}
    (LAB/'private/hidden.json').write_text(json.dumps({'targets': targets, 'thresholds': THRESHOLDS}, indent=1) + '\n')
    values = list(gap.values()) + [ef] + list(ek.values()) + window + [zn, rp]
    (LAB/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    # ---- ablations (realistic wrong paths, from the data)
    G = refs['gamma']
    by = lambda k, n: H[k]['by_size_k0'][str(n)]
    sid = refs['str_identity_basis']['max_corner_state']
    g5 = hb['grid5']
    abl = {
        'strain_from_file_names': {k: (2.0 if 'strained' in k else 0.0) for k in FILES},
        'strain_compressive_from_file_names': {k: (-2.0 if 'strained' in k else 0.0) for k in FILES},
        'gap_at_gamma': {'LiSbZn_orig': 1000 * (G['LiSbZn_orig']['E_eV'][6] - G['LiSbZn_orig']['E_eV'][5]),
                         'LiSbZn_strained_setup': 1000 * (G['LiSbZn_strained_setup']['E_eV'][6] - G['LiSbZn_strained_setup']['E_eV'][5]),
                         'LiBiZn_orig': 0.0},
        'fermi_at_gamma8': 0.0,
        'fermi_coarse_mesh_8': 1000 * (refs['libizn_fermi']['8_0.0'] - g8),
        'fermi_mesh_12': 1000 * (refs['libizn_fermi']['12_0.0'] - g8),
        'hinge_energy_model_zero': {k: 1000 * H[k]['extrapolated_E_eV'] for k in H},
        'hinge_energy_small_prism_N10': {k: 1000 * (by(k, 10)['E_eV'] - H[k]['vbm_eV']) for k in H},
        'hinge_energy_N12': {k: 1000 * (by(k, 12)['E_eV'] - H[k]['vbm_eV']) for k in H},
        'paper_termination': 'X/X', 'paper_zn_s_weight': 0.05,
        'window_grid5': [1000 * (min(r[1] for r in g5 if len(r) > 1) - vbm), 1000 * (max(r[-1] for r in g5 if len(r) > 1) - vbm)],
        'window_clipped_to_bulk_gap': [window[0], gap['LiSbZn_orig']],
        'window_bottom_at_k0': [ek['LiSbZn_orig'], window[1]],
        'window_lower_branch_only': [window[0], 1000 * (max(r[1] for r in hb['scan']['20'] if len(r) > 1) - vbm)],
        'str_hinge_identity_basis': 1000 * (sid['E_eV'] - H['LiSbZn_strained_setup']['vbm_eV']),
        'removal_potential_wrong_corner': None,
        'removal_to_model_zero': removal_to_zero(refs['removal']['20']['levels']),
    }
    rem16 = refs['removal']['16']
    rob = {
        'hinge_energy_N28_direct': {k: 1000 * (by(k, 28)['E_eV'] - H[k]['vbm_eV']) for k in H},
        'hinge_energy_N20_direct': {k: 1000 * (by(k, 20)['E_eV'] - H[k]['vbm_eV']) for k in H},
        'gap_grid_L_point': {k: 1000 * (E[k]['cbm']['grid_E_eV'] - E[k]['vbm']['grid_E_eV']) for k in gap},
        'fermi_mesh_48': 1000 * (refs['libizn_fermi']['48_0.0'] - g8),
        'fermi_mesh_64_lab_max': 1000 * (refs['libizn_fermi']['64_0.5'] - g8),
        'window_N12_refined': [1000 * (hb['bottom_k_pi']['12'] - vbm), 1000 * (hb['top']['12']['E_eV'] - vbm)],
        'window_N12_grid25': [1000 * (min(r[1] for r in hb['scan']['12'] if len(r) > 1) - vbm), 1000 * (max(r[-1] for r in hb['scan']['12'] if len(r) > 1) - vbm)],
        'zn_weight_N16': by('LiSbZn_orig', 16)['zn_s_weight'],
        'removal_N16_linear': rem16['linear_0.2_0.4'],
        'removal_N16_brentq': rem16['brentq'],
    }
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    cal = {'targets': targets, 'thresholds': THRESHOLDS, 'ablations': abl, 'robustness': rob}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1) + '\n')
    print(json.dumps(cal, indent=1))


def main():
    convert()
    models = stage()
    catalog(models)
    if '--reuse' in sys.argv:
        refs = json.loads((HERE/'refs.json').read_text())
    else:
        refs = compute(lab())
        (HERE/'refs.json').write_text(json.dumps(refs, indent=1, default=float) + '\n')
        refs = json.loads((HERE/'refs.json').read_text())
    write(refs)


if __name__ == '__main__':
    main()
