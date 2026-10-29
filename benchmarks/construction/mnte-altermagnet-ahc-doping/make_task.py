"""Build the staged, lab, private and ablation files of mnte-altermagnet-ahc-doping.

    bash benchmarks/construction/mnte-altermagnet-ahc-doping/run_in_image.sh make_task.py

Needs fetch.py (archive members under /data), make_model.py (/data/model/archived.npz) and make_cache.py
(/data/cache/mesh_N.npz for N in MESHES). Runs in the runtime image (numpy 1.21) in about a minute and writes:
  environment/agent/archive/...       staged archive files (unchanged; WT.out without its k-point progress lines)
  environment/agent/catalog.json      staged files, requests, costs, conventions
  environment/lab/data/               archived_hr.npz (H(R) as integer micro-eV, exact) and the mesh caches
  environment/lab/mnte_tb.py          the evaluator (copied unchanged)
  environment/lab/private/hidden.json, hidden_values.json
  solution/ablation_values.json, and calibration.json here.
Graded sigma_xy values are the authors' archived WannierTools output (4_wtools/sigma_ahc.txt, 801 x 801 x 401 mesh,
WannSymm-symmetrized Hamiltonian), read at the chemical potentials where the lab evaluator's 400-mesh electron count
gives the requested hole densities. Hole densities come from the lab evaluator.
"""
import json
import shutil
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import mnte_tb as T

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mnte-altermagnet-ahc-doping'
SRC = Path('/data')
MESHES = [48, 72, 100, 150, 200, 300, 400]
E_FERMI_WT = 6.13012952            # wt.in: the energy zero of sigma_ahc.txt
E_FERMI_OUTCAR = 6.2098            # 1_vasp/OUTCAR, E-fermi of the VASP run (in the gap)
NOCC = 22
VOLCM3 = T.VOL * 1e-24
KT = 0.001
P_TARGETS = {'2.0e21': 2.0e21, '6.0e21': 6.0e21}
N_ELECTRON = 1e20
WINDOW = (1.00, 0.25)              # sign changes counted for mu from VBM - 1.00 to VBM - 0.25 eV
STEP = 0.005                       # ... on a 5 meV grid (151 chemical potentials)
THRESHOLDS = {'sigma_xy_rel': 0.15, 'sign_change_density_rel': 0.08, 'hole_density_rel': 0.05, 'forbidden_abs_S_cm': 0.5}
STAGED = ['1_vasp/INCAR', '1_vasp/POSCAR', '1_vasp/KPOINTS', '1_vasp/IBZKPT', '1_vasp/OSZICAR', '1_vasp/OUTCAR',
          '1_vasp/XDATCAR', '1_vasp/PCDAT', '1_vasp/log', '1_vasp/wannier90.win', '1_vasp/wannier90.wout',
          '2_wannier90/wannier90.win', '2_wannier90/wannier90.wout', '2_wannier90/wannier90.eig',
          '2_wannier90/wannier90.nnkp', '2_wannier90/wannier90_centres.xyz', '2_wannier90/wannier90_band.dat',
          '2_wannier90/wannier90_band.kpt', '2_wannier90/wannier90_band.gnu', '2_wannier90/wannier90_band.labelinfo.dat',
          '3_wannsymm/wannsymm.in', '3_wannsymm/wannsymm.out', '3_wannsymm/symmetries.dat', '3_wannsymm/POSCAR',
          '3_wannsymm/bands_ori.dat', '3_wannsymm/bands_symmed.dat',
          '4_wtools/wt.in', '4_wtools/sigma_ahc.gnu', '4_wtools/WT.out']
FIG6 = ['fig6_script/0_0.450.dat', 'fig6_script/90_0.450.dat', 'fig6_script/Berrycurvature.gnu']


def stage():
    dest = TASK/'environment/agent/archive'
    if dest.exists():
        shutil.rmtree(dest)
    for name in STAGED:
        target = dest/name
        target.parent.mkdir(parents=True, exist_ok=True)
        if name == '4_wtools/WT.out':
            lines = (SRC/'dft_data'/name).read_text().splitlines(keepends=True)
            target.write_text(''.join(l for l in lines if not l.startswith('ik/knv3')))
        else:
            shutil.copy2(SRC/'dft_data'/name, target)
    for name in FIG6:
        (dest/name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SRC/name, dest/name)
    return ['/'.join(p.relative_to(dest).parts) for p in sorted(dest.rglob('*')) if p.is_file()]


def lab_data():
    data = TASK/'environment/lab/data'
    data.mkdir(parents=True, exist_ok=True)
    z = np.load(SRC/'model/archived.npz')
    re, im = np.round(z['H'].real * 1e6), np.round(z['H'].imag * 1e6)
    assert np.abs(re - z['H'].real * 1e6).max() < 1e-6 and np.abs(im - z['H'].imag * 1e6).max() < 1e-6
    np.savez_compressed(data/'archived_hr.npz', R=z['R'].astype(np.int8), deg=z['deg'].astype(np.int16),
                        H_re_ueV=re.astype(np.int32), H_im_ueV=im.astype(np.int32))
    for n in MESHES:
        shutil.copy2(SRC/f'cache/mesh_{n}.npz', data/f'mesh_{n}.npz')
    shutil.copy2(HERE/'mnte_tb.py', TASK/'environment/lab/mnte_tb.py')
    m1, m2 = T.Model(SRC/'model/archived.npz'), T.Model(data/'archived_hr.npz')
    k = np.random.RandomState(0).rand(200, 3) @ T.B
    assert np.abs(m1.bands(k) - m2.bands(k)).max() < 1e-9


def find_vbm(m):
    n = 24
    f = np.arange(n) / n
    kf = np.array(np.meshgrid(f, f, np.arange(12) / 12, indexing='ij')).reshape(3, -1).T
    e = m.bands(kf @ T.B)[:, NOCC - 1]
    best = None
    for i in np.argsort(-e)[:40]:
        r = minimize(lambda x: -m.bands((x @ T.B)[None])[0, NOCC - 1], kf[i], method='Nelder-Mead',
                     options={'xatol': 1e-9, 'fatol': 1e-11, 'maxiter': 6000})
        if best is None or r.fun < best.fun:
            best = r
    return -best.fun, best.x - np.round(best.x)


def path_vbm(m):
    pts = np.array([[1/3, 1/3, 0], [0, 0, 0], [0, 0, .5], [1/3, 1/3, 0], [1/3, 1/3, .5], [.5, 0, 0], [.5, 0, .5]])
    k = np.concatenate([pts[a] + np.linspace(0, 1, 101)[:, None] * (pts[a + 1] - pts[a]) for a in range(6)])
    return float(m.bands(k @ T.B)[:, NOCC - 1].max())


def dft_grid():
    win = (SRC/'dft_data/2_wannier90/wannier90.win').read_text()
    kd = np.array([[float(x) for x in l.split()[:3]] for l in win.split('begin kpoints')[1].split('end kpoints')[0].splitlines() if l.strip()])
    eig = np.loadtxt(SRC/'dft_data/2_wannier90/wannier90.eig')[:, 2].reshape(len(kd), -1)
    return kd, eig


def sign_grid(vbm):
    return vbm - np.round(np.arange(WINDOW[0], WINDOW[1] - 1e-9, -STEP), 6)


def crossings(mu, s):
    sg = np.sign(s)
    j = np.where(sg[1:] * sg[:-1] < 0)[0]
    return [float(mu[i] - s[i] * (mu[i + 1] - mu[i]) / (s[i + 1] - s[i])) for i in j]


def main():
    files = stage()
    lab_data()
    m = T.Model(SRC/'model/archived.npz')
    vbm, kvbm = find_vbm(m)
    vbm_path = path_vbm(m)
    kd, eig = dft_grid()
    w = m.bands(kd @ T.B)
    vbm_dftgrid = float(eig[:, NOCC - 1].max())
    caches = {n: dict(np.load(SRC/f'cache/mesh_{n}.npz')) for n in MESHES}
    holes = lambda n, mu, kT=KT: (NOCC - T.electrons_from_cache(caches[n], mu, kT)) / VOLCM3
    arch = np.loadtxt(SRC/'dft_data/4_wtools/sigma_ahc.txt')
    amu, asx = E_FERMI_WT + arch[:, 0], arch[:, 1]
    sig_arch = lambda mu: float(np.interp(mu, amu, asx))
    # mu(p) on the densest mesh
    mug = np.linspace(vbm - 1.05, vbm + 0.01, 10601)
    def mu_of_p(p, n=400, kT=KT):
        h = holes(n, mug, kT)
        return float(np.interp(p, h[::-1], mug[::-1]))
    mu_p = {k: mu_of_p(p) for k, p in P_TARGETS.items()}
    sigma_targets = {k: sig_arch(v) for k, v in mu_p.items()}
    # sign changes of the archived curve in the window
    g = sign_grid(vbm)
    xa = crossings(g, np.interp(g, amu, asx))
    p0 = float(holes(400, np.array([min(xa)]))[0])
    p030 = float(holes(400, np.array([vbm - 0.30]))[0])
    targets = {'sigma_xy_S_per_cm': sigma_targets, 'sigma_xy_sign_changes': len(xa), 'sign_change_hole_density_cm3': p0,
               'hole_density_0p30_below_vbm_cm3': p030, 'forbidden_components_max_abs_S_per_cm': {'yz': 0.0, 'zx': 0.0},
               'allowed_hall_components': {'neel_y': ['xy'], 'neel_x': []}, 'electron_doped_sigma_xy_S_per_cm': None}

    # --- model-side values (what the lab gives) on every mesh
    mu_scan = np.linspace(vbm - 1.0, vbm, 2001)
    per_mesh = {}
    for n in MESHES:
        s = T.sigma_from_cache(caches[n], mu_scan, KT)
        xs = crossings(g, T.sigma_from_cache(caches[n], g, KT)[2])
        mp = {k: mu_of_p(p, n) for k, p in P_TARGETS.items()}
        per_mesh[n] = {'sigma_xy_at_p': {k: float(T.sigma_from_cache(caches[n], [v], KT)[2][0]) for k, v in mp.items()},
                       'sign_changes': len(xs), 'crossings_eV_below_vbm': [round(vbm - x, 4) for x in xs],
                       'sign_change_density': float(holes(n, np.array([min(xs)]))[0]) if xs else None,
                       'hole_density_0p30': float(holes(n, np.array([vbm - 0.30]))[0]),
                       'max_abs_yz': float(np.abs(s[0]).max()), 'max_abs_zx': float(np.abs(s[1]).max())}
    # smearing variants on the densest mesh
    smear = {}
    for kT in (0.0, 0.0005, 0.0086, 0.025):
        xs = crossings(g, T.sigma_from_cache(caches[400], g, kT)[2])
        smear[kT] = {'sigma_xy_at_p': {k: float(T.sigma_from_cache(caches[400], [mu_of_p(p, 400, kT)], kT)[2][0]) for k, p in P_TARGETS.items()},
                     'sign_changes': len(xs)}
    # --- ablations (realistic wrong paths, from the data)
    A, R = {}, {}
    for coarse in (100, 150, 200):
        A[f'mesh_{coarse}'] = {'sigma_xy_S_per_cm': per_mesh[coarse]['sigma_xy_at_p'],
                               'sigma_xy_sign_changes': per_mesh[coarse]['sign_changes'],
                               'sign_change_hole_density_cm3': per_mesh[coarse]['sign_change_density'],
                               'hole_density_0p30_below_vbm_cm3': per_mesh[coarse]['hole_density_0p30']}
    A['vbm_from_band_path'] = float(holes(400, np.array([vbm_path - 0.30]))[0])
    A['vbm_from_dft_grid'] = float(holes(400, np.array([vbm_dftgrid - 0.30]))[0])
    A['vbm_from_outcar_efermi'] = float(holes(400, np.array([E_FERMI_OUTCAR - 0.30]))[0])
    # spin degeneracy counted twice: reported densities doubled; sigma read where the true density is half
    A['spin_degeneracy_twice'] = {'sigma_xy_S_per_cm': {k: sig_arch(mu_of_p(p / 2)) for k, p in P_TARGETS.items()},
                                  'sign_change_hole_density_cm3': 2 * p0, 'hole_density_0p30_below_vbm_cm3': 2 * p030}
    scan500 = np.linspace(vbm - 1.0, vbm, 500)
    raw = T.sigma_from_cache(caches[400], scan500, KT)
    A['raw_forbidden_components'] = {'yz': float(np.abs(raw[0]).max()), 'zx': float(np.abs(raw[1]).max())}
    raw_deep = T.sigma_from_cache(caches[400], scan500[scan500 <= vbm - 0.25], KT)
    A['raw_forbidden_components_deep'] = {'yz': float(np.abs(raw_deep[0]).max()), 'zx': float(np.abs(raw_deep[1]).max())}
    A['symmetries_dat_literal'] = {'neel_y': [], 'neel_x': []}
    A['inplane_isotropic'] = {'neel_y': ['xy'], 'neel_x': ['xy']}
    ne = NOCC + N_ELECTRON * VOLCM3
    mug_c = np.linspace(7.80, 8.30, 5001)
    ec = T.electrons_from_cache(caches[400], mug_c, KT)
    mu_n = float(np.interp(ne, ec, mug_c))
    A['electron_doped_from_model'] = float(T.sigma_from_cache(caches[400], [mu_n], KT)[2][0])
    A['room_temperature_smearing'] = {'sigma_xy_S_per_cm': smear[0.025]['sigma_xy_at_p'], 'sigma_xy_sign_changes': smear[0.025]['sign_changes']}
    A['wt_out_100K_smearing'] = {'sigma_xy_S_per_cm': smear[0.0086]['sigma_xy_at_p'], 'sigma_xy_sign_changes': smear[0.0086]['sign_changes']}
    R['mesh_300'] = {'sigma_xy_S_per_cm': per_mesh[300]['sigma_xy_at_p'], 'sigma_xy_sign_changes': per_mesh[300]['sign_changes'],
                     'sign_change_hole_density_cm3': per_mesh[300]['sign_change_density'],
                     'hole_density_0p30_below_vbm_cm3': per_mesh[300]['hole_density_0p30']}
    R['mesh_400_model_values'] = {'sigma_xy_S_per_cm': per_mesh[400]['sigma_xy_at_p'],
                                  'sign_change_hole_density_cm3': per_mesh[400]['sign_change_density']}
    R['wt_in_efermi_as_vbm'] = float(holes(400, np.array([E_FERMI_WT - 0.30]))[0])
    R['step_occupation'] = {'sigma_xy_S_per_cm': smear[0.0]['sigma_xy_at_p'], 'sigma_xy_sign_changes': smear[0.0]['sign_changes']}
    # --- write
    hidden = {'targets': targets, 'thresholds': THRESHOLDS,
              'reference': {'vbm_eV': vbm, 'vbm_k_frac': kvbm.round(5).tolist(), 'mu_of_p_eV': mu_p,
                            'archive_crossings_eV': xa, 'archive_energy_zero_eV': E_FERMI_WT}}
    priv = TASK/'environment/lab/private'
    priv.mkdir(parents=True, exist_ok=True)
    (priv/'hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(sigma_targets.values()) + [p0, p030]
    (priv/'hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': R}, indent=1) + '\n')
    cal = {'vbm_eV': vbm, 'vbm_k_frac': kvbm.tolist(), 'vbm_path_eV': vbm_path, 'vbm_dft_grid_eV': vbm_dftgrid,
           'wt_in_E_FERMI': E_FERMI_WT, 'outcar_E_fermi': E_FERMI_OUTCAR,
           'wannier_vs_dft_on_grid_eV': {'valence_max_abs': float(np.abs(w[:, :NOCC] - eig[:, :NOCC]).max()),
                                         'conduction_23_24_max_abs': float(np.abs(w[:, NOCC:NOCC + 2] - eig[:, NOCC:NOCC + 2]).max()),
                                         'dft_cbm': float(eig[:, NOCC].min()), 'wannier_cbm_on_grid': float(w[:, NOCC].min())},
           'targets': targets, 'mu_of_p_eV_below_vbm': {k: vbm - v for k, v in mu_p.items()},
           'archive_sigma_at_p': sigma_targets, 'per_mesh': {str(n): v for n, v in per_mesh.items()},
           'smearing_400': {str(k): v for k, v in smear.items()}, 'ablations': A, 'robustness': R,
           'electron_doping_mu_eV': mu_n, 'staged_files': files}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1, default=float) + '\n')
    write_catalog(files, vbm)
    print(json.dumps({'targets': targets, 'per_mesh': per_mesh, 'ablations': A, 'robustness': R}, indent=1, default=float))


def write_catalog(files, vbm):
    cat = {
        'material': 'hexagonal MnTe (NiAs type, P6_3/mmc), HSE06 + spin-orbit coupling (VASP 5.4.4), in-plane Neel order',
        'source': 'Zenodo record 18863304 (DOI 10.5281/zenodo.18863304, CC BY 4.0): data for K. P. Kluczyk et al., '
                  'Phys. Rev. B 110, 155201 (2024): 1_vasp, 2_wannier90, 3_wannsymm, 4_wtools of dft_data.zip and fig6_script.zip',
        'staged_files': {'root': '/app/archive', 'files': files,
                         'note': 'archive files unchanged, except 4_wtools/WT.out without its k-point progress lines '
                                 '("ik/knv3 ..."); not staged: 2_wannier90/wannier90_hr.dat (evaluated by the service) and '
                                 '4_wtools/sigma_ahc.txt (withheld)'},
        'wannier_model': {'source': '2_wannier90/wannier90_hr.dat', 'num_wann': 32, 'nrpts': 387,
                          'basis': 'Mn d (dz2, dxz, dyz, dx2-y2, dxy) on Mn1, Mn2, then Te p (pz, px, py) on Te1, Te2, spin up '
                                   '(functions 1-16) then spin down (17-32), spin quantised along Cartesian z',
                          'neel': {'+y': 'archived state (MAGMOM 0 4 0 / 0 -4 0)',
                                   '-y': 'time reverse of the archived Hamiltonian, H(R) -> (i s_y) H(R)* (i s_y)^dagger: all moments reversed'}},
        'requests': {
            'bands': {'format': {'kind': 'bands', 'k_frac or k_cart_invA': '[[x, y, z], ...]', 'neel': '+y (default) or -y',
                                 'spin': 'false (default; <= 1000 k-points) or true (<= 500 k-points)'}, 'cost': 1},
            'berry_curvature': {'format': {'kind': 'berry_curvature', 'k_frac or k_cart_invA': '[[x, y, z], ...] (<= 500)',
                                           'neel': '+y or -y'}, 'cost': 1},
            'ahc': {'format': {'kind': 'ahc', 'mesh': 'N, one of 48, 72, 100, 150, 200, 300, 400 (N x N x N/2 mesh)',
                               'mu_eV': '[...] (1-500 chemical potentials, 5.05-8.55 eV)', 'kT_eV': '0-0.025 (0 = step occupation)',
                               'neel': '+y or -y'},
                    'cost': {'48': 1, '72': 1, '100': 1, '150': 2, '200': 3, '300': 5, '400': 8}},
            'electron_count': {'format': {'kind': 'electron_count', 'mesh': 'N as for ahc', 'mu_eV': '[...]', 'kT_eV': '0-0.025'},
                               'cost': 1}},
        'budget': 32,
        'request_size_limit': 'a request body may not exceed 65536 bytes (round coordinates; about 1000 k-points fit)',
        'conventions': {
            'cell_A': T.CELL.tolist(), 'reciprocal_basis_B_invA': T.B.round(8).tolist(),
            'k_frac': 'coordinates in the reciprocal basis b_i of the POSCAR cell (k_cart = k_frac @ B)',
            'k_cart_invA': 'Cartesian k in 1/Angstrom including 2 pi, axes of the POSCAR cell',
            'energies': 'absolute eV on the scale of the Wannier Hamiltonian (the VASP scale)',
            'mesh': 'Gamma-centred N x N x N/2 mesh of fractional coordinates (i/N, j/N, l/(N/2))',
            'berry_curvature': 'Omega_n = curl_k i<u_n|grad_k u_n> (Angstrom^2), Kubo formula, velocity in the atomic '
                               '(Wannier-centre) gauge with the centres at the atomic positions, as WannierTools 2.6.1 does',
            'ahc': 'sigma_ab = (e^2/hbar) (1/(N_k V)) sum_k sum_n f(E_n - mu) Omega^c_n in S/cm, (a, b, c) cyclic, '
                   'Fermi-Dirac f with the given kT; electrons_per_cell counts the 32 Wannier bands'}}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(cat, indent=1) + '\n')
    for p in ('environment/lab/public.json', 'tests/public.json'):
        d = json.loads((TASK/p).read_text())
        d['solver_budget'] = cat['budget']
        (TASK/p).write_text(json.dumps(d, indent=1) + '\n')


if __name__ == '__main__':
    main()
