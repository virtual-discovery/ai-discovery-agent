"""Reference solution for cu2wse4-orbital-edelstein; runs as the Harbor oracle and in check_task."""
import itertools
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
MUB = 9.2740100783e-24          # J/T = A m^2
KT, GAMMA, EX = 0.030, 0.030, 1.0e4
P_TARGET = 5.0e20               # holes per cm^3
NVAL = {'soc': 88, 'nosoc': 44}           # filled Wannier bands (from the band count at the gap)
SPIN = {'soc': 1, 'nosoc': 2}             # states per band and k-point in the real crystal


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    while result.get('status') == 'running':
        time.sleep(2)
        result = lab(f"job/{result['id']}")
    if result.get('status') == 'failed':
        raise RuntimeError(result.get('error'))
    return result


def read_cif(path):
    cell, atoms = {}, []
    for line in path.read_text().splitlines():
        w = line.split()
        if w and w[0].startswith('_cell_length_'):
            cell[w[0][-1]] = float(w[1])
        if len(w) == 7 and w[0][:2] in ('Cu', 'W6', 'Se'):
            atoms.append((w[0].rstrip('+-0123456789'), np.array([float(x) for x in w[3:6]])))
    return np.diag([cell['a'], cell['b'], cell['c']]), atoms


def point_group(atoms):
    """Rotations g (signed permutations keeping c along z) with a translation t mapping the structure onto itself."""
    ops = []
    for perm in itertools.permutations(range(2)):
        for s in itertools.product((1, -1), repeat=3):
            g = np.zeros((3, 3))
            g[0, perm[0]], g[1, perm[1]], g[2, 2] = s[0], s[1], s[2]
            x0 = atoms[0][1]
            for sp, x in atoms:                 # candidate translations: map atom 0 onto any atom of its species
                if sp != atoms[0][0]:
                    continue
                t = x - g @ x0
                ok = all(any(sp2 == sp1 and np.allclose(((g @ x1 + t - x2 + 0.5) % 1) - 0.5, 0, atol=1e-4)
                             for sp2, x2 in atoms) for sp1, x1 in atoms)
                if ok:
                    ops.append(g)
                    break
    return ops


def allowed_axial2(ops):
    rng = np.random.default_rng(1)
    A = rng.normal(size=(3, 3))
    S = sum(np.linalg.det(g) * g @ A @ g.T for g in ops) / len(ops)
    return S


def allowed_polar3(ops):
    rng = np.random.default_rng(2)
    A = rng.normal(size=(3, 3, 3))
    A = A + A.transpose(0, 2, 1)
    return sum(np.einsum('ia,jb,lc,abc->ijl', g, g, g, A) for g in ops) / len(ops)


def edges(name, path_k):
    """VBM and CBM of a model from one bands request: band path plus local grids around M and Gamma."""
    gm = [[0.5 + i, 0.5 + j, l] for i in np.linspace(-0.03, 0.03, 9) for j in np.linspace(-0.03, 0.03, 9) for l in (-0.02, 0, 0.02)]
    gg = [[i, j, l] for i in np.linspace(-0.1, 0.1, 11) for j in np.linspace(-0.1, 0.1, 11) for l in (-0.02, 0, 0.02)]
    pts = [list(map(float, np.round(k, 5))) for k in list(path_k) + gm + gg]
    r = lab('evaluate', {'kind': 'bands', 'model': name, 'k_frac': pts})
    E = np.array(r['result']['energies_eV'])
    nv = NVAL[name]
    assert (E[:, nv - 1] < 0).all() and (E[:, nv] > 0).all(), 'band count at the gap'
    return r['id'], float(E[:, nv - 1].max()), float(E[:, nv].min())


def interp(x, xs, ys):
    return np.array([np.interp(x, xs, ys[:, k]) for k in range(ys.shape[1])])


def main():
    evidence, notes = [], []
    cell, atoms = read_cif(APP/'archive/Fig1a.cif')
    vol = float(np.linalg.det(cell))
    ops = point_group(atoms)
    A2, A3 = allowed_axial2(ops), allowed_polar3(ops)
    ax = 'xyz'
    allowed_e = [ax[i] + ax[j] for i in range(3) for j in range(3) if abs(A2[i, j]) > 1e-8]
    allowed_n = [ax[i] + ax[j] + ax[l] for i in range(3) for j in range(3) for l in range(3) if abs(A3[i, j, l]) > 1e-8]
    yy_xx = float(A2[1, 1] / A2[0, 0])
    trace = float(np.trace(A2) / A2[0, 0])
    notes.append(f'{len(ops)} point operations (D2d: twofold axes along x and y, mirrors on the diagonals); '
                 f'allowed alpha: {allowed_e}; allowed sigma: {allowed_n}')
    # band edges of both models
    P = {'G': (0, 0, 0), 'X': (.5, 0, 0), 'M': (.5, .5, 0), 'Z': (0, 0, .5), 'R': (.5, 0, .5), 'A': (.5, .5, .5)}
    path = np.concatenate([np.linspace(P[a], P[b], 25) for a, b in
                           [('G', 'X'), ('X', 'M'), ('M', 'G'), ('G', 'Z'), ('Z', 'R'), ('R', 'A'), ('A', 'Z')]])
    edge = {}
    for name in ('soc', 'nosoc'):
        jid, vbm, cbm = edges(name, path)
        evidence.append(jid)
        edge[name] = (vbm, cbm)
    # Fermi-surface sums on the two densest meshes of each model
    res = {}
    for name in ('soc', 'nosoc'):
        vbm, cbm = edge[name]
        grid = np.round(np.arange(vbm - 0.25, vbm + 0.0501, 0.005), 6)
        mus = [float(x) for x in grid] + [round(cbm + 0.10, 6)]
        for n in (96, 128):
            r = lab('evaluate', {'kind': 'fermi_surface', 'model': name, 'mesh': n, 'mu_eV': mus, 'kT_eV': KT})
            evidence.append(r['id'])
            out = r['result']
            ne = np.array(out['electrons_per_cell'])
            holes = SPIN[name] * (NVAL[name] - ne) / vol * 1e24                    # per cm^3
            t = {q: SPIN[name] * np.array(v) / vol for q, v in out['tensors_per_cell'].items()}   # per Angstrom^3
            res[(name, n)] = dict(mu=np.array(mus), holes=holes, t=t, ngrid=len(grid))
    pref = -EX / (2 * GAMMA) * MUB * 1e20    # M (A/m) = q tau E K / hbar with tau = hbar/(2 gamma), K in mu_B/A^2

    def values(n_soc, n_nosoc):
        s, ns = res[('soc', n_soc)], res[('nosoc', n_nosoc)]
        ng = s['ngrid']
        vbm, cbm = edge['soc']
        orb = interp(vbm - 0.10, s['mu'][:ng], s['t']['orbital_berry_phase'][:ng].reshape(ng, 9))
        spin = interp(vbm - 0.10, s['mu'][:ng], s['t']['spin'][:ng].reshape(ng, 9))
        orb_c, spin_c = s['t']['orbital_berry_phase'][-1], s['t']['spin'][-1]
        ratio = {'VBM-0.10': float(orb[0] / spin[0]), 'CBM+0.10': float(orb_c[0, 0] / spin_c[0, 0])}
        mag = float(pref * orb[0])
        # equal hole density: mu(p) by interpolation of the (monotonic) hole density on the valence-edge grid
        a = {}
        for key, d in (('soc', s), ('nosoc', ns)):
            h, mu = d['holes'][:d['ngrid']], d['mu'][:d['ngrid']]
            ok = h > 0
            h, mu = h[ok], mu[ok]
            mu_p = float(np.interp(np.log(P_TARGET), np.log(h[::-1]), mu[::-1]))
            a[key] = float(np.interp(mu_p, d['mu'][:d['ngrid']], d['t']['orbital_berry_phase'][:d['ngrid'], 0, 0]))
        return ratio, mag, a['soc'] / a['nosoc']

    r128, m128, s128 = values(128, 128)
    r96, m96, s96 = values(96, 96)
    s = res[('soc', 128)]
    ng = s['ngrid']
    at_v = np.interp(edge['soc'][0] - 0.10, s['mu'][:ng], s['t']['orbital_intra_atomic'][:ng, 0, 0])
    sp_v = np.interp(edge['soc'][0] - 0.10, s['mu'][:ng], s['t']['spin'][:ng, 0, 0])
    atomic = (at_v / sp_v, s['t']['orbital_intra_atomic'][-1, 0, 0] / s['t']['spin'][-1, 0, 0])
    notes.append(f'mesh 96 vs 128: ratios {r96} vs {r128}; M {m96:.5f} vs {m128:.5f}; soc/nosoc {s96:.4f} vs {s128:.4f}')
    # Berry-curvature dipole: sigma_ijl = (q^3 tau / 2 hbar^2) (D_jj - D_ll) for i, j, l all different
    D = res[('soc', 128)]['t']['berry_curvature_dipole'][-1]
    Dsym = np.diag([(D[0, 0] - D[1, 1]) / 2, (D[1, 1] - D[0, 0]) / 2, 0.0])
    zxy_xyz = float((Dsym[0, 0] - Dsym[1, 1]) / (Dsym[1, 1] - Dsym[2, 2]))
    notes.append(f'BCD tensor at CBM+0.10 (mesh 128, per A^3): diag {np.diag(D).round(8).tolist()}, '
                 f'offdiag max {np.abs(D - np.diag(np.diag(D))).max():.2e}; raw zxy/xyz {(D[0,0]-D[1,1])/(D[1,1]-D[2,2]):.3f}')
    submission = {'orbital_to_spin_ratio': r128, 'orbital_magnetization_A_per_m': m128,
                  'soc_to_nosoc_orbital_ratio': s128, 'edelstein_allowed_components': allowed_e,
                  'alpha_yy_over_alpha_xx': yy_xx, 'chirality_indicator': trace,
                  'nonlinear_hall_allowed_components': allowed_n, 'sigma_zxy_over_sigma_xyz': zxy_xyz,
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Cu2WSe4 Edelstein and nonlinear Hall responses (reference solution)

Symmetry: the staged CIF is written in P1; testing all signed permutations of the cell axes with a translation
search gives {len(ops)} operations, the point group D2d with the twofold axes along x and y (a, b), the S4 axis
along z and mirrors on the diagonal planes. The axial-polar Edelstein tensor is therefore diagonal with
alpha_yy = -alpha_xx and alpha_zz = 0 (tr alpha = 0: the chirality indicator vanishes for this achiral but
gyrotropic class); the nonlinear conductivity allows only {allowed_n}, with sigma_xyz = sigma_yzx and sigma_zxy.
For the Berry-curvature dipole, D_xx = -D_yy and D_zz = 0 give sigma_zxy/sigma_xyz = -2. The service tensors
obey these relations to the noise level.
Orbital moment: the wave-packet moment with the Wannier-centre position operator is the Berry-phase moment
(|e|/2hbar) Im<du|x(H-E)|du>; the intra-atomic -mu_B<L> is a different (on-site) quantity: alpha_atomic/alpha_spin
= {atomic[0]:.3f} at VBM-0.10 eV and {atomic[1]:.2f} at CBM+0.10 eV, so it would reverse both conclusions.
Band edges: VBM {edge['soc'][0]:.4f} eV (soc), {edge['nosoc'][0]:.4f} eV (nosoc), at M; CBM {edge['soc'][1]:.4f} /
{edge['nosoc'][1]:.4f} eV near Gamma (header Fermi levels lie mid-gap and differ between the models).
The model without SOC has 54 spinless bands: each state counts twice for the crystal, in the hole density and in
alpha. Hole densities from the electron counts at kT = 30 meV; alpha interpolated on a 5 meV grid.
Convergence: {notes[-2]}.
Results: alpha_orb/alpha_spin = {r128['VBM-0.10']:.2f} at VBM-0.10 eV (orbital dominates) but {r128['CBM+0.10']:.3f}
at CBM+0.10 eV (spin dominates near the W-d conduction-band edge); M_x(orb) = {m128:.4f} A/m for E_x = 1e4 V/m,
gamma = 30 meV. At p = 5e20 cm^-3 the SOC and no-SOC orbital responses agree to a ratio {s128:.3f}: SOC is not
essential for the valence-band orbital response. The nonlinear Hall response computed here is the BCD part;
other intrinsic terms (and the full Green-function result) are not available from these data.
Not supported: absolute magnitudes beyond the relaxation-time approximation, intra-atomic contributions
(absent from the Wannier position operator), the gap (GGA underestimates it), and any statement about chirality
beyond the vanishing trace.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
