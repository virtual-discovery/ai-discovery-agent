"""Reference solution for mnte-altermagnet-ahc-doping; runs as the Harbor oracle and in check_task.

Uses only the staged archive and lab purchases (no literature values):
1. VBM: a coarse band scan of Wannier band 22 and two local refinements (it lies off the high-symmetry lines).
2. Wannier fidelity: Wannier bands at the 288 DFT k-points against wannier90.eig (valence exact, conduction off by
   up to ~1.6 eV, so electron doping is not determined by the model).
3. Hole densities from the 400 mesh electron count; sigma_xy(mu) on the 400 mesh with a 1 meV smearing and the 300
   mesh as a convergence check; the Neel-reversed model on a small mesh; Berry curvature near the VBM.
4. Symmetry: magnetic operations of the archived state from the structure and MAGMOM (not symmetries.dat).
"""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
ARCH = APP/'archive'
CELL = np.array([[4.14, 0.0, 0.0], [-2.07, 3.585345171667576, 0.0], [0.0, 0.0, 6.71]])
VOL_CM3 = abs(np.linalg.det(CELL)) * 1e-24
NOCC = 22
P_KEYS = {'2.0e21': 2.0e21, '6.0e21': 6.0e21}


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request, evidence):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error', 'job failed'))
    evidence.append(r['id'])
    return r['result']


def grid(center, half, n):
    axes = [np.linspace(c - h, c + h, n) for c, h in zip(center, half)]
    return np.array(np.meshgrid(*axes, indexing='ij')).reshape(3, -1).T


def find_vbm(evidence):
    """Band 22 maximum: 10 x 10 x 8 scan, then three local 9^3 refinements around the best point."""
    f, g = np.arange(10) / 10, np.arange(8) / 8
    k = np.array(np.meshgrid(f, f, g, indexing='ij')).reshape(3, -1).T
    e = np.array(evaluate({'kind': 'bands', 'k_frac': k.round(6).tolist()}, evidence)['energies_eV'])[:, NOCC - 1]
    best, half = k[np.argmax(e)], np.array([0.1, 0.1, 0.125])
    for _ in range(3):
        k = grid(best, half, 9)
        e = np.array(evaluate({'kind': 'bands', 'k_frac': k.round(6).tolist()}, evidence)['energies_eV'])[:, NOCC - 1]
        best, half = k[np.argmax(e)], half / 4
    return float(e.max()), best


def path_and_spin(vbm, kvbm, evidence):
    """Band-22 maximum along the K-G-A-K-H-M-L path, and the spin splitting of the top valence bands at k and C2x k
    (moments along y: an altermagnet has momentum-dependent spin splitting without net moment)."""
    pts = np.array([[1/3, 1/3, 0], [0, 0, 0], [0, 0, .5], [1/3, 1/3, 0], [1/3, 1/3, .5], [.5, 0, 0], [.5, 0, .5]])
    kp = np.concatenate([pts[a] + np.linspace(0, 1, 101)[:, None] * (pts[a + 1] - pts[a]) for a in range(6)])
    path_max = float(np.array(evaluate({'kind': 'bands', 'k_frac': kp.round(6).tolist()}, evidence)['energies_eV'])[:, NOCC - 1].max())
    kc = kvbm @ (2 * np.pi * np.linalg.inv(CELL).T)
    ks = [kc.tolist(), [-kc[0], kc[1], kc[2]], [0.3, 0.0, 0.2], [0.0, 0.3, 0.2]]
    sp = evaluate({'kind': 'bands', 'k_cart_invA': np.round(ks, 6).tolist(), 'spin': True}, evidence)
    sy = np.array(sp['spin_expectation']['y'])[:, NOCC - 2:NOCC]
    return path_max, sy


def wannier_fidelity(evidence):
    win = (ARCH/'2_wannier90/wannier90.win').read_text()
    kd = np.array([[float(x) for x in line.split()[:3]] for line in
                   win.split('begin kpoints')[1].split('end kpoints')[0].splitlines() if line.strip()])
    eig = np.loadtxt(ARCH/'2_wannier90/wannier90.eig')
    dft = eig[:, 2].reshape(len(kd), -1)
    w = np.array(evaluate({'kind': 'bands', 'k_frac': kd.round(8).tolist()}, evidence)['energies_eV'])
    return float(np.abs(w[:, :NOCC] - dft[:, :NOCC]).max()), float(np.abs(w[:, NOCC:NOCC + 2] - dft[:, NOCC:NOCC + 2]).max()), \
        float(dft[:, NOCC].min()), float(w[:, NOCC].min())


def main():
    evidence = []
    vbm, kvbm = find_vbm(evidence)
    val_err, cond_err, dft_cbm, wan_cbm = wannier_fidelity(evidence)
    path_max, sy = path_and_spin(vbm, kvbm, evidence)
    # hole densities from the densest mesh (rigid bands, spinor bands: one electron per band and k-point)
    mus = list(np.round(np.linspace(vbm - 1.0, vbm - 0.2, 401), 6)) + [round(vbm - 0.30, 6)]
    ec = evaluate({'kind': 'electron_count', 'mesh': 400, 'mu_eV': mus, 'kT_eV': 0.001}, evidence)
    holes = (NOCC - np.array(ec['electrons_per_cell'])) / VOL_CM3
    p_030 = float(holes[-1])
    order = np.argsort(holes[:-1])
    mu_of_p = {key: float(np.interp(p, holes[:-1][order], np.array(mus[:-1])[order])) for key, p in P_KEYS.items()}
    # sigma(mu) on the 400 and 300 meshes, 1 meV smearing (zero-temperature limit): the 5 meV sign-change grid,
    # a 2 meV scan from VBM - 1 eV to the VBM and the chemical potentials of the requested hole densities
    sgrid = list(np.round(vbm - np.arange(1.0, 0.25 - 1e-9, -0.005), 6))
    scan = list(np.round(np.linspace(vbm - 1.0, vbm, 340), 6))
    req = sgrid + scan + [round(m, 6) for m in mu_of_p.values()]
    res = {n: evaluate({'kind': 'ahc', 'mesh': n, 'mu_eV': req, 'kT_eV': 0.001}, evidence) for n in (400, 300, 200)}
    s400 = {c: np.array(v) for c, v in res[400]['sigma_S_per_cm'].items()}
    s300 = {c: np.array(v) for c, v in res[300]['sigma_S_per_cm'].items()}
    ng, nsc = len(sgrid), len(scan)
    mu = np.array(sgrid)
    sx = s400['xy'][:ng]
    flips = np.where(np.sign(sx[1:]) * np.sign(sx[:-1]) < 0)[0]
    mu0 = [mu[j] - sx[j] * (mu[j + 1] - mu[j]) / (sx[j + 1] - sx[j]) for j in flips]
    e400 = np.array(res[400]['electrons_per_cell'])[:ng]
    p0 = float(np.interp(min(mu0), mu, (NOCC - e400) / VOL_CM3)) if mu0 else float('nan')
    sigma_p = {key: float(s400['xy'][ng + nsc + i]) for i, key in enumerate(mu_of_p)}
    sigma_p_300 = {key: float(s300['xy'][ng + nsc + i]) for i, key in enumerate(mu_of_p)}
    sigma_p_200 = {key: float(res[200]['sigma_S_per_cm']['xy'][ng + nsc + i]) for i, key in enumerate(mu_of_p)}
    sx3 = s300['xy'][:ng]
    flips300 = int((np.sign(sx3[1:]) * np.sign(sx3[:-1]) < 0).sum())
    raw_forbidden = {c: float(np.abs(s400[c][ng:ng + nsc]).max()) for c in ('yz', 'zx')}
    # Neel reversal on a small mesh, Berry curvature near the VBM pocket
    test_mu = [round(vbm - 0.8, 6), round(vbm - 0.5, 6)]
    plus = evaluate({'kind': 'ahc', 'mesh': 48, 'mu_eV': test_mu, 'kT_eV': 0.001, 'neel': '+y'}, evidence)['sigma_S_per_cm']['xy']
    minus = evaluate({'kind': 'ahc', 'mesh': 48, 'mu_eV': test_mu, 'kT_eV': 0.001, 'neel': '-y'}, evidence)['sigma_S_per_cm']['xy']
    kloc = grid(kvbm, [0.03, 0.03, 0.03], 7)
    bc = evaluate({'kind': 'berry_curvature', 'k_frac': kloc.round(6).tolist()}, evidence)
    om = np.array(bc['omega_A2']['z'])
    # a kz = 0.42 plane through the VBM pocket: where the curvature of the top two bands is concentrated
    kpl = np.array([[a, b, kvbm[2]] for a in np.linspace(-0.15, 0.15, 21) for b in np.linspace(-0.15, 0.15, 21)])
    bp = evaluate({'kind': 'berry_curvature', 'k_frac': kpl.round(6).tolist()}, evidence)
    gap = np.array(bp['energies_eV'])[:, NOCC - 1] - np.array(bp['energies_eV'])[:, NOCC - 2]
    hot = np.abs(np.array(bp['omega_A2']['z'])[:, NOCC - 1])
    hot_gap = float(np.median(gap[hot >= np.quantile(hot, 0.95)])) * 1000
    submission = {'sigma_xy_S_per_cm': sigma_p,
                  'sigma_xy_sign_changes': int(len(flips)),
                  'sign_change_hole_density_cm3': p0,
                  'hole_density_0p30_below_vbm_cm3': p_030,
                  'forbidden_components_max_abs_S_per_cm': {'yz': 0.0, 'zx': 0.0},
                  'allowed_hall_components': {'neel_y': ['xy'], 'neel_x': []},
                  'electron_doped_sigma_xy_S_per_cm': None,
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# MnTe anomalous Hall conductivity (reference solution)

Symmetry. The archived state has Mn moments +-y (MAGMOM 0 4 0 / 0 -4 0) on the two Mn sublattices, which are related
by the 6_3 screw (C2z with c/2), not by a translation or inversion. The magnetic point group keeps E, I, C2z, m_z
(with c/2) as unitary operations and C2x, C2y, m_x, m_y combined with time reversal (C2y and m_y with c/2), as the
labels in wannsymm.out state. The rotation matrices printed in symmetries.dat are not symmetries of this structure
(several do not map the atoms onto each other), so they were not used. C2z kills sigma_yz and sigma_zx; C2x T and
C2y T keep sigma_xy: only sigma_xy is allowed (it has the symmetry of a magnetization along z, which is why a weak
moment along c is allowed and appears in the DFT run). For moments along x, C2x is unitary and C2y comes with
c/2 and no time reversal; together with C2z they forbid all three components.

Energy reference. The VBM of the Wannier model is {vbm:.4f} eV at k_frac = {np.round(kvbm, 3).tolist()}, off the
high-symmetry lines (the K-G-A-K-H-M-L path maximum is {path_max:.4f} eV; the VASP E-fermi of 6.2098 eV lies in the
gap and the scf grid misses the maximum). The top valence bands carry <sigma_y> = {np.round(sy, 2).tolist()} at the
VBM, its C2x image and two generic k: momentum-dependent spin splitting of a compensated (altermagnetic) state. Hole densities come from the electron count of the 32 spinor bands (one electron per band and
k-point; 22 filled bands): p(VBM - 0.30 eV) = {p_030:.4g} cm^-3.

Conductivity. sigma_xy from the Kubo formula on the 400 x 400 x 200 mesh with 1 meV smearing: {sigma_p} S/cm at the
requested hole densities (300 mesh: {sigma_p_300}; 200 mesh: {sigma_p_200}, not converged at 2e21). On the 5 meV grid from 1.0 to 0.25 eV below the VBM sigma_xy changes
sign {len(flips)} time(s) (300 mesh: {flips300}), at p = {p0:.4g} cm^-3; coarse meshes show spurious extra sign changes
and spikes from sharp Berry-curvature hot spots of the SOC-split top valence bands (Omega_z of band 22 near the VBM
spans {om[:, NOCC - 1].min():.0f} to {om[:, NOCC - 1].max():.0f} A^2; in the kz plane of the VBM the strongest 5% of
|Omega_z| of band 22 sit where it lies a median {hot_gap:.0f} meV above band 21). Closer to the band edge (realistic
dopings below ~1e21 cm^-3) the mesh results are not converged and are not reported.
The archived Hamiltonian is not symmetrized (the authors symmetrized it with WannSymm before WannierTools); it gives
nonzero sigma_yz and sigma_zx of up to {raw_forbidden['yz']:.2f} and {raw_forbidden['zx']:.2f} S/cm on the 400 mesh,
numerical symmetry breaking: the physical values are zero.

Neel reversal. The time-reversed model gives sigma_xy = {np.round(minus, 2).tolist()} vs {np.round(plus, 2).tolist()} S/cm:
sigma_xy is odd in the Neel vector. In a multidomain film opposite-Neel domains cancel (the three 120-degree
domain orientations give the same sign); a net Hall signal needs an imbalance of +L and -L domains, which the weak
net moment along c (allowed by the same symmetry as sigma_xy) can select in a field (the paper's macrospin model).
The tiny DFT net moment (OSZICAR mag_z ~ 6e-4 muB/cell) cannot itself produce ~100 S/cm: the AHC is a property
of the compensated altermagnetic order.

Not supported. The Wannier conduction bands do not reproduce DFT (up to {cond_err:.2f} eV off on the fitting grid;
DFT CBM {dft_cbm:.3f} eV, Wannier {wan_cbm:.3f} eV; valence bands agree to {val_err * 1000:.2f} meV), so electron-doped
values are not determined. The near-edge (low-doping) AHC is not converged on the available meshes; the WannSymm-
symmetrized Hamiltonian the authors used is not archived; domain populations are outside the model.
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
