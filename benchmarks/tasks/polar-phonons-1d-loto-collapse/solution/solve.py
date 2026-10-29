"""Reference solution for polar-phonons-1d-loto-collapse; runs as the Harbor oracle and in check_task.

Protocol: 1D-boundary-condition data only for the 1D dispersion. The force constants of each tube are
Fourier-interpolated along the axis after subtracting a 1D long-range dipole term (filled dielectric cylinder of
radius t holding the archived polarizability per length), which is added back at the target q. The radius is
calibrated on the released BN(4,4) tube (explicit small-q runs) as a ratio to the atomic-shell radius and
transferred to BN(5,5) and BN(6,6). The 3D-periodic q->0 LO follows from the Gamma run of the 3D_code folder
with the standard nonanalytic term of the periodic supercell.
"""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.special import i0, i1, k0, k1

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
RY_TO_CMM1 = 109737.31570111268
E2 = 2.0
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(file):
    r = lab('evaluate', {'kind': 'dfpt', 'file': file})
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    EVIDENCE.append(r['id'])
    return parse_dyn(r['result']['text'])


def parse_dyn(text):
    lines = text.split('\n')
    h = lines[2].split()
    ntyp, nat, ibrav = int(h[0]), int(h[1]), int(h[2])
    celldm = [float(x) for x in h[3:9]]
    i, lattice = 3, None
    if ibrav == 0:
        lattice = np.array([[float(x) for x in lines[i + 1 + k].split()] for k in range(3)])
        i += 4
    species = {}
    for _ in range(ntyp):
        m = re.match(r"\s*(\d+)\s+'([^']*)'\s+([\d.Ee+-]+)", lines[i])
        species[int(m.group(1))] = float(m.group(3))
        i += 1
    ityp, tau = [], []
    for _ in range(nat):
        p = lines[i].split()
        ityp.append(int(p[1]))
        tau.append([float(x) for x in p[2:5]])
        i += 1
    d = {'nat': nat, 'ibrav': ibrav, 'celldm': celldm, 'lattice': lattice, 'tau_alat': np.array(tau),
         'mass': np.array([species[t] for t in ityp]), 'blocks': []}
    while i < len(lines):
        line = lines[i]
        if 'Dynamical  Matrix in cartesian axes' in line:
            q = [float(x) for x in lines[i + 2].split('(')[1].split(')')[0].split()]
            i += 4
            D = np.zeros((nat, 3, nat, 3), complex)
            for _ in range(nat * nat):
                a, b = [int(x) for x in lines[i].split()]
                for r in range(3):
                    v = [float(x) for x in lines[i + 1 + r].split()]
                    D[a - 1, r, b - 1, :] = np.array(v[0::2]) + 1j * np.array(v[1::2])
                i += 4
            d['blocks'].append((np.array(q), D))
            continue
        if 'Dielectric Tensor' in line:
            d['eps'] = np.array([[float(x) for x in lines[i + 2 + r].split()] for r in range(3)])
            i += 5
            continue
        if 'Effective Charges E-U' in line:
            Z = np.zeros((nat, 3, 3))
            i += 2
            for a in range(nat):
                i += 1
                for r in range(3):
                    Z[a, r] = [float(x) for x in lines[i].split()]
                    i += 1
            d['Z'] = Z
            continue
        if 'Diagonalizing the dynamical matrix' in line:
            break
        i += 1
    return d


def nq_of(folder):
    ph = CATALOG['folders'][folder]['input_files']['ph.in']
    return int(re.search(r'nq3\s*=\s*(\d+)', ph).group(1))


def system(folder, grid):
    g = grid[0]
    alat = g['celldm'][0]
    lat = g['lattice'] * alat if g['ibrav'] == 0 else np.diag([alat, alat, alat * g['celldm'][2]])
    c, nq = lat[2, 2], nq_of(folder)
    Dq = {}
    for d in grid:
        for q, D in d['blocks']:
            Dq[int(round(q[2] * c / alat * nq)) % nq] = D          # q in 2pi/alat -> index on the nq-point axis grid
    assert len(Dq) == nq, (folder, sorted(Dq))
    return {'nq': nq, 'alat': alat, 'lat': lat, 'c': c, 'Dq': Dq, 'tau': g['tau_alat'] * alat, 'Z': g['Z'],
            'eps': g['eps'], 'mass': g['mass'], 'nat': g['nat'], 'area': lat[0, 0] * lat[1, 1]}


def delta(x, eps):
    return 2 * i1(x) * k1(x) / (x * (i0(x) * k1(x) + eps * i1(x) * k0(x)))


def na_1d(s, qz, t):
    """1D dipole-dipole term (Ry/bohr^2): filled cylinder of radius t, eps1d from the polarizability per length,
    Ewald-damped (alpha = (2pi/alat)^2) so that the remainder is short ranged; G along the axis included."""
    eps1d = s['area'] * (s['eps'][2, 2] - 1) / (np.pi * t ** 2)
    D = np.zeros((s['nat'], 3, s['nat'], 3), complex)
    for m in range(-3, 4):
        g = qz + m * 2 * np.pi / s['c']
        gt = g * s['alat'] / (2 * np.pi)
        if abs(g) < 1e-12 or gt ** 2 / 4 > 14:
            continue
        W = (1 - delta(abs(g) * t, eps1d)) / (np.pi * t ** 2 * eps1d * g ** 2) * np.exp(-gt ** 2 / 4)
        zg = g * s['Z'][:, 2, :]
        ph = np.exp(1j * g * (s['tau'][:, None, 2] - s['tau'][None, :, 2]))
        D += 4 * np.pi * E2 / s['c'] * W * np.einsum('ai,bj,ab->aibj', zg, zg, ph)
    return D


class Interp:
    def __init__(self, s, t):
        self.s, self.t = s, t
        nq, c = s['nq'], s['c']
        Dsr = {}
        for m, D in s['Dq'].items():
            qz = m * 2 * np.pi / (nq * c)
            qz = qz - 2 * np.pi / c if qz > np.pi / c + 1e-9 else qz
            Dsr[m] = D - na_1d(s, qz, t)
        self.C = [sum(Dsr[m] * np.exp(2j * np.pi * m * n / nq) for m in range(nq)) / nq for n in range(nq)]
        L = nq * c
        self.ns = np.arange(-2 * nq, 2 * nq + 1)
        self.w = np.zeros((s['nat'], s['nat'], len(self.ns)))
        for a in range(s['nat']):
            for b in range(s['nat']):
                r = self.ns * c + s['tau'][a, 2] - s['tau'][b, 2]
                d0, dm, dp = np.abs(r), np.abs(r - L), np.abs(r + L)
                inside = (dm >= d0 - 1e-6) & (dp >= d0 - 1e-6)
                cnt = 1 + (np.abs(dm - d0) < 1e-6) + (np.abs(dp - d0) < 1e-6)
                self.w[a, b] = np.where(inside, 1.0 / cnt, 0.0)

    def D(self, qz):
        s = self.s
        D = np.zeros((s['nat'], 3, s['nat'], 3), complex)
        for i, n in enumerate(self.ns):
            if self.w[:, :, i].any():
                D += np.einsum('ab,aibj->aibj', self.w[:, :, i], self.C[n % s['nq']]) * np.exp(-1j * qz * n * s['c'])
        return D + na_1d(s, qz, self.t)


def lo(s, D):
    M = s['mass']
    Dm = (D / np.sqrt(M[:, None, None, None] * M[None, None, :, None])).reshape(3 * s['nat'], -1)
    w2, v = np.linalg.eigh(0.5 * (Dm + Dm.conj().T))
    f = np.sign(w2) * np.sqrt(np.abs(w2)) * RY_TO_CMM1
    Pz = np.abs(np.einsum('aj,ajn->n', s['Z'][:, 2, :], v.reshape(s['nat'], 3, -1) / np.sqrt(M)[:, None, None]))
    return float(f[np.argmax(Pz)])


def radius(s):
    xy = s['tau'][:, :2]
    return float(np.mean(np.linalg.norm(xy - xy.mean(0), axis=1)))


def qfrac(s, d):
    return d['blocks'][0][0][2] * 2 * np.pi / s['alat'] / (np.pi / s['c'])


def main():
    # Released BN(4,4): 1D grid plus three explicit small-q 1D runs to calibrate the cylinder radius.
    f44 = 'BN_nanotubes/44/1D_code'
    s44 = system(f44, [buy(f'{f44}/dyn{k}') for k in range(1, 8)])
    explicit44 = [buy(f'{f44}/{n}') for n in ('dyn25', 'dyn8', 'dyn9')]
    pts44 = [(qfrac(s44, d), lo(s44, d['blocks'][0][1])) for d in explicit44]
    R44 = radius(s44)
    tmax = np.sqrt(s44['area'] * (s44['eps'][2, 2] - 1) / np.pi)       # eps1d >= 1
    scan = []
    for t in np.arange(R44, tmax, 0.05):
        ip = Interp(s44, t)
        scan.append((np.sqrt(np.mean([(lo(s44, ip.D(q * np.pi / s44['c'])) - f) ** 2 for q, f in pts44])), t))
    rms44, t44 = min(scan)
    ratio = t44 / R44
    # Hidden tubes: 1D grids, the explicit 0.15 run as a check, and predictions.
    out, checks, info = {}, {}, {}
    targets = {'BN55': ['0.01', '0.025', '0.05', '0.075', '0.125'], 'BN66': ['0.01', '0.025', '0.05', '0.075', '0.1', '0.125']}
    for tube, folder in (('BN55', 'BN_nanotubes/55/1D_code'), ('BN66', 'BN_nanotubes/66/1D_code')):
        s = system(folder, [buy(f'{folder}/dyn{k}') for k in range(1, 8)])
        t = ratio * radius(s)
        ip = Interp(s, t)
        chk = buy(f'{folder}/longwave/dyn15')
        checks[tube] = (lo(s, chk['blocks'][0][1]), lo(s, ip.D(qfrac(s, chk) * np.pi / s['c'])))
        out[tube] = {q: lo(s, ip.D(float(q) * np.pi / s['c'])) for q in targets[tube]}
        info[tube] = {'t': t, 'R': radius(s), 'gamma': lo(s, s['Dq'][0])}
    # 3D-periodic BN(5,5): standard nonanalytic term of the supercell for q -> 0 along the axis.
    g3 = buy('BN_nanotubes/55/3D_code/dyn1')
    alat = g3['celldm'][0]
    lat3 = np.diag([alat, alat, alat * g3['celldm'][2]])
    s3 = {'nat': g3['nat'], 'mass': g3['mass'], 'Z': g3['Z']}
    zq = g3['Z'][:, 2, :]
    na3 = 4 * np.pi * E2 / np.linalg.det(lat3) * np.einsum('ai,bj->aibj', zq, zq) / g3['eps'][2, 2]
    lo3 = lo(s3, g3['blocks'][0][1] + na3)
    lo3_gamma = lo(s3, g3['blocks'][0][1])
    submission = {'BN55_LO_1D_cm-1': out['BN55'], 'BN66_LO_1D_cm-1': out['BN66'], 'BN55_LO_3D_q0_cm-1': lo3,
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    rows = '\n'.join(f'- {tube}: ' + ', '.join(f'q={q}: {v:.1f}' for q, v in out[tube].items()) for tube in out)
    (APP/'report.md').write_text(f"""# Polar LO phonons of BN nanotubes in 1D boundary conditions (reference solution)

Data: the authors' 1D-Coulomb-cutoff DFPT (1D_code) for the 1D dispersion; the 3D-periodic runs (3D_code) only for
the artefact. q is converted from the archived 2pi/alat units (alat = lateral cell) to fractions of |Gamma-Z| = pi/L.
The LO is the mode with the largest longitudinal mode effective charge (Born charges of the Gamma run).

Physics. In 3D the macroscopic field gives a finite LO-TO splitting at Gamma; in 2D the polar shift collapses
linearly in q; in 1D the field lines escape into vacuum and omega_LO^2 - omega_Gamma^2 ~ (2 e^2/L) S^2 q^2 ln(1/(q t))
as q -> 0 (S = mode effective charge per sqrt(mass)), so the shift vanishes at Gamma with a logarithmic overbending.
The leading term does not depend on the dielectric screening or the radius; both enter at the next order through
the cylinder model, in which only the polarizability per length A(eps_zz - 1)/(4 pi) is physical
(eps1d = A(eps_zz - 1)/(pi t^2)).

Method. Short-range force constants from the 12-point 1D grid after subtracting the Ewald-damped 1D cylinder
kernel (full Bessel form, not its small-qt expansion), Fourier interpolation along the axis with Wigner-Seitz
weights, and the kernel added back at the target q. Plain Fourier interpolation misses the q^2 ln q
nonanalyticity (the archived BN(4,4) "1D" interpolated dispersion is of this kind: it lies 30 cm^-1 below the
explicit 1D runs at 0.05-0.10 |Gamma-Z|). The cylinder radius was calibrated on BN(4,4) explicit runs
({', '.join(f'{q:.3f}' for q, _ in pts44)} |Gamma-Z|): t = {t44:.2f} bohr = {ratio:.2f} x atomic-shell radius (rms {rms44:.2f} cm^-1),
transferred as t = {ratio:.2f} R: BN(5,5) t = {info['BN55']['t']:.2f}, BN(6,6) t = {info['BN66']['t']:.2f} bohr. The atomic-shell
radius itself is too small (the electronic polarization extends beyond the atoms). Check at the purchasable
0.15 |Gamma-Z| runs: BN(5,5) DFPT {checks['BN55'][0]:.1f} vs {checks['BN55'][1]:.1f}; BN(6,6) {checks['BN66'][0]:.1f} vs {checks['BN66'][1]:.1f} cm^-1.

Predictions (cm^-1; Gamma LO {info['BN55']['gamma']:.1f} and {info['BN66']['gamma']:.1f}):
{rows}

A Gamma-only model (constant short-range part) works only for q t << 1; beyond ~0.075 |Gamma-Z| the mechanical
dispersion matters, which is why the grids are needed. Interpolation-method spread (Ewald split, radius within
its calibrated range) is about 1-3 cm^-1 at 0.05-0.125 |Gamma-Z|.

3D-periodic artefact. With periodic images the q -> 0 macroscopic field is that of a 3D crystal of tubes:
omega_LO(q->0) = {lo3:.1f} cm^-1 for BN(5,5) (Gamma analytic {lo3_gamma:.1f}), set by the vacuum-diluted eps_zz of the
supercell and its volume; it shrinks only slowly with vacuum (eps_zz - 1 and 1/Omega scale as 1/d^2) and is absent
in the true 1D limit. For q above ~2 pi/d the 3D and 1D data coincide.

Spectroscopy. Raman/IR probe q ~ 2 pi/lambda, far inside the logarithmic region: an isolated tube shows an LO line
near its Gamma value, red-shifted from the bulk-like LO, and the shift depends on radius and wavelength; arrays
with spacing below the probing wavelength recover a 3D-like shift. Not supported: absolute agreement with
experiment (PBE, pseudopotentials), a unique effective radius (the fit is shallow), the hollow-tube geometry.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
