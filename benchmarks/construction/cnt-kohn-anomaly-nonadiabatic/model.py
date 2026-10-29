"""Readers and estimators for cnt-kohn-anomaly-nonadiabatic (construction side).

The archived self-energy tables (Fig2/PiT/<T>) hold, per mode and frequency w, the columns
(1-based) 6 = Re Pi(w), 7 = Re Pi(0) (independent of w) and 8 = Im Pi(w), in meV; their header
names the columns differently (see SOURCE.md). Renormalized frequencies follow the complex-pole
convention of the instruction (the source's Eqs. (6)-(7)).
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
CACHE = Path(os.environ.get('CNT_CACHE', '/cache'))
POOL = CACHE/'pool'
KB_EV = 8.617333262145179e-05
RY_EV = 13.605693122994
CM1_MEV = 0.12398419843320026
T_SIGMA = 0.01 * RY_EV / KB_EV                 # degauss = 0.01 Ry as a Fermi-Dirac temperature (1578.87 K)
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import selfenergy as SE                      # noqa: E402


class Tables:
    """All 37 archived self-energy tables."""

    def __init__(self):
        z = np.load(POOL/'pit.npz')
        self.temps = [float(t) for t in z['temps']]
        self.w = z['w_eV'] * 1000.0
        self.wq = z['wq_eV'] * 1000.0
        self.re, self.static, self.im = z['col6'], z['col7'][:, :, 0], z['col8']

    def get(self, nu, t):
        i = min(range(len(self.temps)), key=lambda j: abs(self.temps[j] - t))
        assert abs(self.temps[i] - t) < 1e-6, t
        return self.wq[i, nu - 1], self.w, self.re[i, nu - 1], self.im[i, nu - 1], self.static[i, nu - 1]


def pole(w, re, im, wq, pdc, gamma_term=True, scale=1.0):
    """Largest positive root of Omega^2 - gamma^2 = wq^2 + 2 wq Re[Pi(Omega) - pdc], gamma = -wq Im Pi(Omega) / Omega.

    Pi is linearly interpolated on the frequency grid w. Returns (Omega, gamma) or (None, None)."""
    def rhs(o):
        r, i = np.interp(o, w, re) * scale, np.interp(o, w, im) * scale
        g = -wq * i / o
        return wq ** 2 + (g ** 2 if gamma_term else 0.0) + 2 * wq * (r - pdc * scale), g
    grid = np.linspace(0.05, 1.6 * wq, 4000)
    vals = np.array([o ** 2 - rhs(o)[0] for o in grid])
    roots = [brentq(lambda o: o ** 2 - rhs(o)[0], grid[i], grid[i + 1]) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0]
    if not roots:
        return None, None
    return roots[-1], rhs(roots[-1])[1]


def authors_iteration(w, wq, re, im, pdc, tol=1.0):
    """The authors' fixed-point iteration (fig2_d_f_h.py): nearest grid point, 1 meV tolerance, signed sqrt."""
    cur = last = wq
    for _ in range(100000):
        idx = int(np.argmin(np.abs(w - cur)))
        o2 = wq ** 2 + im[idx] ** 2 * wq ** 2 / cur ** 2 + 2 * wq * (re[idx] - pdc)
        cur = float(np.sqrt(abs(o2)) * np.sign(o2))
        if abs(cur - last) < tol:
            return cur, float(-wq * im[idx] / cur)
        last = cur
    return float('nan'), float('nan')


def instability_function(wq, static, im_slope, pdc, scale=1.0, gamma_term=True):
    """Omega^2 in the limit Omega -> 0+: wq^2 + gamma0^2 + 2 wq (Pi(0,T) - pdc), gamma0 = -wq dIm Pi/dw(0)."""
    g0 = -wq * im_slope * scale
    return wq ** 2 + (g0 ** 2 if gamma_term else 0.0) + 2 * wq * (static - pdc) * scale


def im_slope(w, im):
    j = int(np.where(w > 0)[0][0])
    return (im[j] - im[j - 1]) / (w[j] - w[j - 1])


def crossing(ts, fs):
    """Highest temperature at which fs changes sign from negative (unstable) to positive, linear in T."""
    ts, fs = list(ts), list(fs)
    for i in range(len(ts) - 2, -1, -1):
        if fs[i] < 0 <= fs[i + 1]:
            return ts[i] + (ts[i + 1] - ts[i]) * (-fs[i]) / (fs[i + 1] - fs[i])
    return None


def couplings(nu, processes=('k+q',)):
    z = np.load(POOL/'coupling.npz')
    out = []
    for p in processes:
        tag = 'p' if p == 'k+q' else 'm'
        out.append(SE.Coupling(z['k'], z[f'e_{nu}_{tag}'], z[f'eq_{nu}_{tag}'], z[f'g_{nu}_{tag}'], z[f'wq_{nu}_{tag}']))
    return out


def threshold_model(a, omega_mev):
    """Run the archived coupling_calculator.c algorithm (rebuilt with gcc, main() parametrised) for one (a, omega)."""
    src = json.loads((POOL/'archive.json').read_text())['files']['Fig1/coupling_calculator.c']
    main = '''int main(int argc, char **argv) {
    double mass = atof(argv[1]), freq = atof(argv[2]);
    double omega_scan[OMEGA_POINTS]; double gamma = 1.0; double g = 100.0;
    for (int i = 0; i < OMEGA_POINTS; i++) omega_scan[i] = freq * i / (OMEGA_POINTS - 1.0);
    for (int i = OMEGA_POINTS - 1; i >= 0; i--) {
        gamma = compute_linewidth(gamma, g, mass, freq, omega_scan[i]);
        g = compute_coupling(gamma, g, mass, N_KPOINTS, freq, omega_scan[i], TEMPERATURE, SIGMA);
        printf("%.6f %.6f\\n", omega_scan[i], g);
        if (isnan(gamma)) gamma = 1.0;
        if (isnan(g)) g = 100.0;
    }
    return 0;
}
'''
    with tempfile.TemporaryDirectory() as d:
        c = Path(d)/'cc.c'
        c.write_text(src[:src.index('int main(void)')] + main)
        subprocess.run(['gcc', '-O2', str(c), '-o', str(Path(d)/'cc'), '-lm'], check=True)
        out = subprocess.run([str(Path(d)/'cc'), repr(float(a)), repr(float(omega_mev))], capture_output=True,
                             text=True, check=True).stdout.split('\n')
    vals = [tuple(map(float, l.split())) for l in out if l.strip()]
    finite = [g for o, g in vals if o > 0 and np.isfinite(g)]
    return finite[-1]
