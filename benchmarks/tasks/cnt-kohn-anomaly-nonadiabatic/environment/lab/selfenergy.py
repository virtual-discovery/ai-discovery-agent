"""Phonon self-energy of the (3,3) carbon nanotube at the nesting wavevector Q (numpy evaluator).

Copied unchanged into the task's lab as environment/lab/selfenergy.py by make_task.py.

Data: the authors' archived EPW output on a uniform 2000-point k grid (Fig2/g/g+q.* and g-q.*): the
two Wannier bands at k and at k+Q (or k-Q), ordered by energy at each k, and |g_mn,nu(k,Q)| in meV for
the band pairs (n, m) and the modes nu = 1, 20, 23.

For a process p (k -> k+Q or k -> k-Q):

    Pi_nu(w, T) = (2 / N_k) sum_k sum_(n,m) |g_mn,nu(k)|^2 [f(e_m,k+pQ) - f(e_n,k)] / [e_m,k+pQ - e_n,k - w - i eta]

summed over the selected processes and band pairs; f is the Fermi-Dirac function at temperature T with
the chemical potential fixed at the archived EPW value. The bands and |g|^2 are linearly interpolated
from the archived grid onto N_k uniform points after the two bands have been followed through their
crossings (branch tracking), so that the interpolation does not mix bands; the band-pair selection
refers to the archive's energy ordering at each interpolated k.
"""
import numpy as np

KB_EV = 8.617333262145179e-05      # Boltzmann constant (eV/K), as in the authors' scripts
MU_EV = -1.0380                    # chemical potential: fermi_energy of the archived EPW self-energy input
MODES = (1, 20, 23)
PROCESSES = {'k+q': '+', 'k-q': '-'}


def track(e):
    """Follow the two energy-ordered bands e (2, N) through crossings on a periodic grid.

    Returns perm (N, 2): branch b at grid index i is the energy-ordered band perm[i, b]."""
    n = e.shape[1]
    perm = np.zeros((n, 2), int)
    perm[0] = perm[1] = (0, 1)
    for i in range(2, n):
        p1, p2 = perm[i - 1], perm[i - 2]
        pred = np.array([2 * e[p1[b], i - 1] - e[p2[b], i - 2] for b in (0, 1)])
        keep = abs(pred[0] - e[0, i]) + abs(pred[1] - e[1, i])
        swap = abs(pred[0] - e[1, i]) + abs(pred[1] - e[0, i])
        perm[i] = (0, 1) if keep <= swap else (1, 0)
    return perm


class Coupling:
    """Archived bands and |g| for one mode and one process, in branch-tracked form."""

    def __init__(self, k, e, eq, g, wq):
        self.k, self.wq = np.asarray(k, float), float(wq)
        e, eq, g = np.asarray(e, float), np.asarray(eq, float), np.asarray(g, float)
        idx = np.arange(len(self.k))
        pi, pf = track(e), track(eq)
        self.e = np.array([e[pi[:, a], idx] for a in (0, 1)])            # initial-state branches (eV)
        self.eq = np.array([eq[pf[:, b], idx] for b in (0, 1)])          # final-state branches (eV)
        self.g2 = np.array([[g[pi[:, a], pf[:, b], idx] ** 2 for b in (0, 1)] for a in (0, 1)])  # meV^2

    def on(self, kd, y):
        kk = np.r_[self.k, self.k[0] + 1.0]
        return np.interp(kd, kk, np.r_[y, y[0]])


def fermi(x, t):
    return 0.5 * (1.0 - np.tanh(x / (2.0 * KB_EV * t)))


def self_energy(couplings, temperature_k, omegas_mev, nk, eta_mev, band_pairs=((1, 1), (1, 2), (2, 1), (2, 2)),
                mu_ev=MU_EV):
    """Pi(w) in meV at the real frequencies omegas_mev (meV), summed over the given Coupling objects.

    band_pairs are (n, m) in the archive's energy ordering (1 = lower band). Returns (re, im) arrays."""
    w = np.asarray(omegas_mev, float) / 1000.0
    eta = eta_mev / 1000.0
    kd = -0.5 + np.arange(nk) / nk
    out = np.zeros(len(w), complex)
    wanted = {tuple(p) for p in band_pairs}
    for c in couplings:
        ei = np.array([c.on(kd, c.e[a]) for a in (0, 1)])
        ef = np.array([c.on(kd, c.eq[b]) for b in (0, 1)])
        rank_i = np.where(ei[0] <= ei[1], 1, 2), np.where(ei[0] <= ei[1], 2, 1)   # energy-order label of branch a
        rank_f = np.where(ef[0] <= ef[1], 1, 2), np.where(ef[0] <= ef[1], 2, 1)
        fi = fermi(ei - mu_ev, temperature_k)
        ff = fermi(ef - mu_ev, temperature_k)
        for a in (0, 1):
            for b in (0, 1):
                mask = np.zeros(nk, bool)
                for n, m in wanted:
                    mask |= (rank_i[a] == n) & (rank_f[b] == m)
                num = (ff[b] - fi[a]) * c.on(kd, c.g2[a, b])
                sel = mask & (np.abs(num) > 1e-14)
                if not sel.any():
                    continue
                num_s, de_s = num[sel], (ef[b] - ei[a])[sel]
                chunk = max(1, int(4e6 // len(num_s)))          # bounds the temporary array to ~64 MB
                for j in range(0, len(w), chunk):
                    ww = w[j:j + chunk]
                    out[j:j + chunk] += (num_s[None, :] / (de_s[None, :] - ww[:, None] - 1j * eta)).sum(axis=1)
    out *= 2.0 / nk / 1000.0          # spin, k average; meV^2 / eV -> meV
    return out.real, out.imag
