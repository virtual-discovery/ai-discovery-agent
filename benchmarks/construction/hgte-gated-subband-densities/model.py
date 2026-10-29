"""Readers for the archived kdotpy runs (pool.npz from extract.py), shared by make_task.py.

Subband densities follow the authors' FigureS3.ipynb: the band-resolved integrated DOS of the main
self-consistent run evaluated at that run's E_F(n) ([ef] in the .kpout); hole pockets from the wedge runs
(kphi 105-165 and 165-285 deg), evaluated at the MAIN run's E_F(n) and divided by 360/(wedge width),
because kdotpy extrapolates a partial angular range to the full circle as if it were isotropic.
Densities are returned in 1e11 cm^-2 (IDOS in nm^-2 times 1e3); holes positive.
"""
import json
import numpy as np

BANDS = [-3, -2, -1, 1, 2, 3, 4, 5, 6]


class Pool:
    def __init__(self, npz='/scratch/pool.npz', meta='/scratch/pool_meta.json'):
        self.P = np.load(npz)
        self.meta = json.load(open(meta))

    def key(self, s, kind, n):
        return f'{s}_{kind}_{n:+.4f}'

    def n_values(self, s, kind):
        return sorted(v['n_nm2'] for v in self.meta['calcs'].values() if v['sample'] == s and v['kind'] == kind)

    def kp(self, s, kind, n):
        return self.meta['calcs'][self.key(s, kind, n)]['kpout']

    def idos_at(self, s, kind, n, energy):
        k = self.key(s, kind, n)
        E, I = self.P[f'E_{k}'], self.P[f'I_{k}']
        i = int(np.abs(E - energy).argmin())
        return {b: I[i, j] * 1e3 for j, b in enumerate(BANDS + ['total'])}

    def band_densities(self, s, n, energy=None):
        """Main-run band densities at E_F(n) (or at a given energy), 1e11 cm^-2 (holes negative)."""
        ef = self.kp(s, 'main', n)['ef'] if energy is None else energy
        return self.idos_at(s, 'main', n, ef)

    def pocket(self, s, n, wedge, own_fermi=False):
        """Single-pocket hole density (positive) from a wedge run: band -1 IDOS / (360 / width)."""
        width = 60 if wedge == 'w105' else 120
        ef = self.kp(s, wedge, n)['ef'] if own_fermi else self.kp(s, 'main', n)['ef']
        return -self.idos_at(s, wedge, n, ef)[-1] / (360 / width)

    def rigid(self, s, n_ref, n_target):
        """Rigid-band filling of the n_ref run: band densities at the energy where its total IDOS = n_target."""
        k = self.key(s, 'main', n_ref)
        E, I = self.P[f'E_{k}'], self.P[f'I_{k}']
        tot = I[:, -1] * 1e3
        i = int(np.abs(tot - n_target).argmin())
        return {b: I[i, j] * 1e3 for j, b in enumerate(BANDS + ['total'])}, E[i]


def interp(x, xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    o = np.argsort(xs)
    return float(np.interp(x, xs[o], ys[o]))
