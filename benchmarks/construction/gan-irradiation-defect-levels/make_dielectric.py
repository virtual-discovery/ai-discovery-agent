"""Write dielectric.json from the construction-time DFPT run (dfpt/out/dynmat.out, dfpt/out/ph.out).

    python3 -B make_dielectric.py

eps_inf is the clamped-ion (electronic) tensor; eps_0 adds the zone-centre polar-mode contribution
(dynmat.x with the crystal acoustic sum rule). Both are rounded to 3 decimals and symmetrised to the
hexagonal form diag(perp, perp, par) of the wurtzite structure (z = c axis).
"""
import json
import re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE/'dfpt/out'


def tensor_after(text, marker):
    i = text.index(marker)
    rows = []
    for line in text[i:].splitlines()[1:]:
        nums = re.findall(r'-?\d+\.\d+', line)
        if len(nums) == 3:
            rows.append([float(x) for x in nums])
            if len(rows) == 3:
                return np.array(rows)
    raise ValueError(marker)


def hexagonal(t):
    t = np.array(t)
    perp = (t[0, 0] + t[1, 1]) / 2
    assert np.abs(t - np.diag(np.diag(t))).max() < 1e-3 and abs(t[0, 0] - t[1, 1]) < 1e-3
    return np.diag([perp, perp, t[2, 2]]).round(3).tolist()


def main():
    dyn = (OUT/'dynmat.out').read_text()
    ph = (OUT/'ph.out').read_text()
    eps_inf = tensor_after(dyn, 'Electronic dielectric permittivity tensor')
    eps_0 = tensor_after(dyn, 'with zone-center polar mode contributions')
    eps_ph = tensor_after(ph, 'Dielectric constant in cartesian axis')
    assert np.abs(eps_inf - eps_ph).max() < 1e-4
    table = dyn[dyn.index('# mode'):].splitlines()[1:13]
    freqs = [round(float(line.split()[1]), 2) for line in table]
    out = {'eps_inf': hexagonal(eps_inf), 'eps_0': hexagonal(eps_0),
           'gamma_frequencies_cm-1': freqs,
           'method': 'DFPT (Quantum ESPRESSO 6.7, PBE, pslibrary 1.0.0 ultrasoft potentials with Ga 3d in valence, '
                     '45/360 Ry, 8x8x5 k-mesh) at the lattice parameters (a = 3.1948 A, c = 5.1886 A) and internal '
                     'parameter (u = 0.3775) of the archive\'s perfect supercell; eps_0 includes the zone-centre '
                     'polar phonons (dynmat.x, crystal acoustic sum rule, lperm). Computed by the task builders, not by the authors.'}
    (HERE/'dielectric.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
