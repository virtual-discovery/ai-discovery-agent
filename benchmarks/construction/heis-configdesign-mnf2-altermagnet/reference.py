"""Reference exchange model for heis-configdesign-mnf2-altermagnet: writes j_reference.json.

    python3 -B reference.py

10-shell sublattice-symmetrised Heisenberg fit (with intercept) to every released spiral whose plane-wave basis
represents the spiral (VASP's 'choose ENMAX > X' with X <= ENCUT), i.e. all 101 spirals except [214] t >= 0.90.
This is the model the classical Monte Carlo reference Neel temperature (tn_mc.py) is computed for.
"""
import json
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
NSHELL = 10


def clean_mask(spirals):
    return np.array([s['enmax_needed_eV'] <= s['encut_eV'] for s in spirals])


def main():
    P = M.POOL['spirals']
    Q = np.array([s['q_reciprocal_lattice'] for s in P])
    Y = np.array([s['energy_meV_per_cell'] for s in P])
    ok = clean_mask(P)
    a, J = M.fit(2 * M.X(Q[ok], NSHELL), Y[ok], intercept=True)
    rmse = float(np.sqrt(np.mean((a + 2 * M.X(Q[ok], NSHELL) @ J - Y[ok]) ** 2)))
    shells = [{'distance_A': s['distance'], 'kind': s['kind'], 'count': s['count'], 'altermagnetic_split': M.altermagnetic(s)}
              for s in M.SHELLS[:NSHELL]]
    out = {'J_meV': J.tolist(), 'intercept_meV_per_cell': a, 'n_points': int(ok.sum()), 'rmse_meV_per_cell': rmse,
           'shells': shells, 'mean_field_TN_K': M.neel_mean_field_K(J)}
    (HERE/'j_reference.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
