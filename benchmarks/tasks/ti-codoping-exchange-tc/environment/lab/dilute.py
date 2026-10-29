"""Random substitutional dopants on the Bi sublattice of Bi2Te3 and classical Monte Carlo (lab side).

Lattice: hexagonal a = 4.38296 A, c = 30.48701 A (archive Bi2Te3.cif); Bi1 at z = 0.4 and Bi2 at z = 0.6 of the
rhombohedral cell. Orthorhombic supercell (7L a) x (4L sqrt(3) a) x (L c), 336 L^3 Bi sites, periodic.
Model: H = -sum_{i<j} J_ij e_i.e_j over spin-carrying dopants, J of the shell (archive shell index) applied to
every neighbour of that shell, nothing beyond the 17 archived shells.
"""
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

A_HEX, C_HEX = 4.382963164015663, 30.48701213585697
RADII = np.array([4.38296, 4.78812, 6.09740, 6.49125, 7.50897, 7.59151, 7.83159, 8.76593, 9.73744, 9.98832,
                  10.47266, 10.67812, 10.90826, 11.35283, 11.59623, 11.75506, 12.16952])
HERE = Path(__file__).resolve().parent
BIN = os.environ.get('MC_BIN', str(HERE/'mc_dilute'))
_SITES = {}


def bi_sites(L):
    if L in _SITES:
        return _SITES[L]
    a, c, s3 = A_HEX, C_HEX, np.sqrt(3)
    base = []
    for fx in (0.0, 0.5):
        for k, (dx, dy) in enumerate([(0, 0), (2 / 3, 1 / 3), (1 / 3, 2 / 3)]):
            x, y = a * (dx - dy / 2), a * s3 / 2 * dy
            for zf in (0.4, 0.6):
                base.append((fx * a + x, fx * a * s3 + y, ((zf + k / 3) % 1) * c))
    base = np.array(base)
    nx, ny, nz = 7 * L, 4 * L, L
    box = np.array([nx * a, ny * a * s3, nz * c])
    g = np.array([(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)], float)
    pos = ((g * np.array([a, a * s3, c]))[:, None, :] + base[None]).reshape(-1, 3) % box
    _SITES[L] = (pos, box)
    return _SITES[L]


def sample(L, conc, J, seed, rcut=12.3, tol=0.02, positions=False):
    """conc: {species: concentration per Bi site}; J: {(a, b): 17 values} for every species pair (both orders).
    Returns species index per spin and a CSR neighbour list (offsets, neighbours, couplings)
    [, spin positions and box when positions=True]."""
    pos, box = bi_sites(L)
    rng = np.random.default_rng(seed)
    names = list(conc)
    counts = [int(round(conc[n] * len(pos))) for n in names]
    chosen = rng.choice(len(pos), size=sum(counts), replace=False)
    spec = np.repeat(np.arange(len(names)), counts)
    p = pos[chosen]
    pairs = cKDTree(p, boxsize=box).query_pairs(rcut, output_type='ndarray')
    n = len(p)
    if len(pairs) == 0:
        out = spec, np.zeros(n + 1, int), np.zeros(0, int), np.zeros(0)
        return out + (p, box) if positions else out
    d = p[pairs[:, 1]] - p[pairs[:, 0]]
    d -= box * np.round(d / box)
    r = np.linalg.norm(d, axis=1)
    k = np.abs(r[:, None] - RADII[None, :]).argmin(axis=1)
    assert (np.abs(r - RADII[k]) < tol).all()
    table = np.zeros((len(names), len(names), len(RADII)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            table[i, j] = J[(a, b)]
    v = table[spec[pairs[:, 0]], spec[pairs[:, 1]], k]
    keep = v != 0
    i, j, v = pairs[keep, 0], pairs[keep, 1], v[keep]
    rows, cols, vals = np.r_[i, j], np.r_[j, i], np.r_[v, v]
    order = np.lexsort((cols, rows))
    rows, cols, vals = rows[order], cols[order], vals[order]
    offs = np.r_[0, np.cumsum(np.bincount(rows, minlength=n))]
    return (spec, offs, cols, vals, p, box) if positions else (spec, offs, cols, vals)


def write_graph(path, spec, offs, cols, vals, pos=None, box=None):
    blocks = [(np.array([len(spec)]), '<i4'), (spec, '<i4'), (offs, '<i4'), (cols, '<i4'), (vals, '<f8')]
    if pos is not None:
        blocks += [(np.asarray(pos).reshape(-1), '<f8'), (np.asarray(box), '<f8')]
    with open(path, 'wb') as f:
        for arr, dt in blocks:
            np.asarray(arr, dt).tofile(f)


def run_sample(L, conc, J, seed, temps, neq, nmeas, nor=2, workdir=None):
    spec, offs, cols, vals, pos, box = sample(L, conc, J, seed, positions=True)
    if len(conc) == 1:
        spec = np.zeros_like(spec)
    with tempfile.TemporaryDirectory(dir=workdir) as d:
        g = Path(d)/'graph.bin'
        write_graph(g, spec, offs, cols, vals, pos, box)
        out = subprocess.run([BIN, str(g), str(seed), str(neq), str(nmeas), str(nor)] + ['%.6f' % t for t in temps],
                             capture_output=True, text=True, check=True).stdout
    data = np.array([[float(x) for x in line.split()] for line in out.strip().split('\n')])
    return {'n_spins': int(len(spec)), 'n_bonds': int(len(cols) // 2), 'data': data}


def workers():
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def run_many(L, conc, J, seeds, temps, neq, nmeas, nor=2, workdir=None, nworkers=None):
    temps = sorted(temps, reverse=True)
    with ThreadPoolExecutor(nworkers or workers()) as ex:
        return list(ex.map(lambda s: run_sample(L, conc, J, s, temps, neq, nmeas, nor, workdir), seeds)), temps
