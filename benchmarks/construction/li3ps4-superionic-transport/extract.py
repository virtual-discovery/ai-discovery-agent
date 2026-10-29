"""Extract the Li3PS4 data pool from the Materials Cloud archive pkn2d-h8190 (see fetch.sh).

    bash benchmarks/construction/li3ps4-superionic-transport/fetch.sh
    python3 benchmarks/construction/li3ps4-superionic-transport/extract.py [DEST]

Reads only author-archived arrays (raw-data-4-figures and the ideal input cells of the ML-MD runs) and writes
pool.json next to this script. No value is modified; time series are decimated as documented below.
"""
import json
import sys
from pathlib import Path
import numpy as np

DEST = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/li3ps4-superionic-transport')
SRC = DEST/'x/MaterialsCloudArchive'
FIG = SRC/'raw-data-4-figures'
OUT = Path(__file__).with_name('pool.json')
FUNCTIONALS = ['PBEsol', 'r2SCAN', 'PBE0']
PHASES = ['alpha', 'beta', 'gamma']
GRID = [200, 250, 300, 350, 400] + list(range(450, 801, 25)) + list(range(850, 1201, 50))   # simulated temperatures


def nominal(t):
    """Nominal grid temperature of a run (the archive stores the run's mean temperature, e.g. 849.95)."""
    g = min(GRID, key=lambda x: abs(x - t))
    assert abs(g - t) < 0.5, t
    return g


def r(x, n=6):
    return float('%.*g' % (n, x))


def conductivity():
    out = {}
    for f in FUNCTIONALS:
        out[f] = {}
        for p in PHASES:
            series = {}
            for kind, key in (('NE', 'nernst_einstein'), ('GK', 'green_kubo')):
                a = np.loadtxt(FIG/f'Fig8-elec-conductivity/{kind}cond_{p}_{f}.txt')
                assert a.shape == (28, 3)
                series[key] = {str(nominal(t)): {'sigma_S_per_m': r(s), 'block_error_S_per_m': r(e), 'mean_temperature_K': r(t, 7)}
                               for t, s, e in a}
            assert list(series['nernst_einstein']) == [str(g) for g in GRID]
            out[f][p] = series
    return out


def fig7(name, f):
    a = np.load(FIG/f'Fig7-transitions-Cp/{name}{f}.npy')
    return [[r(t, 7), r(v), r(e)] for t, v, e in a.T]


def main():
    pool = {'source': 'Materials Cloud Archive pkn2d-h8190, MaterialsCloudArchive-ChemMat2024.zip (md5 e838e8716fb42200c34f06119337e909)',
            'temperature_grid_K': GRID, 'conductivity': conductivity()}
    # Fig. 7: heat capacity (finite difference between consecutive simulated temperatures, 27 midpoints),
    # PS4 polar-angle rms fluctuation (deg), PS4 rotational diffusion (rad^2/ps), P diffusion (A^2/ps).
    pool['heat_capacity'] = {f: fig7('Cp-', f) for f in FUNCTIONALS}
    pool['angle_rms_deg'] = {f: fig7('angles_', f) for f in FUNCTIONALS}
    pool['rotational_diffusion_rad2_per_ps'] = {f: fig7('Domega_', f) for f in FUNCTIONALS}
    pool['p_diffusion_A2_per_ps'] = {f: fig7('DP_', f) for f in FUNCTIONALS}
    for f in FUNCTIONALS:  # the 25-point Cp_ arrays equal the first 25 points of the 27-point Cp- arrays
        short = np.load(FIG/f'Fig7-transitions-Cp/Cp_{f}.npy')
        assert np.allclose(short, np.load(FIG/f'Fig7-transitions-Cp/Cp-{f}.npy')[:, :25])
    # Fig. S12: Li MSD of the ML-PBE0 run started in gamma at 750 K (A^2, every 10 fs, 654 ps). Decimated:
    # every 10 fs to 2 ps, every 0.1 ps to 20 ps, every 1 ps to the end.
    m = np.loadtxt(FIG/'FigS12-MSD-750K/MSD-Li-750K.txt')
    t = np.arange(len(m)) * 0.01
    keep = np.r_[np.arange(0, 201), np.arange(210, 2001, 10), np.arange(2100, len(m), 100), len(m) - 1]
    keep = np.unique(keep)
    pool['li_msd_750K'] = {'t_ps': [r(x, 6) for x in t[keep]], 'msd_A2': [r(x, 6) for x in m[keep, 0]],
                           'block_error_A2': [r(x, 4) for x in m[keep, 1]], 'n_raw_points': int(len(m))}
    # Fig. S13: running GK integral sigma(tau) (S/cm) of the total conductivity and of the Li-S term,
    # in the matrix (P+S centre of mass fixed) and barycentric frames; every second point (20 fs).
    tot = np.loadtxt(FIG/'FigS13-sigma-breakdown/sigmatotal.dat')
    lis = np.loadtxt(FIG/'FigS13-sigma-breakdown/sigmaLiS.dat')
    assert np.allclose(tot[:, 0], lis[:, 0])
    sel = np.arange(0, len(tot), 2)
    pool['current_breakdown'] = {'tau_ps': [r(x, 5) for x in tot[sel, 0]],
                                 'matrix': {'total_S_per_cm': [r(x, 6) for x in tot[sel, 1]], 'Li_S_S_per_cm': [r(x, 6) for x in lis[sel, 1]]},
                                 'barycentric': {'total_S_per_cm': [r(x, 6) for x in tot[sel, 2]], 'Li_S_S_per_cm': [r(x, 6) for x in lis[sel, 2]]}}
    # Fig. 9: integrated Li-A velocity cross-correlation vs distance and Li-A radial distribution functions
    pc = {}
    for pair in ('LiLi', 'LiP', 'LiS'):
        a = np.loadtxt(FIG/f'Fig9-SpatialCorrelation/IntCrossCorr-{pair}.data')
        g = np.loadtxt(FIG/f'Fig9-SpatialCorrelation/RDF-{pair}.data')
        pc[pair] = {'r_A': [r(x, 5) for x in a[:, 0]], 'I_A2_per_ps2': [r(x, 5) for x in a[:, 1]], 'I_error': [r(x, 4) for x in a[:, 2]],
                    'rdf_r_A': [r(x, 5) for x in g[:, 0]], 'rdf': [r(x, 5) for x in g[:, 1]], 'rdf_error': [r(x, 4) for x in g[:, 2]]}
    pool['li_pair_correlation'] = pc
    # Fig. 4: collective variables (s1, s2) of ML-PBE0 runs started in gamma at 700/750/800 K, and reference
    # clusters from alpha/beta/gamma runs below the transition. Every 7th sample of the varying-T sets.
    z = np.load(FIG/'Fig4-CVs/CVsgammavaryingT.npz')
    temps = [int(x) for x in np.load(FIG/'Fig4-CVs/Temperatures.npy')]
    cvs = {str(T): np.round(z[f'CVs{T}K'][::7], 4).tolist() for T in temps}
    ref = {}
    for p in PHASES:
        a = np.load(FIG/f'Fig4-CVs/CVs{p}.npy')
        ref[p] = {'mean': np.round(a.mean(0), 4).tolist(), 'std': np.round(a.std(0), 4).tolist(), 'n_samples': int(len(a))}
    pool['order_parameters'] = {'samples': cvs, 'reference_clusters': ref}
    # Ideal starting cells (768 atoms) of the ML-PBE0 runs
    cells = {}
    gin = (SRC/'MD-simulations/ML-PBE0/gamma-768atoms/input-structures/gamma_LPS_4x4x3.input').read_text().split('\n')
    lo_hi = [float(line.split()[1]) - float(line.split()[0]) for line in gin if line.strip().endswith(('xlo xhi', 'ylo yhi', 'zlo zhi'))]
    assert len(lo_hi) == 3
    cells['gamma'] = {'edges_A': lo_hi, 'volume_A3': r(np.prod(lo_hi))}
    bline = (SRC/'MD-simulations/ML-PBE0/beta-768atoms/input-structures/beta_align_732.xyz').read_text().split('\n')[1]
    abc = [float(x) for x in bline.split(':')[1].split()[:3]]
    cells['beta'] = {'edges_A': abc, 'volume_A3': r(np.prod(abc))}
    aline = (SRC/'MD-simulations/ML-PBE0/alpha-768atoms/input-structures/alpha_align_732.extxyz').read_text().split('\n')[1]
    lat = np.array([float(x) for x in aline.split('"')[1].split()]).reshape(3, 3)
    cells['alpha'] = {'lattice_A': lat.round(6).tolist(), 'volume_A3': r(abs(np.linalg.det(lat)))}
    pool['ideal_cells'] = cells
    OUT.write_text(json.dumps(pool, separators=(',', ':')) + '\n')
    print('wrote', OUT, OUT.stat().st_size, 'bytes; cells', {k: v['volume_A3'] for k, v in cells.items()})


if __name__ == '__main__':
    main()
