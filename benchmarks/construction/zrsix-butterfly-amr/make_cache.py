"""Chambers tables for the lab: per compound and field direction, orbit-averaged contributions of every start point
of the 160 x 160 x 80 Gamma-centred mesh whose energy lies in the compound's window, binned (cloud-in-cell, 1 meV)
by start energy and band, for the sub-meshes 40, 80 and 160 (the 40 and 80 meshes are subsets of the 160 mesh).

    NPROC=4 bash run_in_image.sh make_cache.py ZrSiSe [plane_angle ...]

Writes /data/cache/<name>/<dkey>.npz with
  K[mesh, band, bin, X, 3, 3]  sum over start points of <v_i vbar_j>  ((eV A)^2), per sub-mesh
  W[mesh, band, bin]           number of start points (CIC weights)
  O[mesh, band, bin]           number of start points whose orbit is open (not closed, or closed with G != 0)
Thermal weights -df/dE are applied at lookup from the bin energies.
"""
import numpy as np, sys, os, time, multiprocessing as mp
from pathlib import Path
sys.path.insert(0, '.')
import spline, chambers as C
from common import GROUP, WINDOW, EBIN, XGRID, MESH_TOP, SUBMESH, directions, bvec, dkey

NAME = sys.argv[1]
TOP = int(os.environ.get('ZRSIX_TOP', MESH_TOP))
SUBS = SUBMESH if TOP == MESH_TOP else [TOP]
OUT = Path(f'/data/cache/{NAME}' if TOP == MESH_TOP else f'/data/cache{TOP}/{NAME}'); OUT.mkdir(parents=True, exist_ok=True)
z = np.load(f'/data/spline/{NAME}.npz')
SP = spline.GroupSpline(z['coef'], z['B'], float(z['eref']))
LO, HI = WINDOW[NAME]
NB = int(round((HI - LO) / EBIN)) + 1
N = TOP


def starts():
    g = np.arange(N); g3 = np.arange(N // 2)
    I = np.stack(np.meshgrid(g, g, g3, indexing='ij'), -1).reshape(-1, 3)
    K = I / np.array([N, N, N // 2])
    e = np.concatenate([SP.eig_vel(K[s:s + 200000])[0] for s in range(0, len(K), 200000)])
    ks, bs = np.nonzero((e > LO) & (e < HI))
    member = np.stack([np.all(I[ks] % (N // m) == 0, axis=1) for m in SUBS], 0)     # (nmesh, M)
    return K[ks], bs, e[ks, bs], member


def work(args):
    k0, band, bv = args
    tr = C.integrate(SP, k0, band, bv, 10 * XGRID.max())
    avg, _ = C.orbit_average(tr, XGRID)
    openf = (~tr['closed']) | (np.abs(tr['G']).sum(1) > 0)
    return avg.astype(np.float64), openf, tr['energy']


def main():
    t0 = time.time()
    K, band, e, member = starts()
    print(NAME, 'start points', len(K), 'per mesh', member.sum(1), round(time.time() - t0), 's', flush=True)
    todo = [d for d in directions(NAME) if len(sys.argv) <= 2 or dkey(*d) in sys.argv[2:]]
    print('directions', len(todo), flush=True)
    nproc = int(os.environ.get('ZRSIX_NPROC', '4'))
    chunk = 600
    with mp.Pool(nproc) as pool:
        for d in todo:
            path = OUT/f'{dkey(*d)}.npz'
            if path.exists():
                continue
            t = time.time()
            bv = bvec(*d)
            jobs = [(K[s:s + chunk], band[s:s + chunk], bv) for s in range(0, len(K), chunk)]
            res = pool.map(work, jobs, chunksize=1)
            avg = np.concatenate([r[0] for r in res]); openf = np.concatenate([r[1] for r in res])
            e0 = np.concatenate([r[2] for r in res])
            x = (e0 - LO) / EBIN
            i0 = np.clip(np.floor(x).astype(int), 0, NB - 2); f = np.clip(x - i0, 0, 1)
            nband = len(GROUP[NAME])
            Ksum = np.zeros((len(SUBS), nband, NB, len(XGRID), 3, 3))
            W = np.zeros((len(SUBS), nband, NB)); O = np.zeros_like(W)
            for mi in range(len(SUBS)):
                sel = member[mi]
                for ib, wt in ((i0, 1 - f), (i0 + 1, f)):
                    np.add.at(Ksum[mi], (band[sel], ib[sel]), wt[sel, None, None, None] * avg[sel])
                    np.add.at(W[mi], (band[sel], ib[sel]), wt[sel])
                    np.add.at(O[mi], (band[sel], ib[sel]), wt[sel] * openf[sel])
            np.savez_compressed(path, K=Ksum.astype(np.float32), W=W.astype(np.float32), O=O.astype(np.float32),
                                lo=LO, ebin=EBIN, X=XGRID, mesh=np.array(SUBS), bvec=bv)
            print(dkey(*d), 'orbits', len(K), 'open frac', openf.mean().round(4), round(time.time() - t), 's', flush=True)


if __name__ == '__main__':
    main()
