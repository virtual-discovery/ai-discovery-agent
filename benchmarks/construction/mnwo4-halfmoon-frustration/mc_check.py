"""Monte Carlo cross-check of the archived maps (construction only; ~20 min on 4 CPUs in the runtime image).

    python3 -B mc_check.py CACHE WORKDIR

Runs equilibrium Metropolis + over-relaxation Monte Carlo (mc.c) of the archived Hamiltonian on the archived
36x24x36 box at 14, 20 and 30 K, evaluates S(q) at the centres of 7 x 8 blocks of the Fig.-3 map grid, and prints
the block Pearson correlation of Monte Carlo and of SCGA with the archived maps (SOURCE.md quotes the result).
"""
import subprocess
import sys
from pathlib import Path
import numpy as np
from model import Model, truth, KB, NAX, norm3

HERE = Path(__file__).resolve().parent
CACHE, WORK = Path(sys.argv[1]), Path(sys.argv[2])
WORK.mkdir(parents=True, exist_ok=True)


def blocks7x8(g):
    nk2, nl2 = (g.shape[0] // 7) * 7, (g.shape[1] // 8) * 8
    a, b = (g.shape[0] - nk2) // 2, (g.shape[1] - nl2) // 2
    return g[a:a + nk2, b:b + nl2].reshape(nk2 // 7, 7, nl2 // 8, 8).mean(axis=(1, 3))


def pearson(a, b):
    a = a.ravel() - a.mean(); b = b.ravel() - b.mean()
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def main():
    X = np.load(CACHE/'extracted.npz')
    shells, J, D, _ = truth(CACHE/'model.json')
    M = Model(shells)
    with open(WORK/'bonds.txt', 'w') as f:
        for k, bonds in enumerate(shells):
            for s, s2, n1, n2, n3 in bonds:
                f.write('%d %d %d %d %d %.10f\n' % (s, s2, n1, n2, n3, J[k]))
    K, L = np.meshgrid(X['map_k'], X['map_l'], indexing='ij')
    qc = np.c_[blocks7x8(-L / 2).ravel(), blocks7x8(K).ravel(), blocks7x8(L).ravel()]
    np.savetxt(WORK/'q.txt', qc, fmt='%.8f')
    exe = WORK/'mc'
    subprocess.run(['gcc', '-O3', '-march=x86-64', '-o', str(exe), str(HERE/'mc.c'), '-lm'], check=True)
    procs = {}
    for T in (14, 20, 30):
        procs[T] = subprocess.Popen([str(exe), str(WORK/'bonds.txt'), '36', '24', '36', '%.8f' % (T * KB), '%.10f' % D,
                                     *map(str, NAX), '4000', '12000', '60', str(T), str(WORK/'q.txt'),
                                     str(WORK/('mc_T%d.txt' % T)), '-'], stdout=subprocess.PIPE, text=True)
    for T, p in procs.items():
        p.wait()
        ref = blocks7x8(X['map_T%d' % T][:, :, 3])
        mc = norm3(np.loadtxt(WORK/('mc_T%d.txt' % T)))
        sc = norm3(M.scga(J, D, T, qc, nbz=32))
        print('T = %d K: archived vs MC r = %.3f, archived vs SCGA r = %.3f, MC vs SCGA r = %.3f'
              % (T, pearson(ref, mc), pearson(ref, sc), pearson(mc, sc)))


if __name__ == '__main__':
    main()
