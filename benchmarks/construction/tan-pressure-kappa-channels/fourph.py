"""Four-phonon scattering rates on the 16x16x16 grid by Monte Carlo sampling (construction only).

    python3 -B fourph.py PRESSURE MESH SCALEBROAD NSAMPLE SEED T_K
Run inside ai-theorist-benchmark-runtime:20261007 with the archive mounted at /raw and the output folder at /s
(see README). One run (8192 sampled (q', q'') pairs per irreducible q, all branch combinations) takes ~25-40
min on one core. Output: /s/r4_<P>_<mesh>_<sb>_<nsample>_<seed>_<T>.npz with arrays pp, pm, mm (ps^-1).
"""
import sys, time, numpy as np
sys.path.insert(0, '/repo')
import tanbte as tb
R = '/raw/'          # mount DEST/raw/IFCs-data-pressure-TaN here
P = sys.argv[1]; m = int(sys.argv[2]); sb = float(sys.argv[3]); ns = int(sys.argv[4]); seed = int(sys.argv[5])
T = float(sys.argv[6])
pc = tb.read_poscar(R + P + '/POSCAR')
fc2 = tb.read_fc2(R + P + '/FORCE_CONSTANTS_2ND', 2, (5, 5, 5))
fc4 = tb.read_fcn(R + P + '/FORCE_CONSTANTS_4TH', 4)
cr = tb.Crystal(pc, fc2, fc4=fc4)
g = tb.Grid(cr, (m, m, m))
t = time.time()
out, hits = tb.four_phonon_mc(g, sb, T, ns, seed)
np.savez('/s/r4_%s_%d_%g_%d_%d_%g.npz' % (P, m, sb, ns, seed, T), hits=hits, **out)
print('done', time.time() - t, flush=True)
