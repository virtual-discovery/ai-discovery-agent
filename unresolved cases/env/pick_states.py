import json, sys, numpy as np
d = sys.argv[1]; e1 = np.loadtxt(d + '/eigenvalues_b1.dat', comments='#'); e2 = np.loadtxt(d + '/eigenvalues_b2.dat', comments='#'); E = e1[:, 0]; f1 = e1[:, 1]; f2 = e2[:, 1]
nsave = int(sys.argv[2]) if len(sys.argv) > 2 else len(E)
def brightest(f, lo, hi):
    m = (E >= lo) & (E <= hi) & (np.arange(len(E)) < nsave); return int(np.argmax(np.where(m, f, -1))) if m.any() else None
al = brightest(f1, 0.9, 1.4); L = brightest(f1, E[al] - 0.15, E[al] - 0.02); ap = brightest(f2, E[al] + 0.05, E[al] + 0.45)
bb = sorted([i for i in np.argsort(-np.where((E >= 1.70) & (E <= 1.95) & (np.arange(len(E)) < nsave), f2, -1))[:2]], key=lambda i: E[i])
dd = sorted(sorted([i for i in range(nsave) if E[i] < E[al] - 0.25], key=lambda i: -(f1[i] + f2[i]))[:4], key=lambda i: E[i])
S = {'dd_a': dd[0] + 1, 'dd_b': dd[1] + 1, 'dd_c': dd[2] + 1, 'dd_d': dd[3] + 1, 'L': (L + 1) if L is not None else al + 1, 'alpha': al + 1, 'alpha_prime': ap + 1, 'beta_a': int(bb[0]) + 1, 'beta_b': int(bb[1]) + 1}
lab = {k: '%s (%.3f)' % (k, E[v - 1]) for k, v in S.items()}
json.dump({'STATES': S, 'LABELS': lab}, open(d + '/states.json', 'w'), indent=1); print(S, {k: round(float(E[v - 1]), 3) for k, v in S.items()})
