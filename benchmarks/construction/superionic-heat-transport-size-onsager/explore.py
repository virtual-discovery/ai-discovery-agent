"""Exploration: kappa of every run with several estimators (run in the runtime image)."""
import glob
import json
import sys
import numpy as np
sys.path.insert(0, '/c')
import cepstral as C

RAW = '/raw'
DT = 0.008
out = []
for f in sorted(glob.glob(f'{RAW}/' + (sys.argv[1] if len(sys.argv) > 1 else '**') + '/*_log.npy', recursive=True)):
    a = np.load(f).astype(float)
    parts = f.split('/')
    mat, ens, cell, T = parts[-5], parts[-4], parts[-3], int(parts[-2].split('_')[1])
    V = a[0, 5] ** 3
    Tm = a[:, 3].mean()
    J, Vo = a[:, 6:9], a[:, 9:12]
    pref = C.kappa_prefactor(1.0, Tm) / V / 2           # S0 of extensive J -> kappa
    fr, S, nl = C.spectra([J, Vo], DT)
    res = {'mat': mat, 'ens': ens, 'cell': cell, 'T': T, 'Tmean': round(Tm, 1), 'n': len(a),
           'drift_eV_per_ns': float(np.polyfit(a[:, 0], a[:, 1], 1)[0] * 1000)}
    for fcut in (5.0, 10.0, 20.0):
        k = int(fcut / fr[1])
        s0, e0, p0, _ = C.cepstral(S[:, 0, 0], 3, 1, f_cut_index=k)
        s1, e1, p1, _ = C.cepstral(C.reduced(S), 3, 2, f_cut_index=k)
        res[f'plain_{fcut:g}'] = (round(s0 * pref, 3), round(e0 * pref, 3), p0)
        res[f'mc_{fcut:g}'] = (round(s1 * pref, 3), round(e1 * pref, 3), p1)
    # O/F current zero-frequency value (ionic conductivity proxy)
    s, e, p, _ = C.cepstral(S[:, 1, 1], 3, 1, f_cut_index=int(10 / fr[1]))
    res['SVV0'] = round(s, 3)
    t, I = C.running_integral(J, DT, 20.0)
    res['direct_2ps'] = round(I[int(2 / DT)] * pref * 2, 3)
    res['direct_10ps'] = round(I[int(10 / DT) - 1] * pref * 2, 3)
    out.append(res)
    print(json.dumps(res), flush=True)
json.dump(out, open('/s/explore.json', 'w'), indent=0)
