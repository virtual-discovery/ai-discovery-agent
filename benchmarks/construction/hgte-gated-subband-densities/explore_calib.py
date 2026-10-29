"""Exploration: truth ridges and their spread under defensible SdH analysis variants (prints only)."""
import json, sys
import numpy as np
import sdh
P = np.load('/scratch/pool.npz')
RPRE = {'26': 1000, '70': 557, '107': 10000}
s = '70'
vg, B, D = P[f'tr_{s}_vg'], P[f'tr_{s}_B'], P[f'tr_{s}_data']
# Hall line, authors' convention
nH = np.array([sdh.hall_density(B, D[g, :, 0], D[g, :, 2], RPRE[s]) for g in range(len(vg))])
m = (vg > 0.05) & (vg <= 0.29 + 1e-9)
print('hall line', np.polyfit(vg[m], nH[m], 1), 'authors 9.81968542 3.20604992')
T = {'E1_top_p020': (0.20, 2.81, 0.3), 'E2_bot_p005': (0.05, 0.83, 0.15), 'E3_cb1_p030': (0.30, 1.22, 0.15),
     'E4_top_m045': (-0.45, 0.45, 0.1), 'E5_bot_m045': (-0.45, 0.63, 0.1), 'E6_hole_m070': (-0.70, 2.49, 0.2)}
variants = {
    'authors_avg0.02': dict(kw=dict(), hw=0.02, q='sxx'),
    'single': dict(kw=dict(), hw=0.0, q='sxx'),
    'avg0.01': dict(kw=dict(), hw=0.01, q='sxx'),
    'block0.024': dict(kw=dict(), hw=0.024, q='sxx'),
    'bmax1.0_pad4': dict(kw=dict(bmax=1.0, pad=4), hw=0.02, q='sxx'),
    'bmax0.6_pad4': dict(kw=dict(bmax=0.6, pad=4), hw=0.02, q='sxx'),
    'rxx_hann': dict(kw=dict(window='hann'), hw=0.02, q='rxx'),
    'bmax1.5_bmin0.2_pad4': dict(kw=dict(bmax=1.5, bmin=0.2, pad=4), hw=0.02, q='sxx'),
    'bmax2.0_bmin0.15_pad4': dict(kw=dict(bmax=2.0, bmin=0.15, pad=4), hw=0.02, q='sxx'),
}
for name, v in variants.items():
    n, pw = sdh.chart(vg, B, D, RPRE[s], quantity=v['q'], **v['kw'])
    row = []
    for t, (V, g, w) in T.items():
        L = sdh.averaged_log(vg, pw, V, v['hw'])
        pk = sdh.peak_near(n, L, g, w)
        row.append('%s=%s' % (t, 'None' if pk is None else '%.3f' % pk))
    print('%-22s' % name, ' '.join(row))
