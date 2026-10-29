import numpy as np, sdh
P = np.load('/scratch/pool.npz'); s='70'
vg, B, D = P[f'tr_{s}_vg'], P[f'tr_{s}_B'], P[f'tr_{s}_data']
n, pw = sdh.chart(vg, B, D, 557, bmax=1.0, pad=4)
L = np.log10(pw)
for blocks in ([-0.45,-0.5,-0.6,-0.7],[-0.5,-0.6,-0.7],[-0.6,-0.7],[-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4]):
    m = np.zeros(len(vg),bool)
    for c in blocks: m |= np.abs(vg-c) <= 0.025+1e-9
    best=None
    for sl in np.arange(2.0,10.01,0.05):
        for v0 in np.arange(-0.45,-0.0999,0.005):
            nn = sl*(v0-vg[m])
            ok = (nn>0.9)&(nn<4.0)
            if ok.sum()<20: continue
            val = np.mean([np.interp(x, n, L[i]) for x,i in zip(nn[ok], np.where(m)[0][ok])])
            if best is None or val>best[0]: best=(val,sl,v0)
    val,sl,v0=best
    Lv=sdh.averaged_log(vg,pw,-0.70,0.02); pk=sdh.peak_near(n,Lv,sl*(v0+0.70),0.15)
    print(blocks, 'slope %.2f v0 %.3f line@-0.70 %.3f peak %.3f'%(sl,v0,sl*(v0+0.7),pk))
