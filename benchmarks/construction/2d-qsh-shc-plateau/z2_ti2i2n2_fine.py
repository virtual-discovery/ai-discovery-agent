import sys, json
import numpy as np
sys.path.insert(0, '/repo/benchmarks/construction/2d-qsh-shc-plateau')
import qsh_tb, build_models, reference  # run with run_in_image.sh; output in z2_ti2i2n2_fine.txt
from pathlib import Path
f = Path('/data/x/Data/Wannier_Hamiltonians/insulators/Ti2I2N2_Pmm2')
m = qsh_tb.Model('/repo/benchmarks/tasks/2d-qsh-shc-plateau/environment/lab/data/Ti2I2N2_Pmm2.npz', build_models.cell_of(f))
nocc = 32
for axis in (1, 0):
    for nk in (151, 201, 401, 801):
        k = np.linspace(0, 0.5, nk)
        w = m.wcc(nocc, k, 120, axis)
        z, g = reference.z2_count(w)
        # largest jump of the reference line
        print('axis', axis, 'nk', nk, 'z2', z, 'mingap', round(g, 3), flush=True)
    # detailed look near the gap minimum
    k = np.linspace(0.0, 0.08, 161) if axis == 1 else np.linspace(0, 0.08, 161)
    w = m.wcc(nocc, k, 160, axis)
    def refline(a):
        a = np.sort(a); gg = np.diff(np.r_[a, a[0] + 1]); j = np.argmax(gg); return (a[j] + gg[j] / 2) % 1
    zl = [refline(x) for x in w]
    # print WCC close to 0 and 0.5 (the moving ones)
    for i in range(0, 161, 8):
        near = np.sort(w[i])
        print(round(k[i], 4), 'ref %.3f' % zl[i], np.round(near[(near < 0.08) | (near > 0.92) | ((near > 0.42) & (near < 0.58))], 3).tolist())
