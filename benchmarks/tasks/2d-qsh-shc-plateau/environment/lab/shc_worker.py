"""Worker process of the 2d-qsh-shc-plateau lab (copied unchanged into the task lab).

    python3 shc_worker.py WORKDIR PHASE INDEX NPROC

WORKDIR/request.json holds {material, npz, cell_A, gammas, eta, mus, convention, mesh, refine}.
PHASE 'coarse': mesh rows j with j % NPROC == INDEX -> WORKDIR/coarse-INDEX.npy (rows, N1, n_gamma, n_mu).
PHASE 'refine': WORKDIR/candidates.npy (coarse k-points to refine), entries INDEX::NPROC -> WORKDIR/refine-INDEX.npy.
"""
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import qsh_tb  # noqa: E402


def main():
    work, phase, index, nproc = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    req = json.loads((work/'request.json').read_text())
    m = qsh_tb.Model(req['npz'], req['cell_A'])
    tau = m.centres if req['convention'] == 'wannier_centres' else None
    gammas, eta, mus = req['gammas'], req['eta'], np.array(req['mus'], float)
    n1, n2 = req['mesh']
    if phase == 'coarse':
        k1s = np.arange(n1) / n1
        out = []
        for j in range(index, n2, nproc):
            row = []
            for s in range(0, n1, qsh_tb.CHUNK):
                H, dH = m.row(j / n2, k1s[s:s + qsh_tb.CHUNK])
                row.append(m.contributions(H, dH, gammas, eta, mus, tau))
            out.append(np.concatenate(row))
        np.save(work/f'coarse-{index}.npy', np.array(out).reshape(-1, n1, len(gammas), len(mus)))
    else:
        f = req['refine'][0]
        cand = np.load(work/'candidates.npy')[index::nproc]
        off = (np.arange(f) + 0.5) / f - 0.5
        sub = np.array([(a / n1, b / n2) for a in off for b in off])
        per = max(1, qsh_tb.CHUNK // (f * f))
        res = []
        for s in range(0, len(cand), per):
            block = cand[s:s + per]
            kk = (block[:, None, :] + sub[None]).reshape(-1, 2)
            H, dH = m.hk_dhk(np.c_[kk, np.zeros(len(kk))])
            cs = m.contributions(H, dH, gammas, eta, mus, tau).reshape(len(block), f * f, len(gammas), len(mus))
            res.append(cs.mean(1))
        np.save(work/f'refine-{index}.npy', np.concatenate(res) if res else np.zeros((0, len(gammas), len(mus))))


if __name__ == '__main__':
    main()
