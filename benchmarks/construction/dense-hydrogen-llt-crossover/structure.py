"""Structural observables from the sampled frames (run inside the runtime image).

    python3 -B structure.py     # reads <cache>/frames, writes <cache>/pool_structure.json

Per run (all sampled frames, i.e. six chunks of about ten consecutive dumps spread over the run):
  frames:  time (ps), max_k S(k)/N over |k| <= 5.5 1/A (k = 2 pi n / L), the |k| of that maximum, and the
           instantaneous fraction of atoms that have a neighbour closer than 1.05 A;
  stable_molecular_fraction: per window of 5 consecutive dumps (80 fs), 2 x (pairs closer than 1.05 A in all
           five dumps) / N;
  gofr:    g(r) per chunk (0.02 A bins up to min(L/2, 4 A));
  msd:     mean-square displacement (unwrapped coordinates) between sampled frames, binned by lag.
"""
import json
import os
from pathlib import Path
import numpy as np
import traj

CACHE = Path(os.environ.get('HCACHE', os.path.expanduser('~/.cache/ai-theorist/dense-hydrogen-llt-crossover')))
KMAX, RC, DR = 5.5, 1.05, 0.02
LAGS = [1, 2, 5, 10, 20, 50, 100, 200]


def close_pairs(f, rc):
    L, x = f['L'], f['xyz'] % f['L']
    out = set()
    n = len(x)
    for i0 in range(0, n, 512):
        d = x[i0:i0 + 512, None, :] - x[None, :, :]
        d -= L * np.round(d / L)
        r2 = (d ** 2).sum(-1)
        ii, jj = np.nonzero(r2 < rc * rc)
        ii = ii + i0
        keep = jj > ii
        out |= set(zip(ii[keep].tolist(), jj[keep].tolist()))
    return out


def gofr(frames, rmax):
    edges = np.arange(0, rmax + DR / 2, DR)
    hist = np.zeros(len(edges) - 1)
    for f in frames:
        L, x = f['L'], f['xyz'] % f['L']
        for i0 in range(0, len(x), 512):
            d = x[i0:i0 + 512, None, :] - x[None, :, :]
            d -= L * np.round(d / L)
            r = np.sqrt((d ** 2).sum(-1))
            ii = np.arange(i0, min(i0 + 512, len(x)))[:, None]
            jj = np.arange(len(x))[None, :]
            r = r[jj > ii]
            hist += np.histogram(r[r < rmax], edges)[0]
    n, L = len(frames[0]['xyz']), frames[0]['L']
    shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    return (hist / (len(frames) * n / 2 * (n / L ** 3) * shell)).round(4).tolist()


def analyse(rsdir):
    chunks = []
    for k in range(6):
        b = rsdir/f'chunk{k}.bin'
        if not b.exists():
            continue
        m = json.loads((rsdir/f'chunk{k}.json').read_text())
        fr = traj.frames_from_chunk(b.read_bytes(), m['base'], m['offset'])
        if fr:
            chunks.append(sorted(fr, key=lambda f: f['step']))
    allf = {f['step']: f for c in chunks for f in c}
    frames = [allf[s] for s in sorted(allf)]
    n, L = len(frames[0]['xyz']), frames[0]['L']
    rows = []
    for f in frames:
        s, kb = traj.max_sk_over_n(f, KMAX)
        pairs = close_pairs(f, RC)
        bonded = len({i for p in pairs for i in p}) / n
        rows.append([round(f['step'] * 2e-4, 3), round(s, 5), round(kb, 3), round(bonded, 4)])
    stable = []
    for c in chunks:
        for b in traj.bursts(c, 5):
            for s0 in range(0, len(b) - 4, 5):
                keep = None
                for f in b[s0:s0 + 5]:
                    p = close_pairs(f, RC)
                    keep = p if keep is None else keep & p
                stable.append([round(b[s0]['step'] * 2e-4, 3), round(2 * len(keep) / n, 4)])
    rmax = round(min(L / 2, 4.0) - 1e-9, 2)
    g = [{'t_ps_first': round(c[0]['step'] * 2e-4, 2), 't_ps_last': round(c[-1]['step'] * 2e-4, 2),
          'n_frames': len(c), 'g': gofr(c, rmax)} for c in chunks]
    t = np.array([f['step'] * 2e-4 for f in frames])
    x = np.array([f['xyz'] for f in frames])
    msd = []
    for lag in LAGS:
        vals = []
        for i in range(len(t)):
            j = int(np.argmin(np.abs(t - (t[i] + lag))))
            if abs(t[j] - t[i] - lag) <= 0.15 * lag:
                vals.append(float(((x[j] - x[i]) ** 2).sum(-1).mean()))
        if vals:
            msd.append([lag, round(float(np.mean(vals)), 4), len(vals)])
    return {'n_atoms': n, 'box_A': round(L, 6), 'frames': rows, 'stable_molecular_fraction': stable,
            'gofr_r_max_A': rmax, 'gofr_dr_A': DR, 'gofr': g, 'msd': msd}


def main():
    """Optional arguments SHARD NSHARDS: process every NSHARDS-th run into pool_structure.SHARD.json, then
    `python3 -B structure.py merge` folds the shard files into pool_structure.json."""
    import sys
    out_path = CACHE/'pool_structure.json'
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    if sys.argv[1:] == ['merge']:
        for f in sorted(CACHE.glob('pool_structure.*.json')):
            out.update(json.loads(f.read_text()))
            f.unlink()
        out_path.write_text(json.dumps(dict(sorted(out.items()))))
        print(len(out), 'runs')
        return
    shard, nshard = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) == 3 else (0, 1)
    done = set(out)
    if nshard > 1:
        out_path = CACHE/f'pool_structure.{shard}.json'
        out = json.loads(out_path.read_text()) if out_path.exists() else {}
    todo = []
    for tard in sorted((CACHE/'frames').glob('N_*')):
        n, t = tard.name.split('_')[1], tard.name.split('_T')[1].split('K')[0]
        for rsdir in sorted(tard.glob('rs*')):
            key = f'{n}/{t}/{rsdir.name[2:]}'
            if key not in done and key not in out and len(list(rsdir.glob('chunk*.bin'))) == 6:
                todo.append((key, rsdir))
    for key, rsdir in todo[shard::nshard]:
        if True:
            out[key] = analyse(rsdir)
            print(key, len(out[key]['frames']), 'frames', flush=True)
            out_path.write_text(json.dumps(out))
    print(len(out), 'runs')


if __name__ == '__main__':
    main()
