"""Sample LAMMPS frames out of the nested traj/trajectory.tar members by range reads, and analyse them.

A run's trajectory.tar stores one member '<step>.lammpstrj' per dump (every 100 steps = 20 fs; columns
id type xu yu zu fx fy fz, unwrapped coordinates in Angstrom). Members are read in contiguous chunks at
spread positions of the nested tar, so each chunk yields short bursts of consecutive frames.
"""
import re
import numpy as np
import remotetar as RT

MAGIC = b'ustar'


def frames_from_chunk(blob, base_off, chunk_off):
    """Parse every complete member inside blob (blob starts at absolute offset chunk_off; members 512-aligned to base_off)."""
    out = []
    start = (-(chunk_off - base_off)) % 512
    p = start
    while p + 512 <= len(blob):
        h = blob[p:p + 512]
        name = h[:100].rstrip(b'\0')
        if h[257:262] == MAGIC and re.fullmatch(rb'\d+\.lammpstrj', name):
            size = int(h[124:136].rstrip(b'\0 ').strip() or b'0', 8)
            if p + 512 + size > len(blob):
                break
            out.append(parse_frame(blob[p + 512:p + 512 + size].decode()))
            p += 512 + ((size + 511) // 512) * 512
        else:
            p += 512
    return out


def parse_frame(text):
    lines = text.split('\n')
    step = int(lines[1])
    n = int(lines[3])
    lo, hi = [float(x) for x in lines[5].split()[:2]]
    L = hi - lo
    data = np.array([[float(x) for x in l.split()[:5]] for l in lines[9:9 + n]])
    order = np.argsort(data[:, 0])
    return {'step': step, 'L': L, 'xyz': data[order, 2:5]}


def sample_frames(tar, entry, n_chunks, chunk_bytes, skip_frac=0.0):
    base, size = entry['offset'], entry['size']
    frames = {}
    for k in range(n_chunks):
        rel = skip_frac + (1 - skip_frac) * (k + 0.5) / n_chunks
        off = base + int(rel * (size - chunk_bytes))
        off -= (off - base) % 512
        blob = RT.fetch(tar, off, min(off + chunk_bytes, base + size) - 1)
        for f in frames_from_chunk(blob, base, off):
            frames[f['step']] = f
    return [frames[s] for s in sorted(frames)]


def kvectors(L, kmax):
    nmax = int(np.floor(kmax * L / (2 * np.pi)))
    r = np.arange(-nmax, nmax + 1)
    n = np.array(np.meshgrid(r, r, r, indexing='ij')).reshape(3, -1).T
    # half space (k and -k give the same S)
    keep = (n[:, 0] > 0) | ((n[:, 0] == 0) & (n[:, 1] > 0)) | ((n[:, 0] == 0) & (n[:, 1] == 0) & (n[:, 2] > 0))
    n = n[keep]
    k = 2 * np.pi * n / L
    kn = np.linalg.norm(k, axis=1)
    return n[kn <= kmax], kn[kn <= kmax]


def max_sk_over_n(frame, kmax=4.0):
    L, x = frame['L'], frame['xyz']
    N = len(x)
    n, kn = kvectors(L, kmax)
    nmax = np.abs(n).max()
    ph = np.exp(2j * np.pi * x[:, :, None] * np.arange(-nmax, nmax + 1)[None, None, :] / L)   # (N,3,2nmax+1)
    best, kbest = 0.0, 0.0
    for c in range(0, len(n), 2000):
        nn = n[c:c + 2000] + nmax
        rho = (ph[:, 0, nn[:, 0]] * ph[:, 1, nn[:, 1]] * ph[:, 2, nn[:, 2]]).sum(axis=0)
        s = np.abs(rho) ** 2 / N
        i = int(np.argmax(s))
        if s[i] > best:
            best, kbest = float(s[i]), float(kn[c + i])
    return best / N, kbest


def pair_distances(frame, rmax):
    L, x = frame['L'], frame['xyz'] % frame['L']
    d = x[:, None, :] - x[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d ** 2).sum(-1))
    iu = np.triu_indices(len(x), 1)
    r = r[iu]
    return r, iu


def gofr(frames, rmax, dr=0.02):
    edges = np.arange(0, rmax + dr / 2, dr)
    hist = np.zeros(len(edges) - 1)
    for f in frames:
        r, _ = pair_distances(f, rmax)
        hist += np.histogram(r[r < rmax], edges)[0]
    N, L = len(frames[0]['xyz']), frames[0]['L']
    rho = N / L ** 3
    shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / (len(frames) * N / 2 * rho * shell)
    return 0.5 * (edges[1:] + edges[:-1]), g


def bursts(frames, min_len=5):
    """Groups of consecutive frames (step difference 100) at least min_len long (>= 80 fs)."""
    out, cur = [], [frames[0]]
    for f in frames[1:]:
        if f['step'] - cur[-1]['step'] == 100:
            cur.append(f)
        else:
            if len(cur) >= min_len:
                out.append(cur)
            cur = [f]
    if len(cur) >= min_len:
        out.append(cur)
    return out


def stable_molecular_fraction(frames, rc=1.05, min_len=5):
    vals = []
    for b in bursts(frames, min_len):
        for s in range(0, len(b) - min_len + 1, min_len):
            keep = None
            for f in b[s:s + min_len]:
                r, iu = pair_distances(f, rc)
                close = set(zip(iu[0][r < rc].tolist(), iu[1][r < rc].tolist()))
                keep = close if keep is None else keep & close
            vals.append(2 * len(keep) / len(b[0]['xyz']))
    return float(np.mean(vals)) if vals else None, len(vals)


def msd_table(frames, dt_ps=0.02 * 1e-3 * 1e3, lags_ps=(1, 2, 5, 10, 20, 40)):
    t = np.array([f['step'] * 2e-4 for f in frames])
    x = np.array([f['xyz'] for f in frames])
    out = {}
    for lag in lags_ps:
        vals = []
        for i in range(len(t)):
            j = np.argmin(np.abs(t - (t[i] + lag)))
            if abs(t[j] - t[i] - lag) <= 0.1 * lag:
                vals.append(((x[j] - x[i]) ** 2).sum(-1).mean())
        if vals:
            out[lag] = float(np.mean(vals))
    return out
