"""Reference solution for heis-configdesign-mnf2-altermagnet; runs as the Harbor oracle and in check_task."""
import json
import os
import re
import subprocess
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 0.08617333262          # meV/K
HB2M = 3.80998212           # hbar^2 / 2 m_e in eV A^2


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


# ---------------------------------------------------------------- structure and shells
CELL = np.array(CATALOG['structure']['cell_A'])
FRAC = np.array(CATALOG['structure']['fractional_positions'])
TAU1, TAU2 = FRAC[0], FRAC[1]            # Mn1 (1/2,1/2,1/2), Mn2 (0,0,0)
BREC = 2 * np.pi * np.linalg.inv(CELL).T


def shells(span=4, rmax=12.0):
    bonds = []
    for n in np.array(np.meshgrid(*[np.arange(-span, span + 1)] * 3, indexing='ij')).reshape(3, -1).T:
        for kind, d in (('same', n.astype(float)), ('inter', n + TAU2 - TAU1)):
            if np.allclose(d, 0):
                continue
            bonds.append((round(float(np.linalg.norm(d @ CELL)), 4), kind, d))
    bonds.sort(key=lambda b: (b[0], b[1]))
    out = []
    for dist, kind, d in bonds:
        if dist > rmax:
            break
        if out and abs(out[-1]['distance'] - dist) < 2e-3 and out[-1]['kind'] == kind:
            out[-1]['vectors'].append(d)
        else:
            out.append({'distance': dist, 'kind': kind, 'vectors': [d]})
    return [dict(s, vectors=np.array(s['vectors'])) for s in out]


SHELLS = shells()


def design(q, n):
    """Uniform two-sublattice spiral: per-cell energy = 2 * sum_s J_s sum_{d in s} (1 - cos 2 pi q.d)."""
    q = np.atleast_2d(np.asarray(q, float))
    return 2 * np.array([(1 - np.cos(2 * np.pi * q @ s['vectors'].T)).sum(1) for s in SHELLS[:n]]).T


def site_split(shell):
    """Does a same-sublattice shell split into Mn1-inequivalent bond classes? Mn1 site symmetry mmm with axes
    [110], [1-10], [001] (rutile P4_2/mnm); the 4_2 operation maps Mn1 bonds d to Mn2 bonds R_90 d."""
    ops = []
    for sx in (1, -1):
        for sz in (1, -1):
            ops += [np.diag([sx, sx, sz]), np.array([[0, sx, 0], [sx, 0, 0], [0, 0, sz]])]
    v = shell['vectors']
    orbit = {tuple(np.round(o @ v[0], 6)) for o in ops}
    return shell['kind'] == 'same' and len(orbit) < len(v)


# ---------------------------------------------------------------- Monte Carlo (classical, unit vectors)
def mc_bonds(J):
    out = []
    for Js, s in zip(J, SHELLS):
        for d in s['vectors']:
            if s['kind'] == 'same':
                out += [(0, 0, *np.round(d).astype(int), Js), (1, 1, *np.round(d).astype(int), Js)]
            else:
                out += [(0, 1, *np.round(d + TAU1 - TAU2).astype(int), Js), (1, 0, *np.round(-d - TAU1 + TAU2).astype(int), Js)]
    return out


def mc_run(args):
    L, T_K, bonds, neq, nmeas, seed = args
    rng = np.random.default_rng(seed)
    xs, ys, zs = [a.ravel() for a in np.meshgrid(np.arange(L), np.arange(L), np.arange(L), indexing='ij')]
    classes = []
    for s in (0, 1):
        bl = [b for b in bonds if b[0] == s]
        Js = np.array([b[5] for b in bl])
        nbr = np.stack([((b[1] * L + (xs + b[2]) % L) * L + (ys + b[3]) % L) * L + (zs + b[4]) % L for b in bl], 1)
        own = ((s * L + xs) * L + ys) * L + zs
        key = (xs % 2) * 8 + (ys % 2) * 4 + zs % 4          # no bond joins two sites of one class
        classes += [(own[key == c], nbr[key == c], np.broadcast_to(Js, (int((key == c).sum()), len(Js)))) for c in range(16)]
    n = L ** 3
    S = np.zeros((2 * n, 3))
    S[:n, 2], S[n:, 2] = 1, -1
    beta = 1 / (T_K * KB)
    m2s = []
    for k in range(neq + nmeas):
        for own, nbr, Jb in classes:                               # heat bath
            h = 2 * np.einsum('nb,nbk->nk', Jb, S[nbr])
            hn = np.linalg.norm(h, axis=1)
            u = h / hn[:, None]
            a = beta * hn
            r = rng.random(len(a))
            c = np.clip(1 + np.log(r + (1 - r) * np.exp(-2 * a)) / a, -1, 1)
            ph = 2 * np.pi * rng.random(len(a))
            v = np.where(np.abs(u[:, 2:3]) < 0.9, np.stack([-u[:, 1], u[:, 0], 0 * a], 1), np.stack([0 * a, -u[:, 2], u[:, 1]], 1))
            v /= np.linalg.norm(v, axis=1)[:, None]
            w = np.cross(u, v)
            sn = np.sqrt(1 - c * c)
            S[own] = c[:, None] * u + sn[:, None] * (np.cos(ph)[:, None] * v + np.sin(ph)[:, None] * w)
        for _ in range(2):                                         # over-relaxation
            for own, nbr, Jb in classes:
                h = np.einsum('nb,nbk->nk', Jb, S[nbr])
                u = h / np.linalg.norm(h, axis=1)[:, None]
                e = S[own]
                S[own] = 2 * (e * u).sum(1)[:, None] * u - e
        if k >= neq:
            m = (S[:n].sum(0) - S[n:].sum(0)) / (2 * n)
            m2s.append(m @ m)
    m2s = np.array(m2s)
    blocks = m2s[: len(m2s) // 10 * 10].reshape(10, -1)
    g = np.mean(m2s ** 2) / np.mean(m2s) ** 2
    gb = [(b ** 2).mean() / b.mean() ** 2 for b in blocks]
    return {'L': L, 'T_K': T_K, 'g': float(g), 'g_err': float(np.std(gb) / 3)}


def neel_temperature(J, t_mf):
    """Binder-ratio crossing of L = 8 and 12; the window starts at 0.70-0.80 T_MF and is shifted if needed."""
    bonds = mc_bonds(J)
    factors = [0.70, 0.725, 0.75, 0.775, 0.80]
    for attempt in range(4):
        temps = [round(t_mf * f, 2) for f in factors]
        jobs = [(L, T, bonds, 300, 1500, 17 * L + i + 100 * attempt) for L in (8, 12) for i, T in enumerate(temps)]
        with Pool(min(2, os.cpu_count() or 1)) as pool:
            res = pool.map(mc_run, jobs)
        g8 = np.array([r['g'] for r in res if r['L'] == 8])
        g12 = np.array([r['g'] for r in res if r['L'] == 12])
        diff = g8 - g12                       # > 0 in the ordered phase, < 0 above T_N
        cross = np.where((diff[:-1] > 0) & (diff[1:] <= 0))[0]
        if len(cross):
            k = int(cross[0])
            T = np.array(temps)
            return float(T[k] + (T[k + 1] - T[k]) * diff[k] / (diff[k] - diff[k + 1])), res
        shift = 0.1 if diff[-1] > 0 else -0.1     # crossing above (below) the window
        factors = [f + shift for f in factors]
    raise RuntimeError('no Binder crossing found')


# ---------------------------------------------------------------- main
def main():
    keys = {(s['direction'], round(s['t'], 2)): s['key'] for s in CATALOG['spirals_available']}
    qof = {s['key']: np.array(s['q_reciprocal_lattice']) for s in CATALOG['spirals_available']}
    bought, evidence = {}, []

    def buy(key):
        r = lab('evaluate', {'kind': 'spiral', 'key': key})
        bought[key] = r['result']
        evidence.append(r['id'])
        return r['result']

    g = buy(keys[('Gamma', 0.0)])
    neel = buy(keys[('[001]', 1.0)])
    # The archive documents energies as meV/atom; compare with the raw total energies of the 6-atom cell.
    per_cell = 1000 * (neel['vasp_total_energy_eV'] - g['vasp_total_energy_eV'])
    assert abs(per_cell - neel['archive_csv_energy']) < 0.01      # archive numbers are per 6-atom cell (2 Mn)
    # Spiral basis: VASP asks for ENMAX > X (X grows with |q|); spirals with X > ENCUT are basis-truncated.
    enini, encut = float(g['incar']['ENINI']), float(g['incar']['ENCUT'])

    def needed(q):
        return (np.sqrt(enini) + np.sqrt(HB2M) * np.linalg.norm(q @ BREC) / 2) ** 2

    for res in (g, neel):
        x = float(re.findall(r'choose ENMAX >\s*([\d.e+]+)\s*eV', ' '.join(res['outcar_warnings']))[0])
        assert abs(needed(np.array(res['q_reciprocal_lattice'])) - x) < 0.01 * x + 2
    plan = {'[001]': [0.15, 0.35, 0.55, 0.75, 0.9], '[121]': [0.1, 0.2, 0.3, 0.4, 0.5],
            '[214]': [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85]}
    excluded = [k for (d, t), k in keys.items() if needed(qof[k]) > encut]
    for d, ts in plan.items():
        for t in ts:
            assert keys[(d, t)] not in excluded
            buy(keys[(d, t)])
    q = np.array([bought[k]['q_reciprocal_lattice'] for k in bought])
    y = np.array([1000 * (bought[k]['vasp_total_energy_eV'] - g['vasp_total_energy_eV']) for k in bought])
    # shell count by leave-one-out cross-validation (in-sample error always falls with more shells)
    cv = {}
    for n in range(3, 13):
        A = np.c_[np.ones(len(y)), design(q, n)]
        H = A @ np.linalg.pinv(A)
        cv[n] = float(np.sqrt(np.mean(((y - H @ y) / np.clip(1 - np.diag(H), 1e-9, None)) ** 2)))
    n = min(cv, key=cv.get)
    A = np.c_[np.ones(len(y)), design(q, n)]
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    c0, J = coef[0], coef[1:]
    preds = {p['id']: float(c0 + design(p['q_reciprocal_lattice'], n)[0] @ J) for p in CATALOG['spirals_to_predict']}
    # Neel temperature of the fitted model (classical MC, Binder crossing L = 8 / 12)
    t_mf = 2 / 3 * sum(Js * len(s['vectors']) * (1 if s['kind'] == 'same' else -1) for Js, s in zip(J, SHELLS)) / KB
    tn, mc = neel_temperature(J, t_mf)
    split = next(i for i, s in enumerate(SHELLS) if site_split(s))
    submission = {'spiral_predictions_meV_per_cell': preds, 'J1_meV': float(J[0]), 'J2_meV': float(J[1]), 'J3_meV': float(J[2]),
                  'neel_temperature_K': tn,
                  'altermagnetic_shell': {'distance_A': SHELLS[split]['distance'], 'sublattice_mean_J_meV': float(J[split]),
                                          'sublattice_difference_J_meV': None},
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    shell_txt = ', '.join(f"{s['distance']:.2f} A {s['kind']} J={j:+.3f}" for s, j in zip(SHELLS, J))
    (APP/'report.md').write_text(f"""# MnF2 exchange model from released spin spirals (reference solution)

Conventions: H = -sum_(i!=j) J_ij e_i.e_j over Mn unit vectors, ordered pairs. The archive documents its energies as
meV/atom, but 1000 x (difference of the VASP total energies of the 6-atom cell) reproduces them exactly, so they are
per primitive cell (2 Mn): E_cell(q) - E_cell(0) = 2 sum_s J_s sum_(d in s) (1 - cos 2 pi q.d).

Spiral geometry: the INCARs set the in-cell angle of Mn2 to -2 pi q.tau(Mn1), so with VASP's +2 pi q.R rotation
between cells every Mn carries the azimuth 2 pi q.r: a uniform spiral on the body-centred Mn lattice. q = 0 is the
ferromagnet; q = (0,0,1) (any integer q with h+k+l odd) is the Neel state, the lowest energy found
({neel['archive_csv_energy']:.2f} meV per cell). E(q) has the periodicity of the bct reciprocal lattice, so [121]
spirals are symmetric about t = 0.5.

Data screening: every OUTCAR prints VASP's 'choose ENMAX > X' for the spiral basis; X follows
(sqrt(ENINI) + sqrt(hbar^2/2m)|q|/2)^2 and exceeds ENCUT = {encut:.0f} eV for {', '.join(sorted(excluded))}. Those
spirals are basis-truncated (the [214] t = 1 spiral is the Neel state yet differs from the [001] Neel energy by
~1.5 meV) and were not bought or used; the deviation is numerical, not physics.

Purchases: Gamma, [001], [121] (t <= 0.5 by symmetry) and [214] up to t = 0.85 ({len(evidence)} spirals). [001] alone
aliases the 3.89 A and 8.01 A inter-sublattice shells and cannot see in-plane bonds; [001] + [121] without [214] is
rank deficient beyond 8 shells; commensurate q grids alias long-range shells. Shell count by leave-one-out CV
(errors {json.dumps({k: round(v, 3) for k, v in cv.items()})}): {n} shells. In-sample error barely improves from 3 to
6 shells; the 7.0 and 7.8 A same-sublattice shells are needed. Model: {shell_txt}.
J1 = {J[0]:.3f}, J2 = {J[1]:.3f}, J3 = {J[2]:.3f} meV.

Neel temperature: classical heat-bath + over-relaxation Monte Carlo of the fitted model (unit vectors), Binder ratio
of the staggered magnetisation for L = 8 and 12 (2 L^3 spins): crossing at T_N = {tn:.0f} K (statistical and
finite-size uncertainty ~3%). Mean field gives {t_mf:.0f} K and overestimates. T_N is insensitive (second order) to
the undetermined altermagnetic difference below.

Altermagnetism: the first Mn-Mn shell whose bonds are sublattice-inequivalent is the same-sublattice a<110> shell at
{SHELLS[split]['distance']:.3f} A (Mn1 site symmetry mmm with [110]/[1-10] axes; the 4_2 operation swaps [110] and [1-10]
between sublattices). In any uniform spiral both sublattices see the same q, so the energy contains only
J(Mn1,d) + J(Mn2,d): the sublattice mean ({J[split]:+.3f} meV) is determined, the difference - and hence the magnon
chirality splitting - is not (null), not zero. Resolving it needs configurations in which the two sublattices see
different spin patterns along [110], e.g. collinear supercells with one sublattice ferromagnetic and the other
modulated along [110], or two-q spirals with sublattice-dependent q. These data do not support any statement about
the size or sign of the altermagnetic splitting, nor a validated comparison with the experimental T_N (PBE, classical
spins).
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
