"""Reference solution for amaranta-2d-triangular-helimagnets; runs as the Harbor oracle and in check_task.

Per compound: enumerate bond classes by length up to ~2.05 x the nearest-neighbour distance
(this separates the two nearest-neighbour lengths of the distorted cells), pick the cheapest
diagonal supercell in which no representative pair has a periodic image within that range,
map every class with collinear four-state energies sharing one reference spin, find the
Luttinger-Tisza minimum, classify it, and predict the hidden configurations. The VI2 easy axis
is fitted to stripe energies along nine directions (off-axis probes included).
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
MAXC = 6
evidence, spent = [], [0]


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    for _ in range(3):
        out = subprocess.run(cmd, capture_output=True, text=True).stdout
        try:
            result = json.loads(out)
            break
        except json.JSONDecodeError:
            time.sleep(2)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def energies(compound, n1, n2, configs):
    out = []
    for s in range(0, len(configs), MAXC):
        r = lab('evaluate', {'kind': 'energies', 'compound': compound, 'supercell': [n1, n2],
                             'spins': [np.round(c, 6).tolist() for c in configs[s:s + MAXC]]})
        while r['status'] == 'running':
            time.sleep(1)
            r = lab('job/%d' % r['id'])
        if r['status'] != 'complete':
            raise RuntimeError(r)
        evidence.append(r['id'])
        out += [e * 1000.0 for e in r['result']['total_energies_eV']]   # meV
    return out


class Lattice:
    def __init__(self, comp):
        self.cell = np.array(comp['cell_A'], float)[:2, :2]
        self.basis = np.array(comp['magnetic_sites'], float)[:, :2]
        self.nb = len(self.basis)
        self.B = 2 * np.pi * np.linalg.inv(self.cell).T
        self.frac = [np.linalg.solve(self.cell.T, b - self.basis[0]) for b in self.basis]
        d = [np.linalg.norm(b + i * self.cell[0] + j * self.cell[1] - self.basis[0]) for b in self.basis
             for i in range(-2, 3) for j in range(-2, 3)]
        self.nn = min(x for x in d if x > 0.1)
        self.rcut = 2.05 * self.nn

    def vectors(self, k0=0):
        out = []
        for k, b in enumerate(self.basis):
            for m1 in range(-5, 6):
                for m2 in range(-5, 6):
                    d = b + m1 * self.cell[0] + m2 * self.cell[1] - self.basis[k0]
                    if 0.1 < np.linalg.norm(d) <= self.rcut:
                        out.append((k, m1, m2, d))
        return out

    def classes(self):
        lengths = []
        for k, m1, m2, d in self.vectors():
            L = np.linalg.norm(d)
            if not any(abs(L - x) < 0.005 for x in lengths):
                lengths.append(L)
        return sorted(lengths)

    def cls(self, L, lengths):
        for i, x in enumerate(lengths):
            if abs(L - x) < 0.005:
                return i
        return None

    def primitive_G(self):
        return [np.array([h, k], float) for h in range(-6, 7) for k in range(-6, 7)
                if all(abs((h * f[0] + k * f[1] + 0.5) % 1.0 - 0.5) < 1e-3 for f in self.frac)]


def choose_supercell(lat, lengths):
    """Cheapest diagonal supercell in which every class has an alias-free representative."""
    vecs = lat.vectors()
    best = None
    for n1 in range(1, 37):
        for n2 in range(1, 37):
            N = lat.nb * n1 * n2
            if N > 36:
                continue
            L1, L2 = n1 * lat.cell[0], n2 * lat.cell[1]
            reps = []
            for ci, L in enumerate(lengths):
                ok = None
                for k, m1, m2, d in vecs:
                    if abs(np.linalg.norm(d) - L) > 0.005:
                        continue
                    if all(np.linalg.norm(d + u * L1 + v * L2) > lat.rcut + 0.05
                           for u in range(-3, 4) for v in range(-3, 4) if (u, v) != (0, 0)):
                        ok = (k, m1, m2)
                        break
                if ok is None:
                    break
                reps.append(ok)
            if len(reps) < len(lengths):
                continue
            cost = math.ceil(N / 12) * (2 + 2 * len(lengths))
            key = (cost, -N)
            if best is None or key < best[0]:
                best = (key, n1, n2, reps)
    return best[1], best[2], best[3], best[0][0]


def map_exchange(name, lat):
    lengths = lat.classes()
    n1, n2, reps, cost = choose_supercell(lat, lengths)
    N = lat.nb * n1 * n2
    idx = lambda k, m1, m2: (k * n1 + m1 % n1) * n2 + m2 % n2
    up = np.tile([0.0, 0.0, 1.0], (N, 1))
    def flipped(sites):
        c = up.copy()
        c[list(sites), 2] = -1.0
        return c
    js = [idx(*r) for r in reps]
    configs = [up, flipped([0])] + [flipped([j]) for j in js] + [flipped([0, j]) for j in js]
    E = energies(name, n1, n2, configs)
    nc = len(js)
    J = [(E[0] - E[1] - E[2 + c] + E[2 + nc + c]) / 4 for c in range(nc)]
    return lengths, J, (n1, n2)


def Jq(lat, lengths, J, q):
    e = 0.0
    for k, m1, m2, d in lat.vectors():
        c = lat.cls(np.linalg.norm(d), lengths)
        e += 0.5 * J[c] * math.cos(q @ d)
    return e


def ground_state(lat, lengths, J):
    best = None
    for h in np.linspace(0, 1, 61):
        for k in np.linspace(0, 2, 121):
            e = Jq(lat, lengths, J, np.array([h, k]) @ lat.B)
            if best is None or e < best[0]:
                best = (e, h, k)
    r = minimize(lambda x: Jq(lat, lengths, J, np.array(x) @ lat.B), best[1:], method='Nelder-Mead',
                 options={'xatol': 1e-9, 'fatol': 1e-12})
    G = lat.primitive_G()
    q = min((np.array(r.x) - g for g in G), key=lambda v: np.linalg.norm(v @ lat.B))   # shortest equivalent
    def near(mult):
        return min(np.linalg.norm((g / mult - q) @ lat.B) for g in G)
    def snap(mult):
        g = min(G, key=lambda g: np.linalg.norm((g / mult - q) @ lat.B))
        return g / mult
    tol = 0.01
    if near(1) < tol:
        return 'ferromagnet', [0.0, 0.0], r.fun
    if near(2) < tol and Jq(lat, lengths, J, snap(2) @ lat.B) <= r.fun + 1e-6:
        return 'collinear stripe', snap(2).tolist(), r.fun
    if near(3) < tol and Jq(lat, lengths, J, snap(3) @ lat.B) <= r.fun + 1e-6:
        return '120-degree', snap(3).tolist(), r.fun
    return 'incommensurate spiral', q.tolist(), r.fun


def predict(lat, lengths, J, p):
    n1, n2 = p['supercell']
    S = np.array([[math.sin(math.radians(t)) * math.cos(math.radians(f)), math.sin(math.radians(t)) * math.sin(math.radians(f)),
                   math.cos(math.radians(t))] for t, f in p['spin_angles_deg']])
    N = len(S)
    def energy(spins):
        X = spins.reshape(lat.nb, n1, n2, 3)
        e = 0.0
        for k0 in range(lat.nb):
            for k, m1, m2, d in lat.vectors(k0):
                c = lat.cls(np.linalg.norm(d), lengths)
                Y = np.roll(np.roll(X[k], -m1, axis=0), -m2, axis=1)
                e += 0.5 * J[c] * np.sum(X[k0] * Y)
        return e
    fm = np.tile([0.0, 0.0, 1.0], (N, 1))
    return (energy(S) - energy(fm)) / N


def vi2_axis():
    lat = Lattice(CAT['compounds']['VI2'])
    n1, n2 = 2, 3
    N = lat.nb * n1 * n2
    pattern = np.array([1.0] * (N // 2) + [-1.0] * (N // 2))   # corner sublattice up, centred sublattice down
    dirs = [np.array(v, float) / np.linalg.norm(v) for v in
            ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1], [1, 1, 0], [1, -1, 0])]
    E = energies('VI2', n1, n2, [np.outer(pattern, v) for v in dirs])
    X = np.array([[1, v[0] ** 2, v[1] ** 2, 2 * v[0] * v[1], 2 * v[0] * v[2], 2 * v[1] * v[2]] for v in dirs])
    c = np.linalg.lstsq(X, np.array(E), rcond=None)[0]
    Mm = np.array([[c[1], c[3], c[4]], [c[3], c[2], c[5]], [c[4], c[5], 0.0]]) / N
    w, V = np.linalg.eigh(Mm)
    n = V[:, 0] * np.sign(V[2, 0])
    return n, dict(zip('x y z x+z x-z y+z y-z x+y x-y'.split(), np.round(np.array(E) / N - min(E) / N, 4).tolist())), w


def main():
    results, preds, notes = {}, {}, []
    for name, comp in CAT['compounds'].items():
        lat = Lattice(comp)
        lengths, J, sc = map_exchange(name, lat)
        typ, q, e = ground_state(lat, lengths, J)
        results[name] = {'type': typ, 'q_rlu': [round(x, 6) for x in q]}
        for pid, p in CAT['to_predict'].items():
            if p['compound'] == name:
                preds[pid] = float(predict(lat, lengths, J, p))
        notes.append(f"- **{name}** (supercell {sc[0]}x{sc[1]}): " + ', '.join(f"J({L:.3f} A) = {j:+.3f}" for L, j in zip(lengths, J))
                     + f" meV -> {typ}, q = ({q[0]:.4f}, {q[1]:.4f}) r.l.u., E = {e:.3f} meV/site")
    axis, probe, w = vi2_axis()
    tilt = math.degrees(math.acos(abs(axis[2])))
    submission = {'ground_states': results, 'energy_predictions_meV_per_site': preds,
                  'VI2_moment_direction': [round(float(x), 6) for x in axis], 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Triangular dihalide monolayers (reference solution)

Convention: E = sum over pairs (each once) of e_i.J.e_j for unit spins, plus single-ion terms.
Exchange classes were defined by bond length up to 2.05 x the nearest-neighbour distance: the
rectangular cells of VBr2, VI2 and MnCl2 are distorted and their six nearest neighbours split
into 4 + 2 bonds of different length and very different exchange. Each class was mapped with
collinear four-state energies (shared reference spin) in the cheapest diagonal supercell in
which no pair has a periodic image within the coupling range (a 2a pair in a 4-periodic cell
would be counted twice).

{chr(10).join(notes)}

Ground states follow from the Luttinger-Tisza minimum of the isotropic exchange. Ni dihalides:
ferromagnetic J1 with antiferromagnetic J3 > |J1|/4 gives an incommensurate spiral along the
bond direction. VCl2, MnBr2, MnI2: 120-degree order; in MnBr2/MnI2 J2 exceeds J1/8 (stripe in
a J1-J2 model) but the antiferromagnetic J3 stabilizes the K point. VBr2 and MnCl2: anisotropic
triangular lattice (J' < 2J along b) gives an incommensurate spiral; VI2 (J' > 2J) a collinear
stripe.

VI2 anisotropy: stripe energies along nine directions (meV/site, relative): {probe}. The fitted
quadratic form has its minimum along {np.round(axis, 3).tolist()}, tilted {tilt:.1f} deg from c
in the ac plane: the symmetric off-diagonal J_xz of the nearest-neighbour tensor survives the
mirror average and tilts the axis away from the single-ion easy axis (c). Probing only x, y, z
would miss it.

Finite temperature: the exchange is isotropic to within ~0.1 meV, so by Mermin-Wagner there is
no finite-temperature order without the anisotropy, and any ordering (Ising-like or BKT) is set
by these tiny anisotropies, logarithmically; ordering temperatures cannot be determined
reliably from these data (mean-field values overestimate them). Anisotropies of 0.01-0.1 meV are
at or below the source's stated precision (0.1 meV).
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
