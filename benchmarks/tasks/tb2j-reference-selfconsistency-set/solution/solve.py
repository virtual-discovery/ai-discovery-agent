"""Reference solution for tb2j-reference-selfconsistency-set; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
WINDOWS = list(CATALOG['requests']['exchange']['windows_A'])
evidence = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    while result.get('status') == 'running':
        time.sleep(1)
        result = lab('job/%d' % result['id'])
    return result


def buy(calc, windows):
    meta, pairs = None, []
    for w in windows:
        r = lab('evaluate', {'kind': 'exchange', 'calculation': calc, 'window': w})
        evidence.append(r['id'])
        meta = r['result']['metadata']
        pairs += r['result']['pairs']
    return {'meta': meta, 'pairs': pairs}


def upto(windows, hi):
    return [w for w in windows if float(w.split('-')[1]) <= hi]


class Model:
    """Isotropic TB2J model (E = -sum_{i!=j} J e_i.e_j, ordered pairs) in a calculation's own cell."""

    def __init__(self, d, keep=lambda p: True):
        mag = [a for a in d['meta']['atoms'] if a['name'].rstrip('0123456789') in ('Ni', 'Fe')]
        self.names = [a['name'] for a in mag]
        self.n = len(mag)
        self.cell = np.array(d['meta']['cell_A'])
        self.pos = np.array([a['position_A'] for a in mag])
        m = np.array([[a['M(x)'], a['M(y)'], a['M(z)']] if 'M(z)' in a else [0, 0, a['w_magmom']] for a in mag])
        self.e = m/np.linalg.norm(m, axis=1, keepdims=True)
        idx = {n: k for k, n in enumerate(self.names)}
        P = [p for p in d['pairs'] if keep(p)]
        self.setup(np.array([idx[p[0]] for p in P]), np.array([idx[p[1]] for p in P]),
                   np.array([p[2:5] for p in P], float), np.array([p[5] for p in P]))

    def setup(self, I, J, R, V):
        self.I, self.J, self.R, self.V = I, J, R.reshape(-1, 3), V
        self.A = sp.csr_matrix((np.ones(len(V)), (np.arange(len(V)), I*self.n + J)), shape=(len(V), self.n**2))

    def K(self, k):
        k = np.atleast_2d(np.asarray(k, float))
        ph = np.exp(2j*np.pi*(k @ self.R.T))*self.V
        return np.asarray((self.A.T @ ph.T).T).reshape(-1, self.n, self.n)

    def K0(self):
        return self.K([0, 0, 0])[0].real

    def h(self, k, e=None):
        """Bosonic LSWT matrix (unit spin length) around directions e (Toth-Lake local frames)."""
        e = self.e if e is None else e
        u = []
        for z in e:
            t = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
            x = t - (t @ z)*z
            x /= np.linalg.norm(x)
            u.append(x + 1j*np.cross(z, x))
        u = np.array(u)
        K0 = self.K0()
        C = np.diag(-(K0*(e @ e.T)).sum(1))
        out = []
        for Kq, Kmq in zip(self.K(k), self.K(-np.atleast_2d(k))):
            A, B, Am = -0.5*Kq*(u @ u.conj().T), -0.5*Kq*(u @ u.T), -0.5*Kmq*(u @ u.conj().T)
            out.append(2*np.block([[A - C, B], [B.conj().T, Am.conj() - C]]))
        return np.array(out)


def grid(n):
    g = np.arange(n)/n
    return np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T


def kmin(model):
    """Minimum over k of the lowest eigenvalue of h(k) around the ferromagnetic configuration."""
    fm = np.tile([0, 0, 1.0], (model.n, 1))
    low = lambda k: float(np.linalg.eigvalsh(model.h(np.atleast_2d(k), fm))[0, 0])
    ks = grid(10)
    ev = np.linalg.eigvalsh(model.h(ks, fm))[:, 0]
    starts = list(ks[np.argsort(ev)[:6]]) + [np.array([0, 0.03, 0.2]), np.array([0, 0.5, 0])]
    best = min((minimize(low, s, method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 4000})
                for s in starts), key=lambda r: r.fun)
    x = best.x - np.round(best.x)
    return np.abs(x), best.fun


def site_map(fm, sc):
    """Primitive site and FM-lattice translation of each supercell site (supercell matrix by rounding)."""
    Msc = np.round(sc.cell @ np.linalg.inv(fm.cell)).astype(int)
    inv_eff = np.linalg.inv(np.linalg.inv(Msc) @ sc.cell)
    pf = fm.pos @ np.linalg.inv(fm.cell)
    pmap, tvec = [], []
    for f in sc.pos @ inv_eff:
        d = f - pf
        a = int(np.argmin(np.linalg.norm(d - np.round(d), axis=1)))
        pmap.append(a)
        tvec.append(np.round(f - pf[a]).astype(int))
    return np.array(pmap), np.array(tvec), Msc, inv_eff, pf


def fm_couplings_on_supercell(fm, sc):
    pmap, tvec, Msc, _, _ = site_map(fm, sc)
    Minv = np.linalg.inv(Msc)
    key = lambda a, t: (int(a), tuple(np.round(((t @ Minv) % 1.0)*1000).astype(int) % 1000))
    lookup = {key(a, t): i for i, (a, t) in enumerate(zip(pmap, tvec))}
    K0 = np.zeros((sc.n, sc.n))
    for i in range(sc.n):
        sel = fm.I == pmap[i]
        for b, R, V in zip(fm.J[sel], fm.R[sel].astype(int), fm.V[sel]):
            K0[i, lookup[key(b, tvec[i] + R)]] += V
    return K0


def energy(K0, e):
    return float(-np.einsum('ij,ij->', K0, e @ e.T)/len(e))


def stabilization(fm_d, sc_d):
    fm, sc = Model(fm_d), Model(sc_d)
    ones = np.tile([0, 0, 1.0], (sc.n, 1))
    Kf, Ka = fm_couplings_on_supercell(fm, sc), sc.K0()
    return energy(Kf, ones) - energy(Kf, sc.e), energy(Ka, ones) - energy(Ka, sc.e)


def long_axis_aliased(d):
    """True when every pair of a supercell list has R = 0 along the longest cell vector."""
    ax = int(np.argmax(np.linalg.norm(d['meta']['cell_A'], axis=1)))
    return all(p[2 + ax] == 0 for p in d['pairs']), ax


def minimum_image(d, ax):
    cell = np.array(d['meta']['cell_A'])
    inv = np.linalg.inv(cell)
    out = []
    for p in d['pairs']:
        v = np.array(p[6:9])
        n = -int(np.round((v @ inv)[ax]))
        v = v + n*cell[ax]
        R = list(p[2:5])
        R[ax] += n
        out.append([p[0], p[1], *R, p[5], *v.tolist(), float(np.linalg.norm(v))])
    return {'meta': d['meta'], 'pairs': out}


def folded(fm, sc_d):
    """Supercell couplings mapped onto the primitive cell; bonds equivalent under the supercell translations averaged."""
    sc = Model(sc_d)
    pmap, tvec, Msc, inv_eff, pf = site_map(fm, sc)
    idx = {n: k for k, n in enumerate(sc.names)}
    acc = defaultdict(list)
    for p in sc_d['pairs']:
        a, b = pmap[idx[p[0]]], pmap[idx[p[1]]]
        T = np.round(np.array(p[6:9]) @ inv_eff - (pf[b] - pf[a])).astype(int)
        acc[(int(a), int(b), tuple(T))].append(p[5])
    F = Model.__new__(Model)
    F.n = fm.n
    keys = list(acc)
    F.setup(np.array([k[0] for k in keys]), np.array([k[1] for k in keys]), np.array([k[2] for k in keys], float),
            np.array([np.mean(acc[k]) for k in keys]))
    return F, np.bincount([len(v) for v in acc.values()]).tolist()


def relax(K0, e):
    e = e.copy()
    Kx = K0 - np.diag(np.diag(K0))
    for _ in range(5000):
        for i in range(len(e)):
            f = Kx[i] @ e
            e[i] = f/np.linalg.norm(f)
    return e


def main():
    W = WINDOWS
    d = {}
    # insulators: FM and updated references with the short/intermediate windows
    d['NiO/FM'] = buy('NiO/FM', upto(W, 8))
    d['NiO/AFM'] = buy('NiO/AFM', upto(W, 8))
    d['FeCl2/FM'] = buy('FeCl2/FM', upto(W, 14))
    d['FePS3/DFT/FM'] = buy('FePS3/DFT/FM', upto(W, 18))
    d['FePS3/DFT/AFM'] = buy('FePS3/DFT/AFM', upto(W, 18))
    d['FePS3/DFT+U/AFM'] = buy('FePS3/DFT+U/AFM', upto(W, 11))
    # provenance check of the DFT+U FM reference against the U = 0 one (first window)
    d['FePS3/DFT+U/FM'] = buy('FePS3/DFT+U/FM', upto(W, 5))
    a, b = d['FePS3/DFT+U/FM'], d['FePS3/DFT/FM']
    sig = lambda pairs: sorted(tuple(p[:6]) for p in pairs if p[9] < 5)
    duplicate = (sig(a['pairs']) == sig(b['pairs']) and a['meta']['atoms'] == b['meta']['atoms']
                 and a['meta']['cell_A'] == b['meta']['cell_A'])
    moments = {k: max(abs(x.get('M(z)', x.get('w_magmom', 0))) for x in d[k]['meta']['atoms']) for k in ('FePS3/DFT/FM', 'FePS3/DFT+U/FM', 'FePS3/DFT/AFM', 'FePS3/DFT+U/AFM')}
    # FeCl2 AFM: the supercell list is computed with R = 0 along c only; bonds across the boundary appear at the far
    # image and long archived distances, so every window up to the supercell length is needed
    d['FeCl2/AFM'] = buy('FeCl2/AFM', upto(W, 14))
    alias_fecl2, ax_fecl2 = long_axis_aliased(d['FeCl2/AFM'])
    if alias_fecl2:
        L = np.linalg.norm(d['FeCl2/AFM']['meta']['cell_A'][ax_fecl2])
        more = buy('FeCl2/AFM', [w for w in W if 14 <= float(w.split('-')[0]) < L])
        d['FeCl2/AFM']['pairs'] += more['pairs']
    # FeP: FM reference with an exchange-range convergence test, then the spiral (AFM) reference in full
    d['FeP/FM'] = buy('FeP/FM', upto(W, 30))
    conv = {}
    for hi in (14, 18, 23, 30):
        conv[hi] = kmin(Model(d['FeP/FM'], keep=lambda p, hi=hi: p[9] < hi))
    d['FeP/AFM'] = buy('FeP/AFM', upto(W, 40))
    alias_fep, ax_fep = long_axis_aliased(d['FeP/AFM'])

    k_fm, low_fm = conv[30]
    fm = Model(d['FeP/FM'])
    spiral = minimum_image(d['FeP/AFM'], ax_fep) if alias_fep else d['FeP/AFM']
    F, copies = folded(fm, spiral)
    k_conv, low_conv = kmin(F)
    F_listed, copies_listed = folded(fm, d['FeP/AFM'])
    k_listed, _ = kmin(F_listed)

    # stabilisation energies, E(FM) - E(converged configuration) per magnetic ion, each exchange set
    pairs_of = {'NiO': ('NiO/FM', 'NiO/AFM'), 'FeCl2': ('FeCl2/FM', 'FeCl2/AFM'), 'FePS3_DFT': ('FePS3/DFT/FM', 'FePS3/DFT/AFM'),
                'FePS3_DFT+U': ('FePS3/DFT+U/FM', 'FePS3/DFT+U/AFM'), 'FeP': ('FeP/FM', 'FeP/AFM')}
    E = {}
    for mat, (f, s) in pairs_of.items():
        if mat == 'FePS3_DFT+U':
            E[mat] = {'fm_reference': None, 'converged_reference': stabilization(d['FePS3/DFT/FM'], d[s])[1]}
            continue
        ef, ea = stabilization(d[f], d[s])
        E[mat] = {'fm_reference': ef, 'converged_reference': ea}
    if not duplicate:
        E['FePS3_DFT+U']['fm_reference'] = stabilization(d['FePS3/DFT+U/FM'], d['FePS3/DFT+U/AFM'])[0]
    ratio = {m: (None if e['fm_reference'] is None else e['converged_reference']/e['fm_reference']) for m, e in E.items()}
    short = {'meta': d['FeCl2/AFM']['meta'], 'pairs': [p for p in d['FeCl2/AFM']['pairs'] if p[9] < 30]}
    fecl2_loss = 1 - stabilization(d['FeCl2/FM'], short)[1]/E['FeCl2']['converged_reference']
    nn = lambda k: np.mean([p[5] for p in d[k]['pairs'] if 2.6 < p[9] < 2.7])
    nn_ratio = nn('FeP/AFM')/nn('FeP/FM')
    netsum = {k: Model(d[k]).K0().sum()/Model(d[k]).n for k in ('FeP/FM', 'FeP/AFM')}

    # stability of each calculation at its own reference; FeP spiral: torque and relaxation
    stab = {}
    for k in d:
        if k == 'FePS3/DFT+U/FM':
            continue
        m = Model(spiral if k == 'FeP/AFM' else d[k])
        stab[k] = float(np.linalg.eigvalsh(m.h(grid(5))).min())
    sp_model = Model(spiral)
    K0 = sp_model.K0()
    torque = np.linalg.norm(np.cross(sp_model.e, K0 @ sp_model.e), axis=1)
    e_rel = relax(K0, sp_model.e)
    stab_relaxed = float(np.linalg.eigvalsh(sp_model.h(grid(6), e_rel)).min())
    stab_relaxed_listed = float(np.linalg.eigvalsh(Model(d['FeP/AFM']).h(grid(6), e_rel)).min())

    submission = {'fm_reference_kmin_FeP_rlu': [float(x) for x in k_fm],
                  'converged_exchange_kmin_FeP_rlu': [float(x) for x in k_conv],
                  'stabilization_energy_meV_per_ion': E, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    fmt = lambda v: 'null' if v is None else '%.4g' % v
    table = '\n'.join('| %s | %s | %s | %s |' % (m, fmt(e['fm_reference']), fmt(e['converged_reference']), fmt(ratio[m]))
                      for m, e in E.items())
    (APP/'report.md').write_text(f"""# Reference-state dependence of TB2J exchange (reference solution)

## Conventions and handling of the calculations
Archived TB2J convention E = -sum_(i!=j) J_ij e_i.e_j (ordered pairs, unit vectors, J > 0 FM). Each calculation is
used in its own cell with its own moments. Stabilisation energies dE = E(FM) - E(converged configuration) per
magnetic ion, classical, one exchange set at a time; the FM-reference couplings are carried onto the converged
supercell by lattice translation (the small structural differences between the calculations are ignored). LSWT uses
the bosonic matrix h(k) (Toth-Lake local frames, unit spin length); k_min minimises its lowest eigenvalue around the
FM state.

## Self-consistency tests
Lowest h eigenvalue at each calculation's own moments (meV): {json.dumps({k: round(v, 3) for k, v in stab.items()})}.
All FM references are unstable; the AFM references of NiO, FeCl2 and FePS3 (both U) are positive semi-definite
(Goldstone zero only), so they close the loop: AF-II (L) in NiO, FM layers stacked antiferromagnetically in FeCl2,
zigzag chains with antiparallel stacking in FePS3. The archived FeP spiral moments are not a stationary point of their
own exchange (torques {torque.min():.1f}-{torque.max():.1f} meV), which produces a spurious negative eigenvalue at Gamma
(also present in the authors' archived h(k)); after relaxing the directions in their exchange field (same 1x1x5 cell)
the lowest eigenvalue is {stab_relaxed:.3f} meV with minimum-image couplings ({stab_relaxed_listed:.1e} meV as listed).

## Exchange range
FeP FM-reference k_min (|components|, r.l.u. of the four-Fe cell) with pairs below 14/18/23/30 A:
{[np.round(conv[h][0], 4).tolist() for h in (14, 18, 23, 30)]}. Couplings beyond ~15 A set the helix (a 10-11 A list gives
a b-axis modulation instead); the energies converge only with the 18-30 A windows. Final k_min = {np.round(k_fm, 4).tolist()}
(lowest eigenvalue {low_fm:.2f} meV), close to but not equal to the commensurate 1x1x5 approximant (0, 0, 0.2).

## Reference dependence
| material | dE FM-reference J | dE converged J | ratio |
|---|---|---|---|
{table}

The insulators NiO and FeCl2 are reference independent to 1-3 %; FePS3 (U = 0) changes by about 11 %, but its AFM
calculation also used a different (zigzag-distorted) structure, so reference and lattice effects are confounded.
Metallic FeP is strongly reference dependent: with the spiral-reference couplings the helix is stabilised about
{ratio['FeP']:.1f} times more strongly; the net coupling sum per Fe changes from {netsum['FeP/FM']:+.1f} to {netsum['FeP/AFM']:+.1f} meV and
the 2.66 A bond becomes {nn_ratio:.1f} times more antiferromagnetic.

## FeP self-consistency
The FeP and FeCl2 AFM supercell lists were computed with R = 0 along the long axis only (aliased: {alias_fep}, {alias_fecl2}):
bonds across the cell boundary are listed at the far image (e.g. the -4.2 meV 3.81 A bond appears at 25.7 A), and the
authors' archived FeP h(k) is exactly flat along Gamma-Z. For energies every listed pair counts (FeCl2 loses {100*fecl2_loss:.0f} % of\nits stacking energy without the 30-40 A window). For the folded k_min each pair was moved to its minimum image and
the five copies of each bond averaged (copies per bond: {copies}); the converged exchange then gives
k_min = {np.round(k_conv, 4).tolist()} (as listed it would give {np.round(k_listed, 4).tolist()}, an artefact). The
spiral-reference exchange prefers a slightly longer period than its own 1x1x5 reference (q_c ~ 0.19 versus 0.2 and the
FM-reference 0.218); the residual instability of the relaxed 1x1x5 spiral ({stab_relaxed:.3f} meV) is below the paper's
0.5 meV window, so the commensurate loop is closed.

## Data problems found
The FePS3 DFT+U FM file duplicates the U = 0 FM calculation (same cell, moments {moments['FePS3/DFT+U/FM']:.3f} muB, identical
couplings; the DFT+U AFM has {moments['FePS3/DFT+U/AFM']:.3f} muB and a different cell): duplicate = {duplicate}. The FM-reference
energy at U = 3 eV is therefore not determined (null), and the DFT+U AFM versus 'DFT+U FM' difference is a U effect, not a
reference effect.

## Not supported by these data
Isotropic couplings only (FeP canting needs DMI/anisotropy); no DFT total energies to test the mapping directly;
classical energies, no transition temperatures; reference effects are confounded with the structural relaxation of each
reference; the incommensurate FeP pitch beyond the 1x1x5 approximant is not resolved by the supercell data.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
