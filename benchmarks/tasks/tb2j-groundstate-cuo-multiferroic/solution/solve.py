"""Reference solution for tb2j-groundstate-cuo-multiferroic; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
CONV = np.array(CATALOG['conventional_cell']['cell_A'])
S = 0.5                                   # Cu2+ spin length for the spin waves
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
    """Purchase distance windows of one calculation; return (metadata, pairs)."""
    meta, pairs = None, []
    for w in windows:
        r = lab('evaluate', {'kind': 'exchange', 'calculation': calc, 'window': w})
        evidence.append(r['id'])
        meta = r['result']['metadata']
        pairs += r['result']['pairs']
    return meta, pairs


class Model:
    """Isotropic model of one calculation in its own cell: K_ij(k) = sum_R J_ij(R) e^{2 pi i k.R}."""

    def __init__(self, meta, pairs, rmax=None):
        cu = [a for a in meta['atoms'] if a['name'].startswith('Cu')]
        self.names = [a['name'] for a in cu]
        self.n = len(cu)
        self.cell = np.array(meta['cell_A'])
        self.M = np.round(self.cell @ np.linalg.inv(CONV)).astype(int)        # supercell of the conventional cell
        self.spins = np.sign([a.get('M(z)', a.get('w_magmom')) for a in cu])
        self.fconv = np.array([a['position_A'] for a in cu]) @ np.linalg.inv(self.cell) @ self.M
        idx = {n: k for k, n in enumerate(self.names)}
        sel = [p for p in pairs if rmax is None or p[9] < rmax]
        self.I = np.array([idx[p[0]] for p in sel])
        self.J = np.array([idx[p[1]] for p in sel])
        self.R = np.array([p[2:5] for p in sel], float)
        self.V = np.array([p[5] for p in sel])

    def K(self, k):
        k = np.atleast_2d(np.asarray(k, float))
        ph = np.exp(2j*np.pi*(k @ self.R.T))*self.V
        out = np.zeros((len(k), self.n, self.n), complex)
        for a in range(self.n):
            for b in range(self.n):
                m = (self.I == a) & (self.J == b)
                if m.any():
                    out[:, a, b] = ph[:, m].sum(1)
        return out

    def h(self, k, spins=None):
        """Bosonic LSWT matrix (Holstein-Primakoff, unit spin) around a collinear state; H = -sum J e.e."""
        s = self.spins if spins is None else np.asarray(spins, float)
        K0 = self.K([0, 0, 0])[0].real
        C = np.diag(-(K0*np.outer(s, s)).sum(1))
        out = []
        for Kq in self.K(k):
            A = -0.5*Kq*(1 + np.outer(s, s))
            B = -0.5*Kq*(1 - np.outer(s, s))
            out.append(2*np.block([[A - C, B], [B.conj().T, A - C]]))
        return np.array(out)


def grid(n, shift=0.0):
    g = (np.arange(n) + shift)/n
    return np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T


def lowest_h(model, n=8):
    return float(np.linalg.eigvalsh(model.h(grid(n)))[:, 0].min())


def kmin(model):
    """Wave vector (conventional r.l.u.) minimising the lowest eigenvalue of the FM reference's h(k)."""
    assert (model.M == np.eye(3)).all() and (model.spins > 0).all()
    lowest = lambda k: float(np.linalg.eigvalsh(model.h(k))[0, 0])
    ks = grid(16)
    k0 = ks[np.argmin(np.linalg.eigvalsh(model.h(ks))[:, 0])]
    best = min((minimize(lowest, k0 + d, method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 8000})
                for d in (np.zeros(3), np.array([0.01, 0, -0.01]), np.array([-0.01, 0.01, 0.01]))), key=lambda r: r.fun)
    return best.x - np.round(best.x), best.fun


def magnons(model, n=8):
    """Positive eigenvalues of g.h divided by S on an n^3 grid of the calculation's own Brillouin zone; refined maximum."""
    g = np.diag(np.r_[np.ones(model.n), -np.ones(model.n)])

    def omega(k):
        ev = np.linalg.eigvals(g @ model.h(k))
        assert np.abs(ev.imag).max() < 1e-6
        return np.sort(ev.real, axis=1)[:, model.n:]/S
    ks = grid(n, 0.5)
    w = omega(ks)
    i = np.unravel_index(np.argmax(w), w.shape)[0]
    r = minimize(lambda k: -omega(k).max(), ks[i], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-8})
    return float(w.mean()), float(max(w.max(), -r.fun))


def main():
    near = ['0-4', '4-6', '6-8', '8-10', '10-12', '12-14']
    # 1. FM references: instability and k_min, with an exchange-range convergence check
    fm_meta, fm_pairs = buy('FM/ligand', near + ['14-17'])
    fm14, fm17 = Model(fm_meta, fm_pairs, rmax=14.0), Model(fm_meta, fm_pairs)
    k_lig14, _ = kmin(fm14)
    k_lig, low_lig = kmin(fm17)
    k_lig6, _ = kmin(Model(fm_meta, fm_pairs, rmax=6.0))
    nl_meta, nl_pairs = buy('FM/no_ligand', near[:5])
    k_nol, low_nol = kmin(Model(nl_meta, nl_pairs))
    # 2. AFM reference with the ligand correction: its own (super)cell, moments and pair list
    af_meta, af_pairs = buy('AFM/ligand', near)
    af = Model(af_meta, af_pairs)
    low_af = lowest_h(af)
    K0 = af.K([0, 0, 0])[0].real
    w, v = np.linalg.eigh(K0)
    overlap = abs(v[:, -1] @ af.spins)/np.sqrt(af.n)                        # is the reference the LT optimum in its cell?
    G = np.array(np.meshgrid(range(-2, 3), range(-2, 3), range(-2, 3), indexing='ij')).reshape(3, -1).T
    qs = np.unique(np.round((G @ np.linalg.inv(af.M).T) % 1, 6), axis=0)
    # Fourier amplitude summed over the four conventional sublattices (each a Bravais lattice); the moments have a
    # single component. (Summing all sites together would cancel it: the two plaquette orientations are antiphase.)
    ideal = np.round(af.fconv*4)/4                                         # the AFM cell is strained by up to 0.4 %
    sub = np.round((ideal % 1)*4).astype(int) % 4
    groups = [np.where((sub == t).all(1))[0] for t in np.unique(sub, axis=0)]
    amp = [sum(abs(np.sum(af.spins[g]*np.exp(2j*np.pi*ideal[g] @ q))) for g in groups) for q in qs]
    assert sorted(amp)[-2] < 1e-6*max(amp)
    q_conv = qs[int(np.argmax(amp))]
    sym_split = {}
    for p in af_pairs:
        if 2.9 < p[9] < 3.0:
            i, j = af.names.index(p[0]), af.names.index(p[1])
            sym_split.setdefault('parallel' if af.spins[i]*af.spins[j] > 0 else 'antiparallel', []).append(p[5])
    split = {k: float(np.mean(x)) for k, x in sym_split.items()}
    chain_bonds = [p[5] for p in af_pairs if p[5] < 0.5*min(q[5] for q in af_pairs)]
    chain = -8*float(np.mean(chain_bonds))                                   # S = 1/2, each bond once, J > 0 AFM
    fm_chain = -8*float(np.mean([p[5] for p in fm_pairs if p[5] < 0.5*min(q[5] for q in fm_pairs)]))
    nl_chain = -8*float(np.mean([p[5] for p in nl_pairs if p[5] < 0.5*min(q[5] for q in nl_pairs)]))
    mean_w, max_w = magnons(af)
    # 3. AFM reference without the correction: first window is enough to test it
    an_meta, an_pairs = buy('AFM/no_ligand', ['0-4'])
    an = Model(an_meta, an_pairs)
    low_an = lowest_h(an)
    images = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing='ij')).reshape(3, -1).T @ an.cell
    oxy = np.array([o['position_A'] for o in an_meta['atoms'] if o['name'].startswith('O')])
    cu_o = sorted({round(float(d), 3) for c in an_meta['atoms'] if c['name'].startswith('Cu')
                   for d in np.linalg.norm(oxy[:, None, :] + images[None] - np.array(c['position_A']), axis=2).ravel()
                   if d < 2.1})
    # 4. validate the LSWT input against the authors' archived h(k) at Gamma (row 219 of their path)
    r = lab('evaluate', {'kind': 'author_dispersion', 'calculation': 'AFM/ligand'})
    evidence.append(r['id'])
    gamma = np.sort(r['result']['values'][219])
    mine = np.sort(np.linalg.eigvalsh(Model(af_meta, af_pairs).h([0, 0, 0])[0]))/1000
    gamma_dev = float(np.abs(gamma - mine).max())
    stable = {'FM/ligand': bool(low_lig > -1e-3), 'AFM/ligand': bool(low_af > -1e-3),
              'FM/no_ligand': bool(low_nol > -1e-3), 'AFM/no_ligand': bool(low_an > -1e-3)}
    submission = {'kmin_FM_reference_rlu': {'ligand_corrected': [float(x) for x in k_lig],
                                            'uncorrected': [float(x) for x in k_nol]},
                  'stable_at_own_reference': stable,
                  'converged_propagation_vector_rlu': [float(x) for x in q_conv],
                  'chain_exchange_meV': chain, 'mean_magnon_energy_meV': mean_w, 'magnon_band_maximum_meV': max_w,
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# CuO: self-consistent ground state from TB2J exchange (reference solution)

Conventions: archived TB2J energy E = -sum_(i!=j) J_ij e_i.e_j (ordered pairs, unit vectors, J>0 FM). LSWT uses the
Holstein-Primakoff bosonic matrix h(k) around each calculation's own collinear moments, in its own cell, with its own
pair list; magnon energies are the positive eigenvalues of g.h (paraunitary problem) divided by S = 1/2. The chain
coupling is quoted for H = sum_bonds J S_i.S_j (S = 1/2, each bond once): J = -8 J_TB2J.

FM references are unstable: lowest h eigenvalue {low_lig:.1f} meV (ligand-corrected) and {low_nol:.1f} meV (uncorrected).
k_min (ligand-corrected) = {np.round(k_lig, 4).tolist()} r.l.u. (14 A cutoff: {np.round(k_lig14, 4).tolist()}; 6 A cutoff:
{np.round(k_lig6, 4).tolist()}, so couplings beyond 8 A matter); the commensurate approximant is (1/2, 0, 1/2), i.e. the
2x1x2 AFM cell. Without the correction k_min = {np.round(k_nol, 4).tolist()}, a different (wrong) state.

AFM reference with correction: 2x1x2 supercell, moments with a single Fourier component at q = {q_conv.tolist()}
(AF1). Its own pair list makes h(k) positive semi-definite over the zone (lowest {low_af:.4f} meV, the Goldstone zero),
and the moment pattern is the top eigenvector of its K(0) (overlap {overlap:.5f}): the loop has converged. The archived
h(k) at Gamma agrees with mine to {gamma_dev:.1e}. Bonds that are equivalent in the crystal differ in this list with the
relative spin orientation (2.9-3.0 A: parallel {split.get('parallel', 0):.3f}, antiparallel {split.get('antiparallel', 0):.3f}
meV), a physical reference dependence, so the list is not symmetrised. Chain J = {chain:.1f} meV (FM reference
{fm_chain:.1f}, uncorrected {nl_chain:.1f} meV: the ligand correction lowers it by about a third).
Spin waves (S = 1/2): mean {mean_w:.1f} meV, band maximum {max_w:.1f} meV.

AFM reference without correction: two inequivalent Cu sites (Cu-O {cu_o} A), unequal moments, equivalent bonds split by tens of
meV, and h(k) negative ({low_an:.1f} meV) at its own reference: not converged.

The isotropic converged model selects commensurate collinear AF1. The multiferroic incommensurate AF2 spiral needs
anisotropic exchange/DMI and thermal fluctuations, which these isotropic data cannot address; a Neel temperature is not
supported (mean field overestimates it for this quasi-1D chain magnet).
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
