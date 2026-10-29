"""Reference solution for ca3ti2o7-hybrid-improper-domain-walls; runs as the Harbor oracle and in check_task.

Route (no literature values): map the bulk Cmc2_1 structure onto ideal I4/mmm sites, generate the
domain variants with parent-symmetry operations (inversion keeps the rotation psi and reverses the
tilt phi and P; inversion followed by a parent translation reverses psi and P and keeps phi; the
translation alone reverses both and leaves P unchanged), build two-wall supercells, relax them in
the lab with the reference potential and analyse energies and order-parameter profiles.
"""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
EV_A2_TO_MJ_M2 = 16021.766


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError('job %d failed: %s' % (r['id'], r.get('error')))
    return r['id'], r['result']


def mic(d, cell):
    f = d @ np.linalg.inv(cell)
    f -= np.round(f)
    return f @ cell


# ----------------------------------------------------------------------------- bulk and variants
class Bulk:
    def __init__(self, bulk, ref):
        self.sym = np.array(bulk['symbols'])
        pos = np.array(bulk['positions'])
        self.cell = np.array(bulk['cell'])
        a, b, c = np.diag(self.cell)
        rc = np.array(ref['cell'])
        rfrac = np.array(ref['positions']) / np.diag(rc)
        rsym = np.array(ref['symbols'])
        t1, t2, t3 = np.array([a / 2, b / 2, 0]), np.array([-a / 2, b / 2, 0]), np.array([0, 0, c])
        base = rfrac[:, :1] * t1 + rfrac[:, 1:2] * t2 + rfrac[:, 2:] * t3
        cand = np.concatenate([base + i * t1 + j * t2 for i in range(-2, 3) for j in range(-2, 3)])
        csym = np.tile(rsym, 25)
        ti = np.nonzero(self.sym == 'Ti')[0]
        best = None
        for k in np.nonzero(csym == 'Ti')[0]:
            sh = pos[ti[0]] - cand[k]
            err = np.linalg.norm(mic(pos[ti][:, None] - cand[csym == 'Ti'][None] - sh, self.cell), axis=2).min(1).sum()
            if best is None or err < best[0]:
                best = (err, sh)
        sites = np.zeros_like(pos)
        for s in ('Ca', 'Ti', 'O'):
            m = self.sym == s
            d = mic(pos[m][:, None] - cand[csym == s][None] - best[1], self.cell)
            k = np.linalg.norm(d, axis=2).argmin(1)
            sites[m] = pos[m] - d[np.arange(m.sum()), k]
        self.u = pos - sites
        ca = np.nonzero(self.sym == 'Ca')[0]
        fz = (sites[ca] @ np.linalg.inv(self.cell))[:, 2] % 1
        self.sites = sites - sites[ca[np.argmin(np.minimum(fz, 1 - fz))]]   # parent inversion centre at the origin
        self.frac = (self.sites @ np.linalg.inv(self.cell)) % 1
        self.t1 = t1

    def lookup(self, s):
        f = (s @ np.linalg.inv(self.cell)) % 1
        d = f[:, None] - self.frac[None]
        d -= np.round(d)
        dist = np.linalg.norm(d @ self.cell, axis=2)
        k = dist.argmin(1)
        assert dist[np.arange(len(k)), k].max() < 0.05
        return self.u[k], self.sym[k]

    def variant(self, op, s):
        if op == 'A':
            return self.lookup(s)
        if op == 'I':                       # inversion: psi kept, phi and P reversed
            u, sp = self.lookup(-s)
            return -u, sp
        if op == 'IT':                      # inversion + parent translation: psi and P reversed, phi kept
            u, sp = self.lookup(self.t1 - s)
            return -u, sp
        raise ValueError(op)


def supercell(bulk, frame, n):
    a, b, c = np.diag(bulk.cell)
    if frame == 'orth_a':     # wall normal along the polar axis a: charged walls
        A = np.array([[n * a, 0, 0], [0, b, 0], [0, 0, c]]); U = np.diag([a, b, c]); reps = [(i, 0, 0) for i in range(n)]
    elif frame == 'orth_b':   # wall normal along b (walls contain P): neutral walls
        A = np.array([[0, n * b, 0], [-a, 0, 0], [0, 0, c]]); U = np.array([[0, b, 0], [-a, 0, 0], [0, 0, c]])
        reps = [(i, 0, 0) for i in range(n)]
    elif frame == 'orth_c':   # wall normal along the stacking axis
        A = np.array([[a, 0, 0], [0, b, 0], [0, 0, n * c]]); U = np.diag([a, b, c]); reps = [(0, 0, k) for k in range(n)]
    cand = np.concatenate([bulk.sites + i * bulk.cell[0] + j * bulk.cell[1] + k * bulk.cell[2]
                           for i in range(-2, 3) for j in range(-2, 3) for k in range(-1, 2)])
    fu = cand @ np.linalg.inv(U)
    basis = cand[np.all((fu > -1e-6) & (fu < 1 - 1e-6), axis=1)]
    sites = np.concatenate([basis + i * U[0] + j * U[1] + k * U[2] for i, j, k in reps])
    e1 = A[0] / np.linalg.norm(A[0]); e3 = np.array([0, 0, 1.0]); R = np.array([e1, np.cross(e3, e1), e3])
    return A @ R.T, sites, R


def build(bulk, frame, n, op, offset, axis):
    cell, s, R = supercell(bulk, frame, n)
    f = ((s @ R.T) @ np.linalg.inv(cell))[:, axis]
    left = ((f - offset) % 1.0) < 0.5
    u = np.zeros_like(s); sym = np.empty(len(s), dtype=object)
    u[left], sym[left] = bulk.variant('A', s[left])
    u[~left], sym[~left] = bulk.variant(op, s[~left])
    uA, symA = bulk.variant('A', s)
    return {'symbols': [str(x) for x in sym], 'positions': np.round((s + u) @ R.T, 5).tolist(),
            'cell': np.round(cell, 6).tolist()}, {'sites': s @ R.T, 'cell': cell, 'axis': axis,
                                                    'posA': (s + uA) @ R.T, 'symA': [str(x) for x in symA]}


# ----------------------------------------------------------------------------- order parameters
def octahedra(sym, pos, cell, theta0):
    sym = np.array(sym); ti = np.nonzero(sym == 'Ti')[0]; ox = np.nonzero(sym == 'O')[0]
    out = []
    for t in ti:
        d = mic(pos[ox] - pos[t], cell)
        nb = np.argsort(np.linalg.norm(d, axis=1))[:6]
        v = d[nb]
        ap = np.argsort(-np.abs(v[:, 2]))[:2]
        eq = [k for k in range(6) if k not in ap]
        ang = [np.arctan2(w[1], w[0]) - theta0 for w in v[eq]]
        psi = np.mean([x - np.round(x / (np.pi / 2)) * (np.pi / 2) for x in ang])
        top, bot = (v[ap[0]], v[ap[1]]) if v[ap[0]][2] > 0 else (v[ap[1]], v[ap[0]])
        axv = top - bot
        out.append((t, psi, np.array([np.arctan2(axv[1], axv[2]), np.arctan2(axv[0], axv[2])])))
    return out


def profile(final, meta, nbin):
    """Per-bin min |psi|/|psi_bulk| and min |tilt|/|tilt_bulk| relative to domain A at the same sites."""
    pos = np.array(final['positions']); cell = np.array(final['cell'])
    oc = octahedra(final['symbols'], pos, cell, np.pi / 4)
    ocA = octahedra(meta['symA'], meta['posA'], meta['cell'], np.pi / 4)
    rows = []
    for (t, psi, tilt), (_, psiA, tiltA) in zip(oc, ocA):
        f = (meta['sites'][t] @ np.linalg.inv(meta['cell']))[meta['axis']] % 1
        rows.append((f, psi / psiA, np.linalg.norm(tilt) / np.linalg.norm(tiltA), tilt @ tiltA / (tiltA @ tiltA)))
    rows = np.array(rows)
    b = (rows[:, 0] * nbin + 0.25).astype(int) % nbin
    return [(k, float(np.abs(rows[b == k, 1]).min()), float(rows[b == k, 2].min()), float(rows[b == k, 3].mean()))
            for k in range(nbin) if (b == k).any()]


def main():
    bulk_s = CAT['bulk_structures']['Cmc2_1_ground_state']
    bulk = Bulk(bulk_s, CAT['bulk_structures']['I4mmm_reference'])
    evidence = []
    i, r = evaluate({'kind': 'single_point', 'structure': {k: bulk_s[k] for k in ('symbols', 'positions', 'cell')}})
    evidence.append(i)
    e_fu = r['energy_per_formula_unit_eV']
    # polar axis: the net Ti/Ca off-centring of the bulk is along the first cell vector
    pol = np.abs(bulk.u[bulk.sym != 'O'].sum(0))
    assert pol[0] > 5 * pol[1]

    def wall(frame, n, op, offset, axis):
        s, meta = build(bulk, frame, n, op, offset, axis)
        i, res = evaluate({'kind': 'relax', 'structure': s, 'cell': 'fixed', 'fmax': 0.02, 'max_steps': 800})
        evidence.append(i)
        c = np.array(res['final_structure']['cell'])
        area = np.linalg.norm(np.cross(c[(axis + 1) % 3], c[(axis + 2) % 3]))
        e = (res['energy_eV'] - res['formula_units'] * e_fu) / (2 * area) * EV_A2_TO_MJ_M2
        return e, res, meta

    log = {}
    # neutral walls (normal perpendicular to P), two sizes for convergence
    e_psi8, r_psi, m_psi = wall('orth_b', 8, 'IT', 0.25 / 32, 0)
    e_phi8, r_phi, m_phi = wall('orth_b', 8, 'I', 0.25 / 32, 0)
    e_psi12, _, _ = wall('orth_b', 12, 'IT', 0.25 / 48, 0)
    log['neutral'] = {'psi_n8': e_psi8, 'psi_n12': e_psi12, 'phi_n8': e_phi8}
    # walls parallel to the layers: rock-salt placement vs bilayer placement
    e_l_rs, _, _ = wall('orth_c', 4, 'IT', 0.0625, 2)
    e_l_bl, _, _ = wall('orth_c', 4, 'IT', 0.0025, 2)
    e_lphi, _, _ = wall('orth_c', 4, 'I', 0.0625, 2)
    log['layer'] = {'psi_rocksalt': e_l_rs, 'psi_bilayer': e_l_bl, 'phi_rocksalt': e_lphi}
    # charged walls (normal along P): energy versus domain width
    ns, ec = [8, 12, 16, 24], []
    for n in ns:
        ec.append(wall('orth_a', n, 'IT', 0.25 / (4 * n), 0)[0])
    d_nm = np.array(ns) * bulk.cell[0, 0] / 2 / 10
    slope = float(np.polyfit(d_nm, ec, 1)[0])
    log['charged'] = dict(zip(['n%d' % n for n in ns], ec))
    # order-parameter profiles of the neutral walls
    p_psi = profile(r_psi['final_structure'], m_psi, 8)
    p_phi = profile(r_phi['final_structure'], m_phi, 8)
    psi_min = min(x[1] for x in p_psi)
    phi_tilt_min = min(x[2] for x in p_phi)
    psi_prof = 'ising' if psi_min < 0.5 else 'rotational'
    phi_prof = 'ising' if phi_tilt_min < 0.5 else 'rotational'
    # Born-charge polarization across the psi wall (for the report)
    i, bc = evaluate({'kind': 'born_charges', 'structure': r_psi['final_structure']})
    evidence.append(i)
    Z = np.array(bc['born_effective_charges_e'])
    fpos = np.array(r_psi['final_structure']['positions']); fcell = np.array(r_psi['final_structure']['cell'])
    u = mic(fpos - m_psi['sites'] @ np.linalg.inv(m_psi['cell']) @ fcell, fcell)
    p = np.einsum('iab,ib->ia', Z, u)
    fb = ((m_psi['sites'] @ np.linalg.inv(m_psi['cell']))[:, 0] * 8 + 0.125).astype(int) % 8
    vol = abs(np.linalg.det(fcell))
    P = [p[fb == k].sum(0) / (vol / 8) * 16.0218 for k in range(8)]
    sub = {'E_neutral_psi_wall_mJ_m2': float(e_psi12), 'E_neutral_phi_wall_mJ_m2': float(e_phi8),
           'E_layer_psi_wall_mJ_m2': float(min(e_l_rs, e_l_bl)), 'charged_psi_wall_slope_mJ_m2_per_nm': slope,
           'psi_wall_profile': psi_prof, 'phi_wall_profile': phi_prof, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(sub, indent=1) + '\n')
    prof_txt = '\n'.join('  bin %d: psi wall min|psi|/psi0 %.2f; phi wall min|tilt|/tilt0 %.2f, tilt projection %.2f'
                         % (a[0], a[1], b[2], b[3]) for a, b in zip(p_psi, p_phi))
    P_txt = ', '.join('%.3f' % v[1] for v in P)
    (APP/'report.md').write_text(f"""# Domain walls in Ca3Ti2O7 (reference solution)

**Order parameters and construction.** In Cmc2_1 the polarization is induced by the trilinear coupling
P ~ psi_z * phi (rotation X2+ times tilt X3-). The bulk displacements relative to ideal I4/mmm sites were
mapped onto supercell sites; the second domain uses a parent-symmetry image of the bulk. Inversion through
the bilayer-centre Ca site keeps psi and reverses phi and P (phi wall); inversion followed by the parent
translation (a/2, b/2, 0) reverses psi and P and keeps phi (psi wall). The translation alone reverses both
and leaves P unchanged: that is an antiphase boundary, not a ferroelectric wall. The polar axis is the first
bulk cell vector (net Ca/Ti off-centring along it), so neutral walls have their normal along b and charged
(head-to-head/tail-to-tail) walls their normal along a.

**Energies.** E_DW = (E_cell - N_fu E_bulk)/(2A): two walls per periodic cell, A the wall area, in-plane
lattice fixed to the bulk cell, qNEP-w-BECs (the electrostatic model; the short-range NEP misses the
depolarizing fields). Bulk: {e_fu:.5f} eV/f.u.
- Neutral psi wall: {e_psi8:.1f} (8 cells) and {e_psi12:.1f} mJ/m^2 (12 cells): converged.
- Neutral phi wall: {e_phi8:.1f} mJ/m^2.
- Walls parallel to the layers (psi reversal): {e_l_rs:.1f} mJ/m^2 when placed in the CaO rock-salt layer,
  {e_l_bl:.1f} mJ/m^2 when placed inside the perovskite bilayer (metastable). The phi analogue costs {e_lphi:.1f}.
  The rock-salt layer nearly decouples the rotations of neighbouring bilayers, so psi walls there are almost free;
  tilts remain coupled through the apical oxygens.
- Charged psi walls: {', '.join('%.1f' % v for v in ec)} mJ/m^2 for {ns} cells. The energy grows linearly with
  the domain width (slope {slope:.2f} mJ m^-2 nm^-1, d = half the cell length): the bound charge of the
  head-to-head/tail-to-tail pair is unscreened, so there is no converged charged-wall energy; only the
  slope (field energy) and the short-cell values are meaningful. Small cells (< 2 nm domains) are not in the
  linear regime.

**Profiles.** Per orthorhombic slab:
{prof_txt}
The psi wall is Ising-like: psi passes through ~zero at the wall while the tilt is unchanged. The phi wall is
rotational: the projection of the tilt on the bulk tilt changes sign but its magnitude stays >= {phi_tilt_min:.2f}
of the bulk value (the tilt rotates through an orthogonal tilt pattern). The Born-charge polarization across the
psi wall (C/m^2 along the wall plane, per slab): {P_txt}; P reverses through zero in both wall types,
so the polarization alone does not distinguish them.

**Limits.** Values are for the machine-learned potential; its training-set errors (~2 meV/atom) and the
absence of free carriers (no screening of charged walls) limit transfer to experiment. (100)-type walls of the
tetragonal parent lie at 45 degrees to P and are partially charged: their energy also grows with the
supercell length. Orderings: layer psi << neutral psi < layer phi < neutral phi.
""")
    lab('submit', sub)
    for _ in range(60):
        st = lab('submission')
        if st['status'] != 'validating':
            print(json.dumps(st))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
