"""Reference solution for gan-irradiation-defect-levels; runs as the Harbor oracle and in check_task.

Protocol: host energy from the perfect supercell computed like the defect cells (k = (1/2,1/2,1/2), same PAW
set); band edges from the Gamma-only perfect supercell with the same PAW set; eFNV (Kumagai-Oba) corrections
with the static dielectric tensor and core-potential alignment; lowest-energy configuration per charge state.
"""
import itertools
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.special import erfc

APP = Path(os.environ.get('APP_DIR', '/app'))
KE = 14.399645  # eV * Angstrom
NEUTRAL_NELECT = {'Va_Ga1': 381, 'Va_N1': 379, 'Ga_i1': 387, 'Ga_i2': 387, 'N_i1': 389, 'N_i2': 389,
                  'C_N1': 383, 'VV_GaN1': 376}
FAMILIES = {'V_Ga': ['Va_Ga1'], 'V_N': ['Va_N1'], 'Ga_i': ['Ga_i1', 'Ga_i2'], 'N_i': ['N_i1', 'N_i2'],
            'C_N': ['C_N1'], 'VV': ['VV_GaN1']}
# Purchase plan: each family's charge states until the next state's level is clearly outside the gap.
PLAN = {'Va_Ga1': [1, 0, -1, -2, -3], 'Va_N1': [3, 2, 1, 0, -1], 'Ga_i1': [3, 2, 1, 0], 'Ga_i2': [3, 2, 1, 0],
        'N_i1': [3, 2, 1, 0, -1, -2], 'N_i2': [3, 2, 1, 0, -1, -2], 'C_N1': [0, -1], 'VV_GaN1': [-1, -2, -3]}
EF_BELOW_CBM = 0.5


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    return r['id'], r['result']


# ------------------------------------------------------------------ geometry and electrostatics
def cart(rec, key):
    return np.array(rec[key]) @ np.array(rec['cell_A'])


def min_image(d, cell):
    f = d @ np.linalg.inv(cell)
    return (f - np.round(f)) @ cell


def map_to_perfect(rec, perfect):
    cell = np.array(perfect['cell_A'])
    pp, ps = cart(perfect, 'final_fractional_positions'), np.array(perfect['symbols'])
    ds = np.array(rec['symbols'])
    start = cart(rec, 'initial_fractional_positions')
    d = np.linalg.norm(min_image(start[:, None, :] - pp[None, :, :], cell), axis=2)
    pairs, used, extras, subs = {}, set(), [], []
    for i in np.argsort(d.min(axis=1)):
        free = np.where(np.isin(np.arange(len(pp)), list(used)), np.inf, d[i])
        j = int(np.argmin(free))
        if free[j] < 0.7:
            pairs[int(i)] = j
            used.add(j)
            if ds[i] != ps[j]:
                subs.append((int(i), j))
        else:
            extras.append(int(i))
    vac = [j for j in range(len(pp)) if j not in used]
    final = cart(rec, 'final_fractional_positions')
    # an atom displaced at the start that relaxed back onto the site it left is a host atom
    for i in list(extras):
        for j in list(vac):
            if ds[i] == ps[j] and np.linalg.norm(min_image(final[i] - pp[j], cell)) < 0.6:
                vac.remove(j)
                extras.remove(i)
                break
    pts = [pp[j] for j in vac] + [final[i] for i in extras] + [pp[j] for _, j in subs]
    pts = [pts[0] + min_image(p - pts[0], cell) for p in pts]
    return pairs, subs, np.mean(pts, axis=0)


class Ewald:
    """Periodic point charges with a neutralising background in an anisotropic dielectric."""
    def __init__(self, cell, eps, nreal=4, nrec=6):
        self.cell, self.eps = np.array(cell, float), np.array(eps, float)
        self.vol = abs(np.linalg.det(self.cell))
        self.einv, self.sdet = np.linalg.inv(self.eps), math.sqrt(np.linalg.det(self.eps))
        self.g = math.sqrt(np.trace(self.eps) / 3) / self.vol ** (1 / 3)
        rng = range(-nreal, nreal + 1)
        self.R = np.array(list(itertools.product(rng, rng, rng))) @ self.cell
        rec = 2 * np.pi * np.linalg.inv(self.cell).T
        rng = range(-nrec, nrec + 1)
        G = np.array([g for g in itertools.product(rng, rng, rng) if g != (0, 0, 0)]) @ rec
        GeG = np.einsum('ij,jk,ik->i', G, self.eps, G)
        self.G, self.wG = G, 4 * np.pi / self.vol * np.exp(-GeG / (4 * self.g ** 2)) / GeG

    def _real(self, r, exclude_origin):
        d = r[None, :] - self.R
        s = np.sqrt(np.einsum('ij,jk,ik->i', d, self.einv, d))
        if exclude_origin:
            s = s[s > 1e-8]
        return float(np.sum(erfc(self.g * s) / s)) / self.sdet

    def potential(self, r):
        return self._real(r, False) + float(np.sum(self.wG * np.cos(self.G @ r))) - math.pi / (self.vol * self.g ** 2)

    def madelung(self):
        return (self._real(np.zeros(3), True) + float(np.sum(self.wG)) - 2 * self.g / (math.sqrt(math.pi) * self.sdet)
                - math.pi / (self.vol * self.g ** 2))


def efnv(rec, perfect, q, ew):
    """Kumagai-Oba correction: PC energy plus alignment from core potentials outside the inscribed sphere."""
    if q == 0:
        return 0.0, 0.0
    cell = np.array(perfect['cell_A'])
    e_pc = -0.5 * q * q * KE * ew.madelung()
    pairs, subs, center = map_to_perfect(rec, perfect)
    vd, vp = np.array(rec['core_potentials_eV']), np.array(perfect['core_potentials_eV'])
    pos = cart(rec, 'final_fractional_positions')
    vol = abs(np.linalg.det(cell))
    r_in = min(vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)) / 2
    vals = []
    for i, j in pairs.items():
        if (i, j) in subs:
            continue
        r = min_image(pos[i] - center, cell)
        if np.linalg.norm(r) > r_in:
            # VASP core potentials are electron potential energies: phi_DFT = -(V_defect - V_perfect)
            vals.append((vd[i] - vp[j]) + q * KE * ew.potential(r))
    return e_pc + q * float(np.mean(vals)), e_pc


# ------------------------------------------------------------------ thermodynamics
def ground_state(E, vbm, ef):
    return min(((q, E[q] + q * (vbm + ef)) for q in E), key=lambda t: t[1])


def levels(E, vbm, gap):
    qs = sorted(E)
    grid = np.linspace(-1.0, gap + 1.0, 60001)
    gs = [qs[k] for k in np.array([E[q] + q * (vbm + grid) for q in qs]).argmin(axis=0)]
    out = []
    for k in range(1, len(gs)):
        a, b = gs[k - 1], gs[k]
        if a != b:
            x = (E[b] - E[a]) / (a - b) - vbm
            if 0 <= x <= gap:
                out.append({'from': a, 'to': b, 'level_eV': x})
    return out


def band_edges(rec):
    ev = rec['eigenvalues_eV_and_occupations']
    occ = [e for s in ev.values() for e, o in s if o > 0.5]
    emp = [e for s in ev.values() for e, o in s if o < 0.5]
    return max(occ), min(emp)


def main():
    evidence, recs = [], {}
    for name in ('perfect', 'perfect/kpt_gam'):
        i, recs[name] = buy({'kind': 'calculation', 'id': name})
        evidence.append(i)
    i, diel = buy({'kind': 'dielectric'})
    evidence.append(i)
    for fam, qs in PLAN.items():
        for q in qs:
            name = '%s_%d' % (fam, q)
            i, recs[name] = buy({'kind': 'calculation', 'id': name})
            evidence.append(i)
    perfect, gamma = recs['perfect'], recs['perfect/kpt_gam']
    kdef = recs['Va_Ga1_0']['kpoint_reciprocal_coordinates']
    assert perfect['kpoint_reciprocal_coordinates'] == kdef and perfect['potcar_titles'] == recs['Va_Ga1_0']['potcar_titles']
    assert gamma['kpoint_reciprocal_coordinates'] == [0.0, 0.0, 0.0] and gamma['potcar_titles'] == perfect['potcar_titles']
    vbm, cbm = band_edges(gamma)            # Gamma-point edges (the R-point edges of `perfect` are not band edges)
    gap = cbm - vbm
    r_vb, r_cb = band_edges(perfect)
    eps0 = np.array(diel['eps_0'])
    ew = Ewald(perfect['cell_A'], eps0)
    e_host = perfect['energy_sigma_to_0_eV']
    corr = {}
    fam_E = {}
    for fam, qs in PLAN.items():
        fam_E[fam] = {}
        for q in qs:
            rec = recs['%s_%d' % (fam, q)]
            assert NEUTRAL_NELECT[fam] - rec['nelect'] == q
            c, pc = efnv(rec, perfect, q, ew)
            corr['%s_%d' % (fam, q)] = (round(c, 3), round(pc, 3))
            fam_E[fam][q] = rec['energy_sigma_to_0_eV'] + c - e_host
    E = {}
    for d, fams in FAMILIES.items():
        qs = sorted(set().union(*[fam_E[f] for f in fams]))
        E[d] = {q: min(fam_E[f][q] for f in fams if q in fam_E[f]) for q in qs}
    lv = {d: levels(E[d], vbm, gap) for d in ('V_Ga', 'V_N', 'Ga_i', 'N_i', 'C_N')}
    ef = gap - EF_BELOW_CBM
    (q1, f1), (q2, f2), (q3, f3) = (ground_state(E[d], vbm, ef) for d in ('V_Ga', 'V_N', 'VV'))
    binding = f1 + f2 - f3
    (a, fa), (b, fb) = ground_state(E['V_N'], vbm, ef), ground_state(E['N_i'], vbm, ef)
    (c, fc), (dq, fd) = ground_state(E['V_Ga'], vbm, ef), ground_state(E['Ga_i'], vbm, ef)
    site_gap = {d: {q: round(fam_E[f1_][q] - fam_E[f2_][q], 3) for q in fam_E[f1_] if q in fam_E[f2_]}
                for d, (f1_, f2_) in {'Ga_i': ('Ga_i1', 'Ga_i2'), 'N_i': ('N_i1', 'N_i2')}.items()}
    submission = {'band_gap_eV': gap, 'transition_levels': lv, 'divacancy_binding_energy_eV': binding,
                  'frenkel_pair_energy_eV': {'N': fa + fb, 'Ga': fc + fd}, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    fmt = lambda d: ', '.join('(%+d/%+d) %.2f' % (l['from'], l['to'], l['level_eV']) for l in lv[d]) or 'none'
    (APP/'report.md').write_text(f"""# Native defects and C_N in GaN (reference solution)

Host and band edges. All defect cells were computed at the single k-point (1/2,1/2,1/2) with the Ga (3d in core)
PAW set; `perfect` uses the same k-point and set and is the host energy reference ({e_host:.3f} eV). At that
k-point the supercell edges are not the band edges (R-point gap {r_cb - r_vb:.2f} eV). The Gamma-only supercell with
the same PAW set (perfect/kpt_gam) gives VBM {vbm:.3f} eV and CBM {cbm:.3f} eV (gap {gap:.3f} eV); the Ga_d runs
(WAV_ref, kpt_gam/d10, 864 electrons) use another PAW set and their eigenvalue scale is incompatible. The
perfect-cell Fermi energy (Gaussian smearing) lies inside the gap and is not the VBM.

Finite size. eFNV (Kumagai-Oba): anisotropic point-charge energy with the static DFPT tensor eps_0 plus the
alignment from core potentials of atoms outside the inscribed sphere around the defect. For |q| = 3 the
corrections are 1.4-2.1 eV (point-charge part 1.7 eV); eps_inf would roughly double them and make the 2+
states of V_N and Ga_i stable. Corrections (total, PC) in eV: {corr}.

Configurations. Lowest-energy structure per charge state; Ga_i site i2 is lower by about 1 eV for q = +1 and +2,
N_i switches between the split (i2) and the i1 structure with charge (i1 - i2 in eV: {site_gap}).

Levels (eV above the VBM). V_Ga: {fmt('V_Ga')}. V_N: {fmt('V_N')}. Ga_i: {fmt('Ga_i')}. N_i: {fmt('N_i')}.
C_N: {fmt('C_N')}. V_N is negative-U ((3+/+): the 2+ state is never stable); for Ga_i the (3+/2+) and (2+/+)
levels coincide within a few meV, so the data cannot resolve a 2+ window. V_Ga is a triple deep
acceptor that compensates n-type material; N_i is amphoteric with deep levels; Ga_i is a donor up to the CBM;
C_N (0/-) is the shallowest acceptor level near 0.9 eV.

Pairs at E_F = E_CBM - 0.5 eV: V_Ga-V_N binding {binding:.2f} eV (states V_Ga {q1:+d}, V_N {q2:+d}, VV {q3:+d});
Frenkel pairs (isolated) N: {fa + fb:.2f} eV (V_N {a:+d}, N_i {b:+d}), Ga: {fc + fd:.2f} eV (V_Ga {c:+d}, Ga_i {dq:+d}).
Chemical potentials cancel in both.

Limits. One 96-atom supercell at one special k-point (delocalized carriers are penalised: R-point CB 2.3 eV above
the CBM), fixed low-spin moments (NUPDOWN 0/1; several even-electron states are broken-symmetry singlets),
no optical (vertical) data, mixing parameter 0.305 fixed, no elemental references (no absolute formation energies
or concentrations), and dielectric tensors from PBE rather than the hybrid functional.
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
