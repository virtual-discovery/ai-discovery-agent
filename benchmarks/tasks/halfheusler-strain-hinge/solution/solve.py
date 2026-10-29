"""Reference solution for halfheusler-strain-hinge; runs as the Harbor oracle and in check_task."""
import itertools
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
EVIDENCE = []


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
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error', 'job failed'))
    EVIDENCE.append(r['id'])
    return r['result']


def lattice(name):
    return np.array(CAT['models'][name]['lattice_vectors_A'])


def recip(name):
    return 2 * np.pi * np.linalg.inv(lattice(name)).T


# ---------------------------------------------------------------- geometry: strain and the prism basis of each file
def strain(name):
    """delta = (c0 - c)/c0 with c along z: a primitive fcc/bct cell has |z components| = c/2."""
    a0 = 2 * np.abs(lattice(name.split('_')[0] + '_orig')).max()
    c = 2 * np.abs(lattice(name)[:, 2]).max()
    return 100 * (a0 - c) / a0


def prism_basis(name):
    """Integer basis (rows in the file's lattice vectors) of the physical prism P_N, which is defined in the orig cubic frame
    with the orig identity basis.  Try rotations of the orig frame about z; keep the one that maps the orig lattice vectors
    onto integer combinations of this file's vectors and the Zn site onto reduced (1/4, 1/4, 1/4) from X."""
    ref = lattice(name.split('_')[0] + '_orig')
    A = lattice(name)
    cen = np.array(CAT['models'][name]['wannier_centres_A'])
    dzn = cen[0] - cen[2]
    for deg in range(0, 360, 15):
        t = np.radians(deg)
        rot = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
        S = (ref @ rot.T) @ np.linalg.inv(A)
        if not np.allclose(S, np.round(S), atol=1e-4):
            continue
        S = np.round(S).astype(int)
        f = dzn @ np.linalg.inv(S @ A)
        if abs(round(np.linalg.det(S))) == 1 and np.allclose(f - np.floor(f), 0.25, atol=1e-4):
            return S.tolist()
    raise RuntimeError('no basis found')


# ---------------------------------------------------------------- bulk
def gamma_levels(name):
    r = evaluate({'kind': 'bulk_bands', 'model': name, 'k_frac': [[0, 0, 0]]})
    return np.array(r['energies_eV'][0]), np.array(r['zn_s_weight'][0])


def edges(name, band=6, refine_max=False):
    r = evaluate({'kind': 'band_edges', 'model': name, 'band': band, 'mesh': 24})
    vbm = r['band_n_highest'][0]
    cbm = r['band_n_plus_1_lowest'][0]
    # refine the band-(n+1) minimum on a local grid around the grid minimum
    k0 = np.array(cbm['k_frac'])
    d = np.linspace(-0.02, 0.02, 6)
    cloud = [list(k0 + np.array(x)) for x in itertools.product(d, d, d)]
    e = np.array(evaluate({'kind': 'bulk_bands', 'model': name, 'k_frac': cloud})['energies_eV'])
    cbm_e = min(cbm['E_eV'], e[:, band].min())
    return vbm['E_eV'], cbm_e


# ---------------------------------------------------------------- prism
def prism(name, basis, n, ks, window, potential=None):
    req = {'kind': 'prism', 'model': name, 'basis': basis, 'size': n, 'k_phases': list(ks), 'window_eV': window}
    if potential:
        req['potential'] = potential
    return evaluate(req)


def hinge_states(res, corner=None, wmin=0.3):
    out = []
    for r in res['results']:
        hs = []
        for s in r['states']:
            c = corner or max(('c00', 'c0N', 'cN0', 'cNN'), key=lambda x: s['weights'][x])
            if s['weights'][c] > wmin:
                hs.append((s['E_eV'], c, s['zn_s_weight']))
        out.append(hs)
    return out


def main():
    names = list(CAT['models'])
    strains = {n: strain(n) for n in names}
    g_sb, z_sb = gamma_levels('LiSbZn_orig')
    g_bi, z_bi = gamma_levels('LiBiZn_orig')
    # the Zn-4s doublet at Gamma: above the Gamma8 quartet in LiSbZn, below it in LiBiZn
    inverted = {'LiSbZn': bool(np.argmax(z_sb) < 4), 'LiBiZn': bool(np.argmax(z_bi) < 4)}
    gamma8_bi = float(np.median(g_bi[4:8]))
    # six electrons per cell (Li 1 + Zn 2 + X-p 3): bands 1-6 occupied
    vbm, cbm = {}, {}
    for n in ('LiSbZn_orig', 'LiSbZn_strained_setup', 'LiBiZn_orig'):
        vbm[n], cbm[n] = edges(n)
    gaps = {n: 1000 * (cbm[n] - vbm[n]) for n in vbm}
    ef = evaluate({'kind': 'fermi_level', 'model': 'LiBiZn_orig', 'electrons_per_cell': 6, 'mesh': 64, 'offset': 0.5})
    ef_rel = 1000 * (ef['fermi_level_eV'] - gamma8_bi)

    # ---- LiSbZn prism: hinge state at k = 0 and k = pi, size convergence
    name = 'LiSbZn_orig'
    I3 = prism_basis(name)
    v0, c0 = vbm[name], cbm[name]
    win = [v0 + 0.02, c0]
    e12 = hinge_states(prism(name, I3, 12, [0.0], win))[0]
    r20 = prism(name, I3, 20, [0.0, np.pi], win)
    e20, e20pi = hinge_states(r20)
    r28 = prism(name, I3, 28, [0.0], win)
    e28 = hinge_states(r28)[0]
    corner = e28[0][1]
    levels = [e12[0][0], e20[0][0], e28[0][0]]
    q = (levels[2] - levels[1]) / (levels[1] - levels[0])
    e_inf = levels[2] + (levels[2] - levels[1]) * q / (1 - q) if 0 < q < 1 else levels[2]
    # outermost species of each facet from the cell geometry: Zn at reduced +f from X along a side is outermost at the high
    # end (t = N-1) when f > 0 and at the low end (t = 0) when f < 0 (f taken in (-1/2, 1/2))
    zr = np.array(r28['geometry']['zn_site_reduced'])[1:]
    zr = zr - np.round(zr)
    ends = [('N' if f > 0 else '0') for f in zr]                 # end of each side where Zn is outermost
    species = ['Zn' if corner[1 + i] == ends[i] else 'X' for i in range(2)]
    termination = '/'.join(sorted(species, key=lambda x: x != 'Zn'))
    zn_w = e28[0][2]
    # strained setup: the same physical prism in that file's own basis
    sname = 'LiSbZn_strained_setup'
    Ss = prism_basis(sname)
    es = hinge_states(prism(sname, Ss, 28, [0.0], [vbm[sname] + 0.02, cbm[sname]]))[0]
    energy_k0 = {name: 1000 * (e_inf - v0), sname: 1000 * (es[0][0] - vbm[sname])}

    # ---- hinge band over the 1D BZ (N = 12: the band is strongly localized away from k = 0)
    ks = np.linspace(0, np.pi, 16)
    band = []
    for chunk in (ks[:8], ks[8:]):
        res = prism(name, I3, 12, chunk, [0.15, 0.80])
        for k, hs in zip(chunk, hinge_states(res, corner, 0.5)):
            band.append((k, [h[0] for h in hs[:2]]))
    tops = [(max(e), k) for k, e in band if e]
    kmax = max(tops)[1]
    fine = np.linspace(kmax - 0.15, kmax + 0.15, 8)
    res = prism(name, I3, 12, fine, [0.15, 0.80])
    top = max([max(h[0] for h in hs[:2]) for hs in hinge_states(res, corner, 0.5) if hs] + [max(tops)[0]])
    bottom = min(e20pi[0][0], min(min(e) for k, e in band if e))
    window = [1000 * (bottom - v0), 1000 * (top - v0)]

    # ---- hinge potential: attractive V on the 2x2 cells of that corner pushes the level to the VBM
    vs, lv = [-0.30, -0.38, -0.46], []
    for V in vs:
        pot = [{'region': corner, 'depth': 2, 'V_eV': V, 'species': 'all'}]
        hs = hinge_states(prism(name, I3, 12, [0.0], [v0 - 0.08, c0], pot), corner)[0]
        lv.append(hs[0][0])
    coef = np.polyfit(vs, np.array(lv) - v0, 2)
    roots = [r.real for r in np.roots(coef) if abs(r.imag) < 1e-9 and -0.6 < r.real < -0.2]
    removal = float(roots[0]) if roots else None

    sub = {'strain_percent': {n: round(strains[n], 4) for n in names},
           'band_inverted': inverted,
           'fundamental_gap_meV': gaps,
           'libizn_fermi_level_meV': ef_rel,
           'hinge': {'termination': termination, 'energy_k0_meV': energy_k0, 'band_window_meV': window,
                     'zn_s_weight_k0': zn_w, 'removal_potential_eV': removal},
           'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(sub, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# LiSbZn and LiBiZn: hinge states (reference solution)

**Archived models.** The two `strained_setup` files are the equilibrium cubic structures written as body-centred
tetragonal cells (c = sqrt(2) a_t = c0, same volume): strain {strains}. Their Cartesian frame is rotated by 45 degrees
about z and their lattice basis differs, so the prism P_N is built there with basis {Ss} (orig: {I3}); the Zn site then
sits at reduced (1/4, 1/4, 1/4) as in the orig cell. The FPLO energy zero lies at an arbitrary point of the LiSbZn gap
(VBM at {1000 * v0:.1f} meV in orig and {1000 * vbm[sname]:.1f} meV in the strained setup), so LiSbZn energies are
referred to the VBM. The two LiSbZn calculations differ (gaps {gaps['LiSbZn_orig']:.1f} and
{gaps['LiSbZn_strained_setup']:.1f} meV) although their geometry is identical: a calculation difference, not strain.
The archive therefore contains no strained Hamiltonian and the strain dependence cannot be tested.

**Bulk.** Zn-4s weight at Gamma: LiSbZn normal order (Gamma6 above Gamma8), LiBiZn inverted (Gamma6 well below Gamma8).
Six electrons per cell (Li 1, Zn 2, X-p 3; Li has no Wannier function) fill bands 1-6. The LiSbZn conduction-band minimum
lies next to L (full-zone search), giving an indirect gap. LiBiZn is a compensated semimetal: band 7 near L lies
{-gaps['LiBiZn_orig']:.0f} meV below the Gamma8 maximum; the charge-neutral E_F (64^3 mesh) is {ef_rel:.0f} meV below
Gamma8, i.e. at the FPLO zero, not at the Gamma8 touching point.

**Hinge state.** In P_N of LiSbZn one Kramers doublet lies in the gap at k = 0, localized at corner {corner}: with the Zn
site at +1/4 along both sides the X sites are outermost on the t = 0 facets and Zn on the t = N-1 facets, so the hinge is
{termination} (two Zn-terminated {{111}} facets), with Zn-4s weight {zn_w:.2f} (cation dangling-bond character). This
contradicts the paper's Bi/Sb-hinge assignment; the authors' kwant HoppingKind assignment reverses the hopping vectors
relative to the site positions, which mirrors their real-space picture without changing the spectrum. Size convergence at
k = 0: N = 12, 20, 28 give {levels}; the converged level is {energy_k0[name]:.1f} meV above the VBM
({energy_k0[sname]:.1f} meV in the strained-setup model). The hinge band has its minimum at k = pi, two spin-split
branches peaking near k = {kmax:.2f} at up to {window[1]:.0f} meV, inside the projected gap but far above the global CBM,
and never touching the valence band.

**Origin.** The band is detached from the valence band, appears at one of the four hinges only (termination dependent), and
an attractive potential of {removal} eV on the 2x2 hinge cells pushes it into the valence band while new bound states are
pulled out of the conduction band. With a normal-ordered bulk and gapped facets this is a trivial hinge state, not a
higher-order topological one. In LiBiZn, compensated bulk pockets and surface states at the Fermi level hybridize with
hinge-weighted states, so hinge modes are far less visible than in LiSbZn.

**Not supported by these data:** any strain dependence (no strained Hamiltonian is archived), surface reconstructions,
and the fidelity of the 8-band Wannier models to DFT (no DFT bands are archived).
""")
    lab('submit', sub)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
