"""Reference solution for inmnti2-weyl-census; runs as the Harbor oracle and in check_task."""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
ARC = APP/'archive'
B = np.array(CAT['conventions']['reciprocal_basis_B_invA'])
BINV = np.linalg.inv(B)
ALAT = 3.1335920710 * np.sqrt(2)          # |a1| of CELL_PARAMETERS (pw.x alat), Angstrom
TPIBA = 2 * np.pi / ALAT
SHIFTS = np.array([(i, j, k) for i in range(-2, 3) for j in range(-2, 3) for k in range(-2, 3)]) @ B
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    result = json.loads(out)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


def qe_bands(text):
    for marker in ('End of band structure calculation', 'End of self-consistent calculation'):
        if marker in text:
            text = text[text.rfind(marker):]
            break
    ks, es = [], []
    for b in re.split(r'\n\s+k =', text)[1:]:
        head, rest = b.split('\n', 1)
        ks.append([float(x) for x in re.findall(r'-?\d+\.\d+', head.replace('-', ' -'))[:3]])
        body = re.split(r'occupation numbers|the Fermi energy|highest occupied|Writing', rest)[0]
        es.append([float(x) for x in re.findall(r'-?\d+\.\d+', body.replace('-', ' -'))])
    n = min(len(e) for e in es)
    return np.array(ks), np.array([e[:n] for e in es])


def fermi(text):
    return float(re.findall(r'the Fermi energy is\s+(-?\d+\.\d+)', text)[-1])


def symmetry_ops(text):
    ops = []
    for b in re.split(r'\n\s+isym =\s*\d+', text.split('Sym. Ops.')[1].split('Cartesian axes')[0])[1:]:
        tr = re.search(r'Time Reversal\s+(\d)', b)
        rows = re.findall(r'\(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)', re.split(r'cart\.\s+s\(', b)[1])[:3]
        ops.append((np.array(rows, float), int(tr.group(1)) if tr else 0))
    return ops


def fold(k):
    cand = np.asarray(k) + SHIFTS
    return cand[np.argmin(np.linalg.norm(cand, axis=1))]


def same(k1, k2, tol):
    d = (np.asarray(k1) - np.asarray(k2)) @ BINV
    return np.linalg.norm((d - np.round(d)) @ B) < tol


def band_offset(name, scf_text):
    """DFT band index - Wannier band index, from the scf eigenvalues (DFT) and the Wannier bands at the same k."""
    k, e = qe_bands(scf_text)
    w = np.array(evaluate({'kind': 'wannier_bands', 'model': name, 'k_cart_invA': (k[:12] * TPIBA).tolist()})['energies_eV'])
    dev = {off: np.abs(w[:, 8:20] - e[:12, off + 8:off + 20]).max() for off in range(28, 41)}
    best = min(dev, key=dev.get)
    return best, dev[best]


def census(name, n, mesh, max_gap):
    scan = evaluate({'kind': 'gap_scan', 'model': name, 'lower_band': n, 'mesh': mesh})
    seeds = [m['k_cart_invA'] for m in scan['minima'] if m['gap_eV'] < max_gap]
    found = []
    for s in range(0, len(seeds), 8):
        for r in evaluate({'kind': 'find_nodes', 'model': name, 'lower_band': n, 'seeds_cart_invA': seeds[s:s + 8]})['results']:
            k = fold(r['k_cart_invA'])
            if not any(same(k, f['k'], 2e-3) for f in found):
                found.append({'k': k, 'gap': r['gap_eV'], 'E': r['mean_energy_eV']})
    nodes = [f for f in found if f['gap'] < 1e-6]
    return nodes, [f for f in found if f['gap'] >= 1e-6], scan


def orbits(nodes, ops, tol=0.05):
    parent = list(range(len(nodes)))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for i, f in enumerate(nodes):
        for R, t in ops:
            img = (-1 if t else 1) * R @ f['k']
            for j, g in enumerate(nodes):
                if same(img, g['k'], tol):
                    parent[root(i)] = root(j)
    groups = {}
    for i in range(len(nodes)):
        groups.setdefault(root(i), []).append(i)
    return list(groups.values())


def tpiba_b(text):
    lines = text.split('K_POINTS')[1].split('\n')
    pts = [(np.array([float(x) for x in lines[2 + i].split()[:3]]), int(lines[2 + i].split()[3])) for i in range(int(lines[1]))]
    ks = []
    for i, (p, n) in enumerate(pts):
        if i + 1 < len(pts):
            ks += [p + (pts[i + 1][0] - p) * j / n for j in range(n)]
        else:
            ks.append(p)
    n1, n3, n5 = pts[0][1], pts[2][1], pts[4][1]
    segs = [slice(0, n1 + 1), slice(n1 + 1, n1 + n3 + 2), slice(n1 + n3 + 2, n1 + n3 + n5 + 3)]
    return np.array(ks) * TPIBA, segs


def crossings(s, lo, hi, window=0.03, skip=0.004):
    g = hi - lo
    out = []
    for i0 in [i for i in range(1, len(g) - 1) if g[i] < g[i - 1] and g[i] <= g[i + 1] and g[i] < 0.01]:
        left = [i for i in range(len(s)) if skip < s[i0] - s[i] < window]
        right = [i for i in range(len(s)) if skip < s[i] - s[i0] < window]
        fit = lambda idx, y: np.polyfit(s[idx], y[idx], 1)
        a = (fit(left, hi) + fit(right, lo)) / 2
        b = (fit(left, lo) + fit(right, hi)) / 2
        sc = (b[1] - a[1]) / (a[0] - b[0])
        out.append((sc, np.polyval(a, sc), g[i0], (a[0], b[0])))
    return out


def main():
    scf = {m: (ARC/m/'scf/scf.out').read_text() for m in ('hubbard', 'standard')}
    EF = {m: fermi(t) for m, t in scf.items()}
    nocc = 52                                            # electrons per cell (scf.out), DFT bands 52/53
    off = {m: band_offset(m, scf[m]) for m in scf}       # hubbard 36 (34 excluded + 2 below dis_win_min), standard 34
    n_h, n_s = nocc - off['hubbard'][0], nocc - off['standard'][0]

    # DFT+U census in the Wannier model, chirality from Berry flux of bands 1..n_h
    hub, hub_other, hscan = census('hubbard', n_h, 24, 0.3)
    flux = []
    for s in range(0, len(hub), 8):
        flux += evaluate({'kind': 'berry_flux', 'model': 'hubbard', 'occupied_bands': n_h,
                          'spheres': [{'center_cart_invA': f['k'].tolist(), 'radius_invA': 0.01} for f in hub[s:s + 8]]})['results']
    for f, r in zip(hub, flux):
        f['chi'] = int(round(r['flux_over_2pi']))
    ops = symmetry_ops(scf['hubbard'])
    orb = orbits(hub, ops)
    fam = next(g for g in orb if all(abs(hub[i]['k'][2]) < 0.05 for i in g))
    fam_E = [1000 * (hub[i]['E'] - EF['hubbard']) for i in fam]

    # PBE (non-magnetic) census: 24 nodes; the ~7 meV minima near the pair midpoints are avoided crossings
    std, std_other, sscan = census('standard', n_s, 36, 0.05)
    pbe_E = 1000 * (np.mean([f['E'] for f in std]) - EF['standard'])

    # DFT truth on the archived DFT+U lines through one pair of the near-kz=0 family
    d = evaluate({'kind': 'dft_bands', 'run': 'hubbard/nscfXYZ'})
    kc, segs = tpiba_b(d['nscf.in'])
    e = np.array(d['eigenvalues_eV']) - EF['hubbard']
    res = {}
    for name, seg in zip('xyz', segs):
        ks = kc[seg]
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(ks, axis=0), axis=1))]
        res[name] = [(sc, E, g, sl, ks[0] + (ks[-1] - ks[0]) * sc / s[-1]) for sc, E, g, sl in
                     crossings(s, e[seg, nocc - 1], e[seg, nocc])]
    x1, x2 = sorted(res['x'], key=lambda c: c[0])
    E_dft = 1000 * np.mean([x1[1], res['y'][0][1], res['z'][0][1]])
    sep = x2[0] - x1[0]
    zc = res['z'][0][4]
    nearest = min(hub, key=lambda f: np.linalg.norm(f['k'] - zc))
    kz_plus = abs(zc[2]) * (1 if np.sign(zc[2]) == nearest['chi'] else -1)
    # The Wannier model along the same x-line, for comparison
    wx = np.array(evaluate({'kind': 'wannier_bands', 'model': 'hubbard', 'k_cart_invA': kc[segs[0]].tolist()})['energies_eV'])
    wgap = 1000 * (wx[:, n_h] - wx[:, n_h - 1]).min()

    # Kubo: folders are k-mesh sizes (berry_kmesh in kubo/InMnTi2.win, fixed 2 meV smearing); converge 210 vs 300
    comps = ['A_xy', 'A_yz', 'A_zx']
    k300 = evaluate({'kind': 'kubo', 'run': '300', 'components': comps})['rows']
    k210 = evaluate({'kind': 'kubo', 'run': '210', 'components': comps})['rows']
    xy300, xy210 = k300['A_xy'][0][1], k210['A_xy'][0][1]
    raw = {c: k300[c][0][1] for c in comps}

    submission = {'node_count': {'dftu': len(hub), 'pbe': len(std)}, 'orbit_sizes': sorted(len(g) for g in orb),
                  'node_energy_meV': float(E_dft), 'pair_separation_invA': float(sep), 'kz_offset_plus_invA': float(kz_plus),
                  'pbe_node_energy_meV': float(pbe_E), 'ahc_S_per_cm': {'xy': float(xy300), 'yz': 0.0, 'zx': 0.0},
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# InMnTi2 Weyl-node census (reference solution)

**Models.** `standard` is a non-magnetic PBE+SOC run (no starting magnetization, no magnetization in scf.out): its
nodes are time-reversal symmetric. `hubbard` (DFT+U with linear-response U on Mn and Ti) is the magnetic state:
Ti and Mn moments antiparallel along z, net 0.11 muB/cell; its scf keeps 8 operations, 4 of them combined with
time reversal. PBE and DFT+U total energies cannot be compared to choose the ground state.

**Band index.** The cell has 52 electrons. Matching Wannier and DFT eigenvalues at the scf k-points gives DFT band =
Wannier band + {off['hubbard'][0]} (hubbard; 34 excluded bands plus 2 below dis_win_min) and + {off['standard'][0]}
(standard), so the 52/53 crossings are Wannier bands {n_h}/{n_h + 1} and {n_s}/{n_s + 1}.

**Census.** Gap scans plus gap minimisation and Berry flux (bands 1..52) give {len(hub)} nodes in the DFT+U model
(chirality sum {sum(f['chi'] for f in hub)}) and {len(std)} in the PBE model. PBE minima with finite gaps
({sorted(round(1000 * f['gap'], 1) for f in std_other)} meV) are avoided crossings, not nodes; PBE pairs are only
0.155 1/A apart, so a coarse grid merges them into 12 minima. PBE nodes sit at {pbe_E:.1f} meV.

**Symmetry.** Under the magnetic group the 24 DFT+U nodes form orbits of sizes {sorted(len(g) for g in orb)}: one
orbit of 8 near kz = 0 (shifted to kz = +-delta, sign tied to chirality) and four orbits of 4 on the C2T-invariant
planes kx = 0 / ky = 0. Paper-style 'all equivalent' (24) holds only for the cubic non-magnetic group.

**Wannier fidelity.** The DFT+U Wannier model reproduces DFT to 0.4 meV on its 9x9x9 fitting grid (the scf k-points)
but not in between: nodes of one orbit sit at {min(fam_E):.0f} to {max(fam_E):.0f} meV and the model leaves a
{wgap:.1f} meV gap on the archived DFT line where DFT crosses (0.3 meV). Node energy, positions and separation
are therefore taken from the DFT lines: E_W = {E_dft:.1f} meV, pair separation {sep:.4f} 1/A, kz of the chi = +1
member {kz_plus:+.4f} 1/A. Velocities along the three lines ({x1[3][0]:.2f}/{x1[3][1]:.2f}, {res['y'][0][3][0]:.2f}/{res['y'][0][3][1]:.2f},
{res['z'][0][3][0]:.2f}/{res['z'][0][3][1]:.2f} eV A) have opposite signs: type I. Pauli blocking puts the
interband onset near 2|E_W|; the Gamma hole pocket and the Wannier node scatter smear the archived Kubo onset.

**Anomalous Hall.** The kubo folders are k-mesh sizes (berry_kmesh 300^3 in kubo/InMnTi2.win, 2 meV fixed
smearing), not broadenings. sigma_xy(0) = {xy300:.2f} (300^3) and {xy210:.2f} (210^3) S/cm. sigma_yz and sigma_zx
must vanish for m || z (C2z, C2x T); the archived {raw['A_yz']:.2f} and {raw['A_zx']:.2f} S/cm are numerical
symmetry breaking of the Wannier model.

**Not supported.** The electronic structure uses the PBE lattice (3.134 A primitive) where DFT+U is under about
238 kbar; the archived DFT+U vc-relax expands the cell to 3.364 A (+24% volume), so the census at the DFT+U
equilibrium volume is not established. Energies of the four kx/ky-plane orbits have no DFT check in the archive.
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
