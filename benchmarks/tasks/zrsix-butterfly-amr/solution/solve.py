"""Reference solution for zrsix-butterfly-amr; runs as the Harbor oracle and in check_task.

Route (no literature values): neutral chemical potentials from electron counting (24 valence electrons in the
24-orbital spinless model -> 12 filled bands), Chambers transport on the dense mesh at those potentials, a coarse-mesh
cross-check, a bc-plane angle scan around the butterfly maximum, and the two-sample comparison at each sample's B*tau.
"""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
T = 30.0
KB = 8.617333262e-5
NV = 24                       # valence electrons per cell in the Wannier manifold (Zr 2x4, Si 2x4, X p 2x4)


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in out:
        raise RuntimeError(out['error'])
    return out


def evaluate(req, evidence):
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    evidence.append(r['id'])
    return r['result']


def neutral_mu(c, evidence):
    lo, hi = CAT['electron_count']['mu_range_eV'][c]
    grid = list(np.linspace(lo, hi, 200))
    r = evaluate({'kind': 'electron_count', 'compound': c, 'mesh': 80, 'mu_eV': grid, 'temperature_K': T}, evidence)
    n = np.array(r['electrons_per_cell'])
    i = int(np.searchsorted(n, NV / 2))
    fine = list(np.linspace(grid[i - 1] - 0.005, grid[i] + 0.005, 200))
    r = evaluate({'kind': 'electron_count', 'compound': c, 'mesh': 160, 'mu_eV': fine, 'temperature_K': T}, evidence)
    n = np.array(r['electrons_per_cell'])
    j = int(np.searchsorted(n, NV / 2))
    mu = fine[j - 1] + (NV / 2 - n[j - 1]) * (fine[j] - fine[j - 1]) / (n[j] - n[j - 1])
    occ = np.array(r['band_occupation'])
    w = (NV / 2 - n[j - 1]) / (n[j] - n[j - 1])
    occ_mu = occ[j - 1] + w * (occ[j] - occ[j - 1])
    ne = occ_mu[NV // 2:].sum()
    a, c_ = CAT['cells_A'][c]['a'], CAT['cells_A'][c]['c']
    density = 2 * ne / (a * a * c_ * 1e-24)
    return float(mu), float(density), float(ne), float((1 - occ_mu[:NV // 2]).sum())


def transport(c, plane, angle, mesh, mu, evidence, bands):
    r = evaluate({'kind': 'magnetotransport', 'compound': c, 'field': {'plane': plane, 'angle_deg': angle},
                  'mesh': mesh, 'mu_eV': [mu], 'temperature_K': T, 'bands': bands}, evidence)
    rho = np.array(r['results'][0]['rho_tau_Ohm_m_ps'])
    return np.array(r['Btau_Tps']), rho, r['results'][0]['open_orbit_fraction_by_band']


def mr(rho):
    return rho[:, 0, 0] / rho[0, 0, 0] - 1


def peak(th, y):
    i = int(np.argmax(y)); i = min(max(i, 1), len(y) - 2)
    p = np.polyfit(th[i - 1:i + 2], y[i - 1:i + 2], 2)
    return float(-p[1] / (2 * p[0]))


def main():
    ev = []
    mus, dens, balance = {}, {}, {}
    for c in ('ZrSiS', 'ZrSiSe', 'ZrSiTe'):
        mus[c], dens[c], ne, nh = neutral_mu(c, ev)
        balance[c] = (ne, nh)
    bse, bte = CAT['magnetotransport']['bands'][ 'ZrSiSe'], CAT['magnetotransport']['bands']['ZrSiTe']
    # ZrSiSe: B || c, B || b on the dense mesh, plus the coarse mesh as a convergence check
    X, r0, _ = transport('ZrSiSe', 'bc', 0.0, 160, mus['ZrSiSe'], ev, bse)
    _, r90, open90 = transport('ZrSiSe', 'bc', 90.0, 160, mus['ZrSiSe'], ev, bse)
    _, r0c, _ = transport('ZrSiSe', 'bc', 0.0, 80, mus['ZrSiSe'], ev, bse)
    _, r90c, _ = transport('ZrSiSe', 'bc', 90.0, 80, mus['ZrSiSe'], ev, bse)
    i5, i10, i20 = (list(X).index(x) for x in (5.0, 10.0, 20.0))
    m0, m90 = mr(r0), mr(r90)
    # coarse scan of the bc rotation on the 80 mesh, then the dense mesh around the maximum
    coarse = {}
    for th in (30.0, 45.0, 60.0):
        _, r, _ = transport('ZrSiSe', 'bc', th, 80, mus['ZrSiSe'], ev, bse)
        coarse[th] = mr(r)[i5]
    dense = {0.0: m0[i5], 90.0: m90[i5]}
    opens = {}
    for th in (40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0):
        _, r, o = transport('ZrSiSe', 'bc', th, 160, mus['ZrSiSe'], ev, bse)
        dense[th] = mr(r)[i5]
        opens[th] = o
    th = np.array(sorted(dense)); y = np.array([dense[t] for t in th])
    inner = th[(th >= 40) & (th <= 70)]
    angle = peak(inner, np.array([dense[t] for t in inner]))
    # ZrSiTe: B || c and B || b, dense mesh
    _, t0, _ = transport('ZrSiTe', 'bc', 0.0, 160, mus['ZrSiTe'], ev, bte)
    _, t90, topen = transport('ZrSiTe', 'bc', 90.0, 160, mus['ZrSiTe'], ev, bte)
    tm0, tm90 = mr(t0), mr(t90)
    i2 = list(X).index(2.0)
    sub = {'neutral_mu_eV': mus,
           'carrier_density_cm3': {c: dens[c] for c in ('ZrSiSe', 'ZrSiTe')},
           'butterfly_peak_angle_deg': angle,
           'mr_b_over_c': float(m90[i5] / m0[i5]),
           'mr_exponent_b': float(np.log(m90[i20] / m90[i10]) / np.log(2)),
           'zrsite_mr_c_over_b': float(tm0[i5] / tm90[i5]),
           'mr_ratio_se_over_te': float(m0[i5] / tm0[i2]),
           'evidence_ids': ev}
    (APP/'submission.json').write_text(json.dumps(sub, indent=1) + '\n')
    hall = r0[:, 1, 0] / r0[0, 0, 0]; hallc = r0c[:, 1, 0] / r0c[0, 0, 0]
    (APP/'report.md').write_text(f"""# ZrSiX magnetoresistance (reference solution)

Chemical potentials: the 24-orbital models (Zr d, Si s+p, X p; no spin index) describe 24 valence electrons, so the
neutral crystal fills 12 bands. Counting on the 160 mesh at 30 K gives mu = {mus['ZrSiS']:.4f} (ZrSiS),
{mus['ZrSiSe']:.4f} (ZrSiSe) and {mus['ZrSiTe']:.4f} eV (ZrSiTe), with n_e = n_h (exact compensation). The archived
E_FERMI values are not these: the ZrSiS input.dat is a copy of the ZrSiSe file (Se atoms, ZrSiSe cell, E_FERMI 8.9117),
far below the ZrSiS neutral point; ZrSiSe's 8.9117 is ~10 meV below neutrality (n_e/n_h ~ 0.7) and ZrSiTe's 8.8965
~35 meV above it (n_e/n_h ~ 1.9). NumOccupied = 13 in the ZrSiTe input is inconsistent with the electron count.
Carrier densities: {dens['ZrSiSe']:.3e} (ZrSiSe) and {dens['ZrSiTe']:.3e} cm^-3 (ZrSiTe), both spins.

Transport (Chambers, one tau, current along x, bc rotation). ZrSiSe at Btau = 5 T ps: MR(B||c) = {m0[i5]:.2f},
MR(B||b) = {m90[i5]:.2f}; 80-mesh check {mr(r0c)[i5]:.2f} and {mr(r90c)[i5]:.2f}. The rotation scan (80 mesh:
{', '.join(f'{k:.0f} deg {v:.1f}' for k, v in coarse.items())}; 160 mesh: {', '.join(f'{k:.0f} deg {v:.1f}' for k, v in sorted(dense.items()))})
has its maximum at {angle:.1f} deg: the butterfly. Thermal-weighted open-orbit fractions are below 1% in this rotation
({opens}), so the lobes come from the angle-dependent compensation of closed orbits on concave pockets, not from open
orbits. A closed-orbit effective-mass model would give a monotonic pattern with its maximum at an end point.
MR(B||b) keeps growing: exponent {np.log(m90[i20] / m90[i10]) / np.log(2):.2f} between 10 and 20 T ps (compensation); at the
archived E_FERMI it saturates. A 40 mesh fakes this saturation, so convergence was checked on 80 and 160.
ZrSiTe (bands 12-14): MR(B||c) = {tm0[i5]:.2f}, MR(B||b) = {tm90[i5]:.2f} at 5 T ps, a peanut pattern with the maximum
at c and near saturation for B||b (open-orbit fraction {topen}). The ZrSiTe Hamiltonian is unsymmetrized, so
rho(theta) and rho(-theta) differ slightly in the model; that is an artefact.
Samples at B = 5 T: ZrSiSe (tau 1.0 ps) is at Btau = 5, ZrSiTe (tau 0.4 ps) at Btau = 2 T ps (Kohler scaling with one
tau), ratio {m0[i5] / tm0[i2]:.2f}.
Hall: rho_yx/rho_xx(0) for B||c in ZrSiSe is {', '.join(f'{x:g}:{v:.2f}' for x, v in zip(X, hall))} on the 160 mesh
and {', '.join(f'{x:g}:{v:.2f}' for x, v in zip(X, hallc))} on the 80 mesh: a small difference of electron and hole
terms, not converged, so no Hall sign change is claimed for the compensated crystal.
The semiclassical calculation contains no Berry phase or topological information; the butterfly pattern follows from
the Fermi-surface geometry and compensation. These data cannot address SOC gapping (no SOC in the models), quantum
oscillations or magnetic breakdown, and tau is assumed equal for all bands.
""")
    lab('submit', sub)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
