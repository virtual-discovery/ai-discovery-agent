"""Reference solution for hydride-coulomb-mustar; runs as the Harbor oracle and in check_task.

Route: the full-W Tc cannot be computed at the withheld pressures, so the Coulomb input is transferred.
mu = N(E_F) W(E_F,E_F); W(E_F,E_F) varies slowly with pressure while N(E_F) does not, so mu at a withheld
pressure is N(E_F) from its own DOS file times W(E_F,E_F) interpolated between released pressures. A fixed
mu* in the full-bandwidth vDOS solver at omega_c = 7 eV reproduces the full-W Tc when mu* = kappa * mu, with
kappa calibrated against full-W runs at a released pressure of the same compound. Tc(P) of H3S peaks near
180 GPa, so Tc is computed from each pressure's own alpha2F rather than interpolated.
"""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


def tc_list(system, p, solver, mus, wc=7000):
    res = evaluate({'kind': 'eliashberg', 'system': system, 'pressure_GPa': p, 'solver': solver,
                    'mu_star': [round(m, 5) for m in mus], 'omega_c_meV': wc})
    return res, [r['Tc_K'] for r in res['results']]


def invert(mus, tcs, target):
    """mu* at which Tc(mu*) = target, by interpolation of a monotonic decreasing table."""
    order = np.argsort(tcs)
    return float(np.interp(target, np.array(tcs)[order], np.array(mus)[order]))


def nef_at_ef(inputs):
    d = inputs['dos_file']
    return float(np.interp(d['EF_eV'], d['E_eV'], d['dos'])) / 2     # per spin, per cell


def main():
    # 1) full-W runs at released pressures: mu, W(E_F,E_F) and Tc
    W00, MU, TCW = {}, {}, {}
    for system, p in [('Mg2IrH6', 10), ('Mg2IrH6', 20), ('Mg2IrH6', 30), ('H3S', 200)]:
        r = evaluate({'kind': 'eliashberg', 'system': system, 'pressure_GPa': p, 'solver': 'vdos_W'})['results'][0]
        W00[(system, p)], MU[(system, p)], TCW[(system, p)] = r['W_EF_EF_eV'], r['mu_from_W'], r['Tc_K']
    # H3S 170 GPa has W but no alpha2F: W(E_F,E_F) from the kernel (bilinear at E_F, as the solver does)
    k = evaluate({'kind': 'w_kernel', 'system': 'H3S', 'pressure_GPa': 170})
    e, Wm = np.array(k['window_E_eV']), np.array(k['W_window_eV'])
    i = int(np.searchsorted(e, 0.0))
    t = (0 - e[i - 1]) / (e[i] - e[i - 1])
    W00[('H3S', 170)] = float((1 - t) ** 2 * Wm[i - 1, i - 1] + t * (1 - t) * (Wm[i - 1, i] + Wm[i, i - 1]) + t ** 2 * Wm[i, i])
    # 2) kappa = mu*_eff / mu from the full-W Tc at a released pressure of each compound
    kappa = {}
    for system, p, fr in [('Mg2IrH6', 20, (0.94, 0.97, 1.0)), ('H3S', 200, (0.86, 0.9, 0.94, 0.98))]:
        mus = [f * MU[(system, p)] for f in fr]
        _, tcs = tc_list(system, p, 'vdos', mus)
        kappa[system] = invert(mus, tcs, TCW[(system, p)]) / MU[(system, p)]
    # 3) withheld pressures: N(E_F) from their DOS files, W(E_F,E_F) interpolated in pressure
    neigh = {('Mg2IrH6', 15): (10, 20), ('Mg2IrH6', 25): (20, 30), ('H3S', 180): (170, 200), ('H3S', 190): (170, 200)}
    tc, inputs = {}, {}
    for (system, p), (a, b) in neigh.items():
        inputs[(system, p)] = evaluate({'kind': 'inputs', 'system': system, 'pressure_GPa': p})
        w = W00[(system, a)] + (W00[(system, b)] - W00[(system, a)]) * (p - a) / (b - a)
        mu = nef_at_ef(inputs[(system, p)]) * w
        _, tcs = tc_list(system, p, 'vdos', [kappa[system] * mu])
        tc[f'{system}_{p}GPa'] = tcs[0]
    # 4) equivalent constant-DOS mu* (omega_c = 1 eV) reproducing the archived full-W Tc at 20 GPa
    arch = evaluate({'kind': 'archived_run', 'system': 'Mg2IrH6', 'pressure_GPa': 20, 'run': 'vdos_W_1'})
    target = float(re.search(r'Tc = ([\d.]+)', arch['Summary.dat']).group(1))
    mus = [0.45, 0.55, 0.65, 0.75]
    _, tcs = tc_list('Mg2IrH6', 20, 'cdos', mus, wc=1000)
    m0 = invert(mus, tcs, target)
    mus2 = [m0 - 0.02, m0 - 0.007, m0 + 0.007, m0 + 0.02]
    _, tcs2 = tc_list('Mg2IrH6', 20, 'cdos', mus2, wc=1000)
    mu_eq = invert(mus2, tcs2, target)
    # 5) N(E_F) of H3S at 180 GPa per formula unit: the cell holds one S per formula unit
    scf = inputs[('H3S', 180)]['cell_scf_in']
    n_s = len(re.findall(r'^\s*S\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+', scf, re.M))
    nef_fu = nef_at_ef(inputs[('H3S', 180)]) / n_s
    submission = {'tc_full_W_K': tc, 'mu_star_cdos_1eV_Mg2IrH6_20GPa': mu_eq,
                  'N_EF_H3S_180GPa_per_fu_per_spin': nef_fu, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    mu_tab = ', '.join(f'{s} {p} GPa: mu = {MU[(s, p)]:.3f}, Tc_W = {TCW[(s, p)]:.1f} K' for s, p in MU)
    (APP/'report.md').write_text(f"""# First-principles Coulomb repulsion in Mg2IrH6 and H3S (reference solution)

mu = N(E_F) W(E_F,E_F) from the archived static RPA W: {mu_tab}. W(E_F,E_F) is nearly pressure
independent (Mg2IrH6 {W00[('Mg2IrH6', 10)]:.3f}-{W00[('Mg2IrH6', 30)]:.3f} eV; H3S {W00[('H3S', 170)]:.3f}-{W00[('H3S', 200)]:.3f} eV),
so mu tracks N(E_F): Mg2IrH6, with its large N(E_F), has mu ~ 0.55, about twice H3S.

mu is the bare Coulomb parameter at E_F, not a pseudopotential. In the full-bandwidth vDOS solver with a 7 eV
Matsubara cutoff a fixed mu* = kappa*mu reproduces the full-W Tc (kappa = {kappa['Mg2IrH6']:.3f} for Mg2IrH6,
{kappa['H3S']:.3f} for H3S, calibrated at 20 and 200 GPa). A constant-DOS calculation with a conventional 1 eV cutoff needs
mu* = {mu_eq:.3f} to reproduce the full-W Tc of Mg2IrH6 at 20 GPa ({target} K), far above 0.10-0.16; mu* = 0.16 overestimates Tc
by about 40 K. Morel-Anderson rescaling of mu to the cutoff gives a much smaller mu* and does not reproduce the full-W result.

Withheld pressures: mu from each pressure's own N(E_F) and interpolated W(E_F,E_F), Tc from its own alpha2F:
{json.dumps(tc)}. H3S Tc(P) peaks near 180 GPa (lambda grows as the lattice softens), so interpolating Tc between
160 and 200 GPa would miss by 5-10 %. The H3S cell holds {n_s} formula units: N(E_F) at 180 GPa is {nef_fu:.4f} states/eV/spin per H3S.

Caveats: static RPA screening, isotropic Eliashberg, harmonic alpha2F. Below ~170 GPa harmonic H3S phonons are unstable
(lambda becomes very large), and anharmonicity lowers Tc substantially; the 140-160 GPa harmonic Tc are not reliable.
Several archived run folders are mislabelled (e.g. Mg2IrH6 vdos_muStar0.1_1 used mu* = 0.2; H3S vdos_muStar0.2_1 used 0.1):
the Info.txt files, not the folder names, give the settings.
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
