"""Reference solution for mos2-doping-dome-cdw; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import subprocess
import time
from pathlib import Path

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 8.61733e-5
MU = 0.13
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


def allen_dynes(lam, wlog_meV, w2_meV=None):
    wlog = wlog_meV * 1e-3
    tc = wlog / (1.2 * KB) * math.exp(-1.04 * (1 + lam) / (lam - MU * (1 + 0.62 * lam)))
    tc *= (1 + (lam / (2.46 * (1 + 3.8 * MU))) ** 1.5) ** (1 / 3)
    if w2_meV is not None:
        r = w2_meV / wlog_meV
        tc *= 1 + (r - 1) * lam ** 2 / (lam ** 2 + (1.82 * (1 + 6.3 * MU) * r) ** 2)
    return tc


def lerp(x, x0, y0, x1, y1):
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def main():
    # 1. Hidden-doping Tc: Eliashberg (full DOS, FA) at the neighbours, transferred with the authors' lambda/omega_log
    #    table, because lambda jumps between 0.14 and 0.16 (Q-valley filling and the M-point soft mode).
    mom = evaluate({'kind': 'moments_table'})
    lam_fa = {'%.2f' % r[0]: r[3] for r in mom['lambda']['rows']}
    wlog_fa = {'%.2f' % r[0]: r[3] for r in mom['omega_log_meV']['rows']}
    el = {x: evaluate({'kind': 'eliashberg', 'phase': '1x1H', 'label': x, 'flavor': 'FA', 'dos': 'full', 'mu_star': MU})
          for x in ('0.14', '0.16', '0.18')}
    tc_hidden, notes = {}, {}
    for x, (lo, hi) in {'0.15': ('0.14', '0.16'), '0.17': ('0.16', '0.18')}.items():
        ratio = [el[k]['Tc_K'] / allen_dynes(el[k]['lambda'], el[k]['omega_log_meV'], el[k]['omega_2_meV']) for k in (lo, hi)]
        w2 = (el[lo]['omega_2_meV'] + el[hi]['omega_2_meV']) / 2
        ad = allen_dynes(lam_fa[x], wlog_fa[x], w2)
        tc_hidden[x] = ad * sum(ratio) / 2
        notes[x] = (ad, ratio, (el[lo]['Tc_K'] + el[hi]['Tc_K']) / 2)
    # Cross-check at 0.15: neighbours' spectra averaged, rescaled to the tabulated lambda, rigid-band DOS of 0.14.
    s0 = evaluate({'kind': 'spectrum', 'phase': '1x1H', 'label': '0.14', 'flavor': 'FA'})['rows']
    s1 = evaluate({'kind': 'spectrum', 'phase': '1x1H', 'label': '0.16', 'flavor': 'FA'})['rows']
    rows = [[a[0], (a[1] + b[1]) / 2] for a, b in zip(s0, s1) if a[0] > 0]
    lam_avg = 2 * sum((rows[i + 1][0] - rows[i][0]) / 2 * (rows[i][1] / rows[i][0] + rows[i + 1][1] / rows[i + 1][0])
                      for i in range(len(rows) - 1))
    rows = [[w, a * lam_fa['0.15'] / lam_avg] for w, a in rows]
    check = evaluate({'kind': 'eliashberg_custom', 'spectrum_meV': rows[:600], 'dos': {'phase': '1x1H', 'label': '0.14'},
                      'electrons_per_cell': 0.15, 'mu_star': MU})
    # 2. Ground-state boundary: energy of the 2x2 CDW end of the PES relative to 1x1 H (2x2 cells: tot_charge / 4).
    ry_mev = 13.605693122994e3
    dE = {}
    for label in CATALOG['pes_available']:
        e = evaluate({'kind': 'pes', 'label': label})['energies_Ry_per_cell']
        dE[int(label.split('_')[1]) / 100] = (e[-1] - e[0]) * ry_mev
    xs = sorted(dE)
    i = next(k for k in range(len(xs) - 1) if dE[xs[k]] > 0 >= dE[xs[k + 1]])
    xc = lerp(0, dE[xs[i]], xs[i], dE[xs[i + 1]], xs[i + 1]) / 4
    a_A = float(CATALOG['cells_and_doping']['qe_system_blocks']['1x1H (input/1x1H/sc.in)']['celldm(1)']) * 0.529177210903
    per_1e14 = 1 / (math.sqrt(3) / 2 * a_A ** 2 * 1e-16) / 1e14
    soft = evaluate({'kind': 'soft_mode', 'branch': '1x1H'})['rows']
    j = next(k for k, r in enumerate(soft) if r[1] < 0)
    x_soft = (soft[j - 1][0] + soft[j][0]) / 2 / per_1e14
    soft_cdw = evaluate({'kind': 'soft_mode', 'branch': '2x2CDW'})['rows']
    hyst = evaluate({'kind': 'mo_displacement', 'branch': '2x2CDW'})['rows']
    x_cdw_low = next(r[0] for r in hyst if r[1] > 0.005) / per_1e14
    # 3. SOC Tc of the 1x1 H phase at the boundary (Tc rises monotonically, so this is its ground-state maximum).
    soc = {x: evaluate({'kind': 'eliashberg', 'phase': '1x1H-SOC', 'label': x, 'flavor': 'FA', 'dos': 'full', 'mu_star': MU})
           for x in ('0.18', '0.19')}
    tc_soc = lerp(xc, 0.18, soc['0.18']['Tc_K'], 0.19, soc['0.19']['Tc_K'])
    # 4. High-doping CDW: identify dop1-3 by matching Allen-Dynes Tc of their spectra to the authors' table.
    table = evaluate({'kind': 'archived_tc_table', 'phase': '2x2CDW'})['rows']
    cdw, mapping = {}, {}
    for d in ('dop1', 'dop2', 'dop3'):
        r = evaluate({'kind': 'eliashberg', 'phase': '2x2CDW-high', 'label': d, 'flavor': 'FA', 'dos': 'constant', 'mu_star': MU})
        ad = allen_dynes(r['lambda'], r['omega_log_meV'], r['omega_2_meV'])
        row = min(table, key=lambda t: abs(t[3] - ad))
        cdw[d], mapping[d] = (r['Tc_K'], ad, row[3]), row[0] / 4
    by_x = {round(mapping[d], 3): cdw[d][0] for d in cdw}
    tc_cdw = {'0.35': by_x[0.35], '0.40': by_x[0.4]}
    submission = {'tc_1x1_K': tc_hidden, 'critical_doping_e_per_cell': xc, 'tc_1x1_soc_at_critical_doping_K': tc_soc,
                  'tc_cdw_high_K': tc_cdw, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    caption = {'dop1': 3.0, 'dop2': 4.3, 'dop3': 4.6}
    (APP/'report.md').write_text(f"""# Electron-doped MoS2: dome, CDW boundary and Tc (reference solution)

Doping x is electrons per 1x1 cell; 2x2 tot_charge is divided by 4 and n (10^14 cm^-2) by {per_1e14:.3f}
(a = {a_A:.4f} A). Tc is the isotropic Eliashberg Tc with mu* = 0.13 and the fully nonadiabatic (FA) spectra.

**Rise of Tc.** In the 1x1 H phase lambda_FA grows from {lam_fa['0.13']:.2f} (0.13) to {lam_fa['0.19']:.2f} (0.19) while
omega_log falls: the Q valleys fill and the M-point acoustic (A1) mode softens towards the Kohn anomaly. Adiabatic
lambda reaches {mom['lambda']['rows'][-1][1]:.2f} at 0.19, the precursor of the lattice instability; nonadiabatic shifts and
linewidths reduce it by half. With the energy-dependent DOS, Tc (full DOS) = {el['0.14']['Tc_K']:.2f}, {el['0.16']['Tc_K']:.2f}
and {el['0.18']['Tc_K']:.2f} K at 0.14, 0.16 and 0.18.

**Hidden dopings.** lambda jumps between 0.14 and 0.16 (0.77 -> 1.41), so interpolating Tc linearly is wrong at 0.15
(it gives {notes['0.15'][2]:.2f} K). I used the tabulated lambda and omega_log at 0.15 and 0.17 in Allen-Dynes
({notes['0.15'][0]:.2f}, {notes['0.17'][0]:.2f} K) times the Eliashberg/Allen-Dynes ratio of the neighbours
({notes['0.15'][1][0]:.3f}-{notes['0.17'][1][1]:.3f}): Tc = {tc_hidden['0.15']:.2f} K (0.15) and {tc_hidden['0.17']:.2f} K (0.17).
A cross-check with the neighbours' averaged spectrum rescaled to lambda(0.15) and the rigid-band 0.14 DOS gives
{check['Tc_K']:.2f} K. The full DOS matters because the low-temperature chemical potential sits inside the Q-valley edge
(mu0 = {el['0.14']['chemical_potential_noninteracting_meV']:.1f} meV above the 800 K Fermi level at 0.14).

**CDW boundary.** The PES (2x2 cell) gives E(CDW) - E(1x1) = {', '.join(f'{dE[x]:+.3f} meV at {x/4:.4f}' for x in xs)};
the 1x1 H phase stops being the ground state at x_c = {xc:.4f}. The 1x1 soft mode only becomes imaginary at
{x_soft:.4f} and the CDW survives (hysteresis) down to {x_cdw_low:.4f}: the coexistence window brackets x_c, and the
dynamical instability is not the thermodynamic boundary.

**SOC.** Tc with SOC (full DOS) is {soc['0.18']['Tc_K']:.2f} K (0.18) and {soc['0.19']['Tc_K']:.2f} K (0.19), about 25-30 %
below the scalar-relativistic values; at x_c the 1x1 H phase reaches {tc_soc:.2f} K.

**High-doping CDW.** The archive labels dop1-3 carry no doping and the paper's caption quotes
{caption['dop1']}, {caption['dop2']}, {caption['dop3']} x 10^14 cm^-2. Allen-Dynes Tc of the spectra
({', '.join(f'{d}: {cdw[d][1]:.2f}' for d in cdw)} K) matches the authors' table rows at
{', '.join(f'{4 * mapping[d]:.1f}' for d in cdw)} e per 2x2 cell, so n1 is about 4.0, not 3.0 x 10^14 cm^-2 (the archived
table is Allen-Dynes). Constant-DOS Eliashberg gives {tc_cdw['0.35']:.2f} K at 0.35 and {tc_cdw['0.40']:.2f} K at 0.40: Tc
falls as the M* mode hardens and the flat band lowers N(EF).

**Dome.** Tc rises in the 1x1 H phase up to x_c, the archived Allen-Dynes values of the 2x2 CDW in the coexistence
window are similar (30-34 K), the 1T' phase is lower (about 21 K) and the high-doping CDW falls to 6-13 K: a dome due to
the 1x1 -> 2x2 transition. The computed scale exceeds experiment (max about 11 K near 1.5 x 10^14 cm^-2) and sits at
higher doping. mu* was not tuned.

**Not supported.** The 800 K electronic temperature shifts the instability; the gating geometry, anharmonicity,
anisotropy, polaronic and larger-supercell CDWs (not computed in DFT here) and substrate effects are outside these data.
The SOC transfer assumes the PES boundary is unchanged by SOC. The archived A/NA columns of the 1x1 Tc table do not
match the archived A/NA spectra, so I did not use them.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
