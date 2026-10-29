"""Reference solution for bab2-hole-doping-two-gap; runs as the Harbor oracle and in check_task."""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())


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
        time.sleep(1)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    return r['id'], r['result']


def tot_charge(result):
    return float(re.search(r'tot_charge\s*=\s*([\d.]+)', result['scf.in']).group(1))


def soft_mode(q_points):
    """Lowest frequency over the q-points, excluding the three acoustic modes at Gamma; returns (freq, q)."""
    best = (np.inf, None)
    for d in q_points:
        f = sorted(d['frequencies_cm1'])
        if np.allclose(d['q_cartesian_2pi_over_alat'], 0):
            f = f[3:]
        if f[0] < best[0]:
            best = (f[0], d['q_cartesian_2pi_over_alat'])
    return best


def clusters(rows, T):
    a = np.array(rows)
    wgt, gap = a[:, 0] - T, a[:, 1]          # EPW gap0 file: column 1 = T + normalised weight, column 2 = gap (meV)
    keep = wgt > 0
    gap, wgt = gap[keep], wgt[keep]
    o = np.argsort(gap)
    gap, wgt = gap[o], wgt[o]
    k = int(np.argmax(np.diff(gap)))
    return [float(np.sum(g * w) / np.sum(w)) for g, w in ((gap[:k + 1], wgt[:k + 1]), (gap[k + 1:], wgt[k + 1:]))]


def main():
    evidence = []
    # 1. Map archive folders to dopings from each run's own scf.in (the README lists dopings in another order).
    mapping, logs = {}, {}
    for entry in CATALOG['runs']:
        if not entry['run'].startswith('sensitivity/'):
            continue
        jid, res = evaluate({'kind': 'iso_log', 'run': entry['run']})
        evidence.append(jid)
        x = tot_charge(res)
        mapping[x] = entry['run']
        logs[x] = res
    run_of = lambda x: mapping[min(mapping, key=lambda v: abs(v - x))]
    nan_runs = [x for x, r in logs.items() if 'NaN' in r['iso.out (isotropic Eliashberg section)']]
    # 2. Soft mode versus doping from the purchasable DFPT runs bracketing the instability.
    withheld = {w.split()[-1] for w in CATALOG['withheld']}
    soft = {}
    for x in (0.04, 0.06, 0.10):
        assert run_of(x) not in withheld
        jid, res = evaluate({'kind': 'dfpt_phonons', 'run': run_of(x)})
        evidence.append(jid)
        assert abs(tot_charge(res) - x) < 1e-9
        soft[x] = soft_mode(res['q_points'])
    xs = sorted(soft)
    w2 = np.array([np.sign(soft[x][0]) * soft[x][0] ** 2 for x in xs])   # omega^2 is analytic across the instability
    p = np.polyfit(xs, w2, 2)
    roots = [r.real for r in np.roots(p) if abs(r.imag) < 1e-12 and xs[0] < r.real < xs[-1]]
    x_c = float(min(roots))
    w_hom = float(np.sqrt(np.polyval(p, 1 / 9)))
    # 3. Isotropic Eliashberg Tc (mu* = 0.10, 500 meV cutoff) at x = 0.065 and 0.10; validate against the archive first.
    jid, val = evaluate({'kind': 'eliashberg_iso', 'run': run_of(0.065), 'mu_star': 0.13, 'cutoff_meV': 500,
                         'temperatures_K': [15, 21, 27]})
    evidence.append(jid)
    archived = {}
    T = None
    for line in logs[0.065]['iso.out (isotropic Eliashberg section)'].splitlines():
        m = re.match(r'\s*temp\(\s*\d+\)\s*=\s*([\d.]+) K', line)
        if m:
            T = float(m.group(1))
        parts = line.split()
        if len(parts) == 4 and parts[0].isdigit() and parts[3] != 'NaN':
            archived[T] = float(parts[3])
    validation = {r['temperature_K']: (r['gap_meV'], archived[r['temperature_K']]) for r in val['results']}
    iso_tc = {}
    for key, x, coarse in (('0.065', 0.065, [28, 31, 33, 34, 35, 36, 37, 39]), ('0.10', 0.10, [26, 29, 31, 32, 33, 34, 35, 37])):
        settings = {'kind': 'eliashberg_iso', 'run': run_of(x), 'mu_star': 0.10, 'cutoff_meV': 500,
                    'max_iter': 5000, 'conv_thr': 1e-6, 'mixing': 0.8}
        jid, res = evaluate({**settings, 'temperatures_K': coarse})
        evidence.append(jid)
        good = [r for r in res['results'] if r['converged'] and r['gap_meV'] > 0]
        t_hi = max(r['temperature_K'] for r in good)
        fine = [round(t_hi - 1.0 + 0.25 * i, 2) for i in range(8)]
        jid, res2 = evaluate({**settings, 'temperatures_K': fine})
        evidence.append(jid)
        pts = sorted({r['temperature_K']: r['gap_meV'] for r in res['results'] + res2['results']
                      if r['converged'] and r['gap_meV'] > 0}.items())
        (t1, d1), (t2, d2) = pts[-2], pts[-1]
        iso_tc[key] = float(t2 + (t2 - t1) * d2 ** 2 / (d1 ** 2 - d2 ** 2))
    # 4. Anisotropic run at x = 0.065: gap structure at 25 K and closing temperature beyond the archived 60 K.
    aniso = {}
    for T in (25, 50, 55, 60):
        jid, res = evaluate({'kind': 'aniso_gap', 'run': run_of(0.065), 'temperature_K': T, 'axis': 'imaginary'})
        evidence.append(jid)
        line = next(l for l in res['epw_int.out (lines for this temperature)'] if 'Min. / Max.' in l)
        aniso[T] = {'dmax': float(line.split('=')[1].split()[1]), 'clusters': clusters(res['rows'], T)}
    small, large = aniso[25]['clusters']
    ratio = large / small
    d55, d60 = aniso[55]['dmax'], aniso[60]['dmax']
    tc_aniso = float(60 + 5 * d60 ** 2 / (d55 ** 2 - d60 ** 2))
    # 5. Supercells: check what each archived supercell actually contains before using it.
    cells = {}
    for name in ('331', '331_cs'):
        jid, res = evaluate({'kind': 'supercell_phonons', 'run': name})
        evidence.append(jid)
        nel = float(re.search(r'number of electrons\s*=\s*([\d.]+)', ' '.join(res['scf.out (selected lines)'])).group(1))
        folded = [d for d in res['q_points'] if not np.allclose(d['q_cartesian_2pi_over_alat'], 0)][0]
        cells[name] = {'electrons': nel, 'lowest_at_folded_H': min(folded['frequencies_cm1'])}
    assert cells['331']['electrons'] == 144.0      # 9 x 16 valence electrons: the '331' cell is undoped
    w_cs = cells['331_cs']['lowest_at_folded_H']

    submission = {'critical_hole_doping': x_c, 'iso_tc_K': iso_tc, 'aniso_tc_K': tc_aniso, 'gap_ratio_25K': ratio,
                  'soft_mode_cm1': {'homogeneous_one_ninth': w_hom, 'cs_substituted': w_cs}, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    mapping_text = ', '.join(f"{run} -> x = {x:g}" for x, run in sorted(mapping.items(), key=lambda kv: kv[1]))
    (APP/'report.md').write_text(f"""# Hole-doped BaB2 (reference solution)

**Run mapping.** The archive README lists the sensitivity dopings as 0.0, 0.02, 0.04, 0.06, 0.065, 0.08, 0.1, but the
folders are not in that order. From each run's own scf.in (tot_charge) and Fermi energy: {mapping_text}.

**Dynamical stability.** The lowest DFPT mode sits at H = (1/3, 1/3, 1/2) (Ba out-of-plane motion coupled to B in-plane
displacements): {', '.join(f'x = {x:g}: {soft[x][0]:.2f} cm-1' for x in xs)}. omega^2, not omega, is analytic in the doping
near a displacive instability, and omega^2(x) is nearly linear; a quadratic fit of omega^2 gives the stabilisation
threshold x_c = {x_c:.4f} holes per formula unit (omega itself would put it near 0.070). For x < x_c the harmonic alpha2F
is not physical: the coupling at x = 0 (lambda = 2.10) is inflated by the soft branch near H, while N(E_F) changes by only
about 5% across the series, so the low-doping Tc values describe a structure that is not stable. The archived DFPT of the
x = 0.065 and 0.08 runs were not available here.

**Isotropic Eliashberg Tc** (mu* = 0.10, 500 meV Matsubara cutoff, archived alpha2F, new numpy solutions). The solver
reproduces the archived EPW gaps of the x = 0.065 run at mu* = 0.13 ({', '.join(f'{T:g} K: {a:.3f} vs {b:.3f} meV' for T, (a, b) in validation.items())}).
Close to Tc the fixed-point iteration slows down critically; runs that do not converge leave a spurious finite gap above
Tc, so only converged gaps were used and Delta^2 was extrapolated to zero: Tc = {iso_tc['0.065']:.2f} K at x = 0.065 and
{iso_tc['0.10']:.2f} K at x = 0.10. The archived logs (mu* = 0.13, 15-30 K) do not bracket Tc, and the x = {', '.join(f'{x:g}' for x in nan_runs)} log
ends in NaN with "Convergence was not reached" at 30 K, where the gap collapses; the header's Allen-Dynes estimate
underestimates the strong-coupling Eliashberg Tc.

**Anisotropic solution and gaps.** The archived anisotropic run (mu* = 0.13) at x = 0.065 has two well separated gap
clusters on the imaginary axis: at 25 K the weight-averaged gaps are {large:.2f} meV (sigma-like cylinders around Gamma-A)
and {small:.2f} meV (pi-like pockets), a ratio of {ratio:.2f}. The gap is still open at the last archived temperature
(60 K: max gap {d60:.2f} meV); extrapolating Delta_max^2 linearly from 55 and 60 K gives Tc(aniso) = {tc_aniso:.1f} K, about
twice the isotropic value at the same mu*, because the strongly coupled sigma states pair largely by themselves. The
real-axis (Pade) gap distributions become erratic (negative gaps) above about 30 K, so they were not used.

**Realisability.** The archived '331' supercell, described in the README as homogeneously 1-hole doped, contains 144
valence electrons (9 x 16), the undoped Fermi energy (11.3085 eV) and the undoped cell; its folded-H mode
({cells['331']['lowest_at_folded_H']:.2f} cm-1) reproduces the undoped primitive cell. Homogeneous doping of 1/9 h/f.u. was therefore taken from the
primitive-cell series: omega_H = {w_hom:.1f} cm-1 (stable). The Cs-substituted cell (143 electrons, 1/9 h/f.u.) stays unstable
at the folded H point (lowest {w_cs:.2f} cm-1): a local Cs dopant does not act like a homogeneous hole background.

**Not supported by these data:** convergence of alpha2F and of the anisotropic gaps with k/q grids; mu* itself (Tc moves
by about 8% between mu* = 0.10 and 0.13); the anisotropic Tc beyond extrapolation from 60 K; the exact x_c below the
spacing of the purchasable dopings; any Tc for x < x_c; and whether any chemical route reproduces homogeneous doping.
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
