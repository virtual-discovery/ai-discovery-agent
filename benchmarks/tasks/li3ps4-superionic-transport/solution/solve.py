"""Reference solution for li3ps4-superionic-transport; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
GRID = CATALOG['simulated_temperatures_K']
KB_EV, KB, E = 8.617333262e-5, 1.380649e-23, 1.602176634e-19
AMU = 1.66053906660e-27
N_LI = 288
EVIDENCE = []


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
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


CACHE = {}


def sublattice(f, T):
    if (f, T) not in CACHE:
        CACHE[(f, T)] = buy({'kind': 'sublattice_dynamics', 'functional': f, 'T_K': T})
    return CACHE[(f, T)]


def sigma(f, p, T, est):
    r = buy({'kind': 'conductivity', 'functional': f, 'initial_phase': p, 'T_K': T, 'estimator': est})
    return r['sigma_S_per_m'], r['block_error_S_per_m']


def angle(f, T):
    return sublattice(f, T)['ps4_polar_angle_rms_fluctuation_deg']['value']


def arrhenius(T, s, e):
    T, s, e = map(np.asarray, (T, s, e))
    c = np.polyfit(1 / (KB_EV * T), np.log(s * T), 1, w=s / e)
    return float(-c[0]), c


def extrapolate(c, T):
    return float(np.exp(np.polyval(c, 1 / (KB_EV * T))) / T)


def find_transition(f, threshold_deg=12.0):
    """Bisect the simulated grid between 600 and 900 K for the step in the PS4 angle fluctuation."""
    lo, hi = GRID.index(600), GRID.index(900)
    if angle(f, GRID[lo]) > threshold_deg or angle(f, GRID[hi]) < threshold_deg:
        raise RuntimeError('no angle step bracketed for ' + f)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if angle(f, GRID[mid]) > threshold_deg:
            hi = mid
        else:
            lo = mid
    return GRID[lo], GRID[hi]


def main():
    notes = {}
    # 1. Superionic transition of each potential: step of the PS4 angle fluctuation (discrete PS4 flips),
    #    confirmed by the heat capacity of that interval; the higher-T Cp peak coincides with P diffusion (melting).
    brackets = {f: find_transition(f) for f in ('PBEsol', 'PBE0', 'r2SCAN')}
    cp = {}
    for f, (lo, hi) in brackets.items():
        cp[f] = buy({'kind': 'heat_capacity', 'functional': f, 'T_low_K': lo, 'T_high_K': hi})['Cp_kB_per_atom']
    melt = {'PBE0': (850, 900), 'r2SCAN': (900, 950)}
    cp_melt = {f: buy({'kind': 'heat_capacity', 'functional': f, 'T_low_K': a, 'T_high_K': b})['Cp_kB_per_atom']
               for f, (a, b) in melt.items()}
    dp = {f: sublattice(f, 900).get('p_diffusion_coefficient_A2_per_ps', {}).get('value') for f in ('PBEsol', 'PBE0', 'r2SCAN')}
    tc = {f: (lo + hi) / 2 for f, (lo, hi) in brackets.items()}
    # 2. Superionic branch of the gamma-started runs (Green-Kubo, the full charge conductivity) and Haven ratio.
    plan = {'PBE0': {'gk': [775, 850, 1000], 'ne': [775, 850]},
            'r2SCAN': {'gk': [850, 900], 'ne': [850, 900]},
            'PBEsol': {'gk': [650, 700], 'ne': [650, 700]}}
    cell_T = {'PBEsol': 625, 'PBE0': 750, 'r2SCAN': 800}
    predictions, haven, sup = {}, {}, {}
    for f, p in plan.items():
        assert all(t > tc[f] for t in p['gk']), 'superionic points must lie above the transition'
        g = {t: sigma(f, 'gamma', t, 'green_kubo') for t in p['gk']}
        n = {t: sigma(f, 'gamma', t, 'nernst_einstein') for t in p['ne']}
        _, c = arrhenius(list(g), [v[0] for v in g.values()], [v[1] for v in g.values()])
        predictions[f'{f}/gamma/{cell_T[f]}'] = extrapolate(c, cell_T[f])
        haven[f] = float(np.mean([n[t][0] / g[t][0] for t in p['ne']]))
        sup[f] = g
    # 3. Activation energies of sigma (PBE0, gamma start): resolved Green-Kubo points below the transition
    #    (below ~600 K the Green-Kubo values sit at the noise floor), and the superionic branch to 1000 K.
    low = {t: sigma('PBE0', 'gamma', t, 'green_kubo') for t in (625, 675, 725)}
    assert all(v[0] > 3 * v[1] for v in low.values())
    ea_below, _ = arrhenius(list(low), [v[0] for v in low.values()], [v[1] for v in low.values()])
    g = sup['PBE0']
    ea_above, _ = arrhenius(list(g), [v[0] for v in g.values()], [v[1] for v in g.values()])
    # 4. Paddle-wheel test at 800 K (PBE0, superionic, below melting): Li hop rate from the tracer diffusivity vs
    #    PS4 reorientation rate from rotational diffusion.
    ne800, _ = sigma('PBE0', 'gamma', 800, 'nernst_einstein')
    volume = CATALOG['ideal_starting_cells_768_atoms']['gamma']['volume_A3']   # finite-T volume not archived
    d_li = ne800 * KB * 800 * volume * 1e-30 / (N_LI * E ** 2) * 1e8           # A^2/ps
    hop = 6 * d_li / 2.5 ** 2                                                   # per Li, hop length ~2.5 A
    domega = sublattice('PBE0', 800)['ps4_rotational_diffusion_rad2_per_ps']['value']
    flip = 4 * domega / math.radians(70.5) ** 2                                 # tau = phi^2 / (4 D_omega)
    log_ratio = math.log10(hop / flip)
    # 5. Li-S coupling: the barycentric Li-S term is the matrix recoil against the Li centre of mass; the
    #    dynamical part is read in the matrix frame.
    frames = {fr: buy({'kind': 'current_breakdown', 'frame': fr}) for fr in ('matrix', 'barycentric')}
    tau = np.array(frames['matrix']['tau_ps'])
    w = (tau >= 2) & (tau <= 8)
    frac = {fr: float(np.mean(np.array(d['Li_S_S_per_cm'])[w] / np.array(d['total_S_per_cm'])[w])) for fr, d in frames.items()}
    m = CATALOG['masses_amu']
    m_matrix = 96 * m['P'] + 384 * m['S']
    recoil = (2 * 384 * m['Li'] / m_matrix) / (1 + 288 * m['Li'] / m_matrix) ** 2
    submission = {'transition_temperature_K': {'PBE0': tc['PBE0'], 'r2SCAN': tc['r2SCAN']},
                  'sigma_predictions_S_per_m': predictions,
                  'haven_ratio_PBE0_gamma_750K': haven['PBE0'],
                  'activation_energy_eV_PBE0_gamma': {'below_transition': ea_below, 'superionic_to_1000K': ea_above},
                  'li_s_cross_fraction': frac['matrix'],
                  'log10_hop_to_flip_rate_ratio_PBE0_800K': log_ratio,
                  'paddle_wheel_supported': False,
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Li3PS4 transport from the released ML-MD data (reference solution)

**Transition.** For each potential the PS4 polar-angle fluctuation steps up between {brackets['PBEsol']} K (PBEsol),
{brackets['PBE0']} K (PBE0) and {brackets['r2SCAN']} K (r2SCAN): discrete PS4 flips set in. The heat capacity of
those intervals peaks ({cp}), and the gamma-started conductivity jumps there. A second, larger Cp peak
(PBE0 850-900 K: {cp_melt['PBE0']:.2f}; r2SCAN 900-950 K: {cp_melt['r2SCAN']:.2f} kB/atom) coincides with the onset of P
diffusion (D_P at 900 K: {dp}) and free PS4 rotation: melting of the PS4 sublattice, not the superionic transition.
Tc = {tc} K; it depends strongly on the DFT reference (PBEsol < PBE0 < r2SCAN). PBE0 is the most credible reference
(band gap, beta/alpha conductivities closest to experiment), but simulated Tc values carry finite-size and
time-scale limits; the gamma-started runs transform to a mixed beta/alpha arrangement, not pure alpha.

**Conductivity.** sigma is the Green-Kubo (charge-current) conductivity; Nernst-Einstein neglects ion-ion
correlations. In the superionic state of the gamma-started runs H_R = sigma_NE/sigma ~ {haven['PBE0']:.2f} (PBE0),
{haven['r2SCAN']:.2f} (r2SCAN), {haven['PBEsol']:.2f} (PBEsol): Li-Li correlated hops (first-shell Li-Li velocity
correlation) roughly double sigma relative to NE. The withheld cells lie just above each transition, so they were
predicted from Arrhenius fits of the superionic gamma-started branch (not by interpolating across the jump, and
not from alpha-started runs, which stay in a different, more conductive state with lower H_R): {predictions}.
Below ~600 K the Green-Kubo values sit at a noise floor (~0.03 S/m) and were excluded from fits.
E_a(sigma T), PBE0 gamma-start: {ea_below:.2f} eV below the transition, {ea_above:.3f} eV in the superionic phase
(the NE-based value is higher because H_R rises with T).

**Mechanism.** At 800 K (PBE0, superionic, below melting) D_Li ~ {d_li:.2f} A^2/ps gives ~{hop:.2f} hops/ps per Li,
whereas D_omega = {domega:.2e} rad^2/ps gives one PS4 flip per ~{1 / flip:.0f} ps: log10(ratio) = {log_ratio:.1f}.
PS4 librations are fast but are not reorientations. Li hops cannot be assisted one-by-one by PS4 flips: the
paddle-wheel picture is rejected; flips act by changing the Li site network (more available sites, lower barrier).
The barycentric Li-S term ({frac['barycentric']:.3f} of sigma) equals the kinematic recoil of the matrix against the Li
centre of mass ({recoil:.3f}); in the matrix frame it is {frac['matrix']:.3f}, i.e. negligible dynamical Li-S coupling.

**Not supported.** A unique experimental Tc, room-temperature extrapolations, or which pure polymorph forms above Tc.
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
