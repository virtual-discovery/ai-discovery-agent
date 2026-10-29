"""Generate the staged, lab, private and ablation data for li3ps4-superionic-transport from pool.json.

    cd benchmarks/construction/li3ps4-superionic-transport && python3 -B make_task.py

Every reference value is computed from author-archived arrays (released start); wrong paths and robustness
variations are computed from the same arrays with the analysis choices named in each key.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/li3ps4-superionic-transport'
P = json.loads((HERE/'pool.json').read_text())
GRID = P['temperature_grid_K']
F3 = ['PBEsol', 'r2SCAN', 'PBE0']
KB_EV = 8.617333262e-5
KB = 1.380649e-23
E = 1.602176634e-19
AMU = 1.66053906660e-27
N_LI, N_P, N_S = 288, 96, 384
M_LI, M_P, M_S = 6.941, 30.973762, 32.065       # masses used in the authors' LAMMPS input
BUDGET = 60
COSTS = {'nernst_einstein': 1, 'green_kubo': 2, 'sublattice_dynamics': 1, 'heat_capacity': 1, 'li_msd': 2,
         'current_breakdown': 1, 'li_pair_correlation': 1, 'order_parameters': 1}
# Withheld conductivity cells (both estimators): the first simulated temperature above each functional's
# transition in the runs started from gamma.
CELLS = {'PBEsol/gamma/625': ('PBEsol', 'gamma', 625), 'PBE0/gamma/750': ('PBE0', 'gamma', 750),
         'r2SCAN/gamma/800': ('r2SCAN', 'gamma', 800)}
# Beta-initialised series above 625 K are not served: in the archive the r2SCAN and PBE0 beta series coincide at
# 650-1150 K to 0.1 % (z ~ 0.01 against block errors of ~1 %), so at least one of them is mislabelled.
BETA_MAX_T = 625
THRESH = {'transition_K': 30.0, 'sigma_log10': 0.15, 'haven_ratio': 0.12, 'ea_below_eV': 0.35, 'ea_above_eV': 0.06,
          'log10_rate_ratio': 1.0, 'li_s_fraction': 0.10}
HOP_LENGTH_A = 2.5          # Li-site hop length used for the reference rate (alternatives in ROBUSTNESS)
FLIP_ANGLE_DEG = 70.5       # P-S bond reorientation in a +x -> -x PS4 flip (0 -> 70.5 deg and 109.5 -> 180 deg)
PS_BOND_A = 2.05


def series(f, p, est):
    s = P['conductivity'][f][p][est]
    T = np.array([float(t) for t in s])
    return T, np.array([v['sigma_S_per_m'] for v in s.values()]), np.array([v['block_error_S_per_m'] for v in s.values()])


def value(f, p, est, T):
    v = P['conductivity'][f][p][est][str(T)]
    return v['sigma_S_per_m'], v['block_error_S_per_m']


def arrhenius(T, s, e=None):
    """Fit sigma*T = A exp(-Ea/kT); returns (Ea, polyfit coefficients in x = 1/kT)."""
    x, y = 1 / (KB_EV * T), np.log(s * T)
    c = np.polyfit(x, y, 1, w=None if e is None else s / e)
    return float(-c[0]), c


def predict(c, T):
    return float(np.exp(np.polyval(c, 1 / (KB_EV * T))) / T)


def fig7(name, f):
    return np.array(P[name][f])


def melting_onset(f):
    """Lowest simulated T at which P atoms diffuse (D_P > 1e-3 A^2/ps)."""
    d = fig7('p_diffusion_A2_per_ps', f)
    return min(round(t / 25) * 25 for t, v, e in d if v > 1e-3)


def cp_intervals(f):
    cp = fig7('heat_capacity', f)
    out = []
    for (t, v, e) in cp:
        i = min(range(len(GRID) - 1), key=lambda k: abs((GRID[k] + GRID[k + 1]) / 2 - t))
        assert abs((GRID[i] + GRID[i + 1]) / 2 - t) < 0.5
        out.append((GRID[i], GRID[i + 1], v, e))
    return out


def transition(f):
    """Lower heat-capacity peak (below the melting onset) and the largest jump of the gamma tracer conductivity."""
    tm = melting_onset(f)
    lo, hi, v, e = max((c for c in cp_intervals(f) if c[1] < tm), key=lambda c: c[2])
    T, s, err = series(f, 'gamma', 'nernst_einstein')
    ok = (T < tm) & (s > 3 * err)
    jumps = [(np.log(s[i + 1] / s[i]), T[i], T[i + 1]) for i in range(len(T) - 1) if ok[i] and ok[i + 1]]
    j = max(jumps)
    assert (j[1], j[2]) == (lo, hi), (f, j, lo, hi)          # the two signatures agree for all functionals
    ang = fig7('angle_rms_deg', f)
    a = {round(t / 25) * 25: x for t, x, _ in ang}
    assert a[hi] - a[lo] == max(a[GRID[k + 1]] - a[GRID[k]] for k in range(len(GRID) - 1) if GRID[k + 1] < tm and GRID[k + 1] in a and GRID[k] in a)
    melt_peak = max((c for c in cp_intervals(f) if c[0] >= tm - 50), key=lambda c: c[2])
    return {'lo': lo, 'hi': hi, 'Tc': (lo + hi) / 2, 'melting_onset': tm, 'melting_peak': (melt_peak[0] + melt_peak[1]) / 2,
            'cp_peak': v, 'cp_peak_err': e}


def log_interp(f, Tw):
    i = GRID.index(Tw)
    t1, t2 = GRID[i - 1], GRID[i + 1]
    s1, s2 = value(f, 'gamma', 'green_kubo', t1)[0], value(f, 'gamma', 'green_kubo', t2)[0]
    w = (1 / t1 - 1 / Tw) / (1 / t1 - 1 / t2)
    return float(10 ** (math.log10(s1) + w * (math.log10(s2) - math.log10(s1))))


def superionic_points(f, Tc, Tw, upto=1000):
    return [t for t in GRID if Tc < t <= upto and t != Tw]


def li_diffusivity_from_ne(sigma_ne, T, volume_A3):
    """D_Li (A^2/ps) from the Nernst-Einstein conductivity: sigma_NE = N_Li e^2 D / (V k_B T)."""
    return sigma_ne * KB * T * volume_A3 * 1e-30 / (N_LI * E ** 2) * 1e8   # m^2/s -> A^2/ps is x 1e8


def flip_rate(domega, angle_deg=FLIP_ANGLE_DEG):
    """PS4 reorientation rate (1/ps) from rotational diffusion: tau = phi^2 / (4 D_omega)."""
    return 4 * domega / math.radians(angle_deg) ** 2


def at(name, f, T):
    for t, v, e in fig7(name, f):
        if abs(t - T) < 0.5:
            return v
    raise KeyError((name, f, T))


def main():
    tr = {f: transition(f) for f in F3}
    V0 = P['ideal_cells']['gamma']['volume_A3']
    # ---------------- references ----------------
    sigma_ref = {k: value(f, p, 'green_kubo', T)[0] for k, (f, p, T) in CELLS.items()}
    ne_hidden = {k: value(f, p, 'nernst_einstein', T)[0] for k, (f, p, T) in CELLS.items()}
    hr_ref = ne_hidden['PBE0/gamma/750'] / sigma_ref['PBE0/gamma/750']
    Tc0 = tr['PBE0']['Tc']
    T, gk, ge = series('PBE0', 'gamma', 'green_kubo')
    _, ne, nerr = series('PBE0', 'gamma', 'nernst_einstein')
    below = (T < Tc0)
    resolved = below & (gk > 3 * ge)
    start = max(t for t in T[below] if not (gk[T == t] > 3 * ge[T == t])[0]) if (~resolved[below]).any() else T[0]
    win_below = below & (T > start)                   # contiguous resolved run ending at the transition
    ea_below_ref, _ = arrhenius(T[win_below], gk[win_below], ge[win_below])
    win_above = (T > Tc0) & (T <= 1000)
    ea_above_ref, _ = arrhenius(T[win_above], gk[win_above], ge[win_above])
    d800 = li_diffusivity_from_ne(value('PBE0', 'gamma', 'nernst_einstein', 800)[0], 800, V0)
    hop800 = 6 * d800 / HOP_LENGTH_A ** 2
    dom800 = at('rotational_diffusion_rad2_per_ps', 'PBE0', 800)
    ratio_ref = math.log10(hop800 / flip_rate(dom800))
    cb = P['current_breakdown']
    tau = np.array(cb['tau_ps'])
    plateau = (tau >= 2) & (tau <= 8)

    def lis_fraction(frame):
        return float(np.mean(np.array(cb[frame]['Li_S_S_per_cm'])[plateau] / np.array(cb[frame]['total_S_per_cm'])[plateau]))
    lis_ref = lis_fraction('matrix')
    # Pure recoil: in the barycentric frame the matrix moves rigidly against the Li centre of mass, so the released
    # Li-S term (one cross term) over the total is (2 N_S m_Li / M_matrix) / (1 + N_Li m_Li / M_matrix)^2.
    m_matrix = N_P * M_P + N_S * M_S
    recoil = (2 * N_S * M_LI / m_matrix) / (1 + N_LI * M_LI / m_matrix) ** 2
    hidden = {
        'targets': {'transition_temperature_K': {'PBE0': tr['PBE0']['Tc'], 'r2SCAN': tr['r2SCAN']['Tc']},
                    'sigma_S_per_m': sigma_ref, 'haven_ratio_PBE0_gamma_750K': hr_ref,
                    'activation_energy_eV': {'below_transition': ea_below_ref, 'superionic_to_1000K': ea_above_ref},
                    'li_s_cross_fraction': lis_ref, 'log10_hop_to_flip_rate_ratio_PBE0_800K': ratio_ref,
                    'paddle_wheel_supported': False},
        'withheld_tracer_values_S_per_m': ne_hidden,
        'thresholds': THRESH,
        'reference_details': {'transitions': tr, 'ea_below_window_K': [float(T[win_below].min()), float(T[win_below].max())],
                              'ea_above_window_K': [float(T[win_above].min()), float(T[win_above].max())],
                              'D_Li_800K_A2_per_ps': d800, 'hop_rate_800K_per_ps': hop800, 'D_omega_800K': dom800,
                              'flip_rate_800K_per_ps': flip_rate(dom800), 'ideal_gamma_volume_A3': V0,
                              'barycentric_recoil_prediction': recoil}}
    # ---------------- ablations ----------------
    A, R = {}, {}
    A['sigma_tracer'] = dict(ne_hidden)
    A['haven_one'] = 1.0
    ea_above_ne, _ = arrhenius(T[win_above], ne[win_above], nerr[win_above])
    A['ea_tracer'] = {'below_transition': arrhenius(T[win_below & (ne > 3 * nerr)], ne[win_below & (ne > 3 * nerr)],
                                                    nerr[win_below & (ne > 3 * nerr)])[0], 'superionic_to_1000K': ea_above_ne}
    noisy = below & (gk > 0)
    A['ea_noise_floor'] = {'below_transition': arrhenius(T[noisy], gk[noisy])[0], 'superionic_to_1000K': ea_above_ref}
    A['tc_melting_peak'] = {f: tr[f]['melting_peak'] for f in ('PBE0', 'r2SCAN')}
    low_branch = {}
    for k, (f, p, Tw) in CELLS.items():
        Tf, g, e = series(f, 'gamma', 'green_kubo')
        m = (Tf < Tw) & (g > 3 * e)
        idx = np.where(m)[0][-4:]
        low_branch[k] = predict(arrhenius(Tf[idx], g[idx], e[idx])[1], Tw)
    A['sigma_low_branch'] = low_branch
    A['tc_pbesol_for_all'] = {f: tr['PBEsol']['Tc'] for f in ('PBE0', 'r2SCAN')}
    A['tc_experimental_alpha'] = {'PBE0': 746.0, 'r2SCAN': 746.0}
    A['sigma_interpolated_across_jump'] = {k: log_interp(f, Tw) for k, (f, p, Tw) in CELLS.items()}
    A['sigma_alpha_branch'] = {k: value(f, 'alpha', 'green_kubo', Tw)[0] for k, (f, p, Tw) in CELLS.items()}
    A['haven_alpha_branch'] = value('PBE0', 'alpha', 'nernst_einstein', 750)[0] / value('PBE0', 'alpha', 'green_kubo', 750)[0]
    A['li_s_barycentric'] = lis_fraction('barycentric')
    # libration taken for reorientation: harmonic libration frequency from the angle fluctuation at 800 K
    inertia = 8 / 3 * M_S * AMU * (PS_BOND_A * 1e-10) ** 2
    theta = math.radians(at('angle_rms_deg', 'PBE0', 800))
    nu_lib = math.sqrt(KB * 800 / inertia) / theta / (2 * math.pi) * 1e-12        # 1/ps
    A['log_ratio_libration'] = math.log10(hop800 / nu_lib)
    allres = (T >= T[win_below].min()) & (T <= 1000) & (T != 750) & (gk > 3 * ge)
    A['ea_single_arrhenius'] = arrhenius(T[allres], gk[allres], ge[allres])[0]
    # ---------------- robustness ----------------
    ne_branch, sup_fit = {}, {}
    hr3, hr_all = None, None
    for k, (f, p, Tw) in CELLS.items():
        pts = superionic_points(f, tr[f]['Tc'], Tw)
        Tf, g, e = series(f, 'gamma', 'green_kubo')
        _, n, ne_ = series(f, 'gamma', 'nernst_einstein')
        m = np.isin(Tf, pts)
        sup_fit[k] = predict(arrhenius(Tf[m], g[m], e[m])[1], Tw)
        h3 = float(np.mean(n[m][:3] / g[m][:3]))
        ne_branch[k] = predict(arrhenius(Tf[m], n[m], ne_[m])[1], Tw) / h3
        if f == 'PBE0':
            hr3 = h3
            m2 = np.isin(Tf, superionic_points(f, tr[f]['Tc'], Tw, upto=1100))
            hr_all = float(np.mean(n[m2] / g[m2]))
    R['sigma_tracer_branch_over_haven'] = ne_branch
    R['haven_first_three'] = hr3
    R['sigma_superionic_gk_fit'] = sup_fit
    R['haven_all_superionic_to_1100K'] = hr_all
    m875 = (T > Tc0) & (T <= 875)
    R['ea_above_to_875K'] = arrhenius(T[m875], gk[m875], ge[m875])[0]
    R['ea_below_unweighted'] = arrhenius(T[win_below], gk[win_below])[0]
    R['ea_below_tracer'] = A['ea_tracer']['below_transition']
    m_ex = win_above & (T != 750)
    R['ea_above_without_750K'] = arrhenius(T[m_ex], gk[m_ex], ge[m_ex])[0]
    R['tc_upper_bracket'] = {f: float(tr[f]['hi']) for f in ('PBE0', 'r2SCAN')}
    R['tc_lower_bracket'] = {f: float(tr[f]['lo']) for f in ('PBE0', 'r2SCAN')}
    R['log_ratio_domega_as_rate'] = math.log10(hop800 / dom800)
    msd = P['li_msd_750K']
    t, y = np.array(msd['t_ps']), np.array(msd['msd_A2'])
    w = (t >= 20) & (t <= 300)
    slope = np.polyfit(t[w], y[w], 1)[0]
    d750 = slope / 6
    R['D_Li_750K_msd'] = d750
    sigma_ne_750 = d750 * 1e-8 * N_LI * E ** 2 / (V0 * 1e-30 * KB * 750)
    R['sigma_PBE0_750_msd_route'] = sigma_ne_750 / hr3
    R['log_ratio_msd_750_hop3A'] = math.log10((6 * d750 / 3.0 ** 2) / flip_rate(at('rotational_diffusion_rad2_per_ps', 'PBE0', 750)))
    # ---------------- staged files ----------------
    conductivity_available = {}
    lab_cond = {}
    for f in F3:
        for p in ('alpha', 'beta', 'gamma'):
            temps = [t for t in GRID if (p != 'beta' or t <= BETA_MAX_T) and f'{f}/{p}/{t}' not in CELLS]
            conductivity_available[f'{f}/{p}'] = temps
            for est in ('nernst_einstein', 'green_kubo'):
                for t_ in temps:
                    lab_cond[f'{f}/{p}/{t_}/{est}'] = P['conductivity'][f][p][est][str(t_)]
    sub_available, lab_sub = {}, {}
    for f in F3:
        rows = {}
        for name, key in (('angle_rms_deg', 'ps4_polar_angle_rms_fluctuation_deg'),
                          ('rotational_diffusion_rad2_per_ps', 'ps4_rotational_diffusion_rad2_per_ps'),
                          ('p_diffusion_A2_per_ps', 'p_diffusion_coefficient_A2_per_ps')):
            for t_, v, e in fig7(name, f):
                rows.setdefault(round(t_ / 25) * 25, {})[key] = {'value': v, 'block_error': e}
        sub_available[f] = sorted(rows)
        for t_, row in rows.items():
            lab_sub[f'{f}/{t_}'] = row
    cp_available, lab_cp = {}, {}
    for f in F3:
        cp_available[f] = []
        for lo, hi, v, e in cp_intervals(f):
            cp_available[f].append([lo, hi])
            lab_cp[f'{f}/{lo}/{hi}'] = {'Cp_kB_per_atom': v, 'block_error': e}
    catalog = {
        'material': 'Li3PS4 (lithium ortho-thiophosphate); NpT ML-MD of a 768-atom cell (288 Li, 96 P, 384 S) at p = 0, '
                    'released by Gigli, Tisi, Grasselli and Ceriotti (2024)',
        'potentials': {'PBEsol': 'ML potential trained on PBEsol DFT', 'r2SCAN': 'ML potential trained on r2SCAN DFT',
                       'PBE0': 'ML potential trained on PBE0 (hybrid) DFT'},
        'initial_phases': {'alpha': 'started from ideal alpha-Li3PS4 (Cmcm)', 'beta': 'started from ideal beta-Li3PS4 (Pnma)',
                           'gamma': 'started from ideal gamma-Li3PS4 (Pmn2_1), the stable room-temperature polymorph'},
        'simulated_temperatures_K': GRID,
        'run_lengths': 'about 6 ns (gamma), 4 ns (beta) and 2 ns (alpha) per temperature, 2 fs time step',
        'ideal_starting_cells_768_atoms': P['ideal_cells'],
        'ionic_charges': {'Li': 1, 'P': 5, 'S': -2},
        'masses_amu': {'Li': M_LI, 'P': M_P, 'S': M_S},
        'purchase_costs': {'conductivity (nernst_einstein)': 1, 'conductivity (green_kubo)': 2, **{k: v for k, v in COSTS.items()
                                                                                                 if k not in ('nernst_einstein', 'green_kubo')}},
        'budget': BUDGET,
        'requests': {
            'conductivity': {'request': {'kind': 'conductivity', 'functional': '<PBEsol|r2SCAN|PBE0>', 'initial_phase': '<alpha|beta|gamma>',
                                         'T_K': '<simulated temperature>', 'estimator': '<nernst_einstein|green_kubo>'},
                             'returns': 'conductivity in S/m with its block-averaging error and the run mean temperature. '
                                        'nernst_einstein: from the Li self-diffusion coefficient (tracer MSD). '
                                        'green_kubo: from the total charge-current autocorrelation (Helfand-Einstein form).',
                             'available': conductivity_available},
            'sublattice_dynamics': {'request': {'kind': 'sublattice_dynamics', 'functional': '<functional>', 'T_K': '<temperature>'},
                                    'returns': 'whichever of these the archive holds at that temperature: rms fluctuation of the '
                                               'P-S polar angle about its local equilibrium (deg), PS4 rotational diffusion '
                                               'coefficient from the P-S bond angular-velocity autocorrelation (rad^2/ps), and P '
                                               'linear diffusion coefficient (A^2/ps); with block errors',
                                    'source': "the authors' NpT temperature series used for their phase-transition analysis "
                                              '(initial phase not labelled in the archive)',
                                    'available': sub_available},
            'heat_capacity': {'request': {'kind': 'heat_capacity', 'functional': '<functional>', 'T_low_K': '<T>', 'T_high_K': '<next T>'},
                              'returns': 'isobaric heat capacity in k_B per atom from the finite difference of the mean enthalpy '
                                         'between two consecutive simulated temperatures, with block error (same series as '
                                         'sublattice_dynamics)',
                              'available': cp_available},
            'li_msd': {'request': {'kind': 'li_msd'},
                       'returns': 'Li mean-squared displacement (A^2) with block errors versus time (0-654 ps; every 10 fs to 2 ps, '
                                  'every 0.1 ps to 20 ps, then every 1 ps) for the ML-PBE0 run started from gamma at 750 K'},
            'current_breakdown': {'request': {'kind': 'current_breakdown', 'frame': '<matrix|barycentric>'},
                                  'returns': 'running Green-Kubo integral sigma(tau) in S/cm, 0-8 ps every 20 fs, of the total '
                                             'conductivity and of the Li-S cross-correlation term, for one superionic ML-MD run '
                                             '(the archive does not state its functional or temperature). matrix: velocities '
                                             'relative to the centre of mass of the P and S atoms; barycentric: centre of mass '
                                             'of all atoms fixed'},
            'li_pair_correlation': {'request': {'kind': 'li_pair_correlation'},
                                    'returns': 'distance-resolved integral of the Li-Li, Li-P and Li-S velocity cross-correlation '
                                               '(A^2/ps^2, matrix frame) and the Li-A radial distribution functions, from 190 ps of '
                                               'an alpha-phase run at 650 K'},
            'order_parameters': {'request': {'kind': 'order_parameters', 'T_K': '<700|750|800>'},
                                 'returns': 'samples of the PS4 orientational order parameters (s1, s2) along the ML-PBE0 run '
                                            'started from gamma at that temperature, plus the centroids of reference samples '
                                            'from alpha, beta and gamma runs below the transition. Ideal structures: '
                                            'alpha (-1, 0), beta (0, 0), gamma (1, 1)'}},
        'to_predict': list(CELLS),
        'notes': ['Results are author-archived analyses of the released simulations; no new simulation is run.',
                  'Beta-initialised conductivities are served only up to 625 K: above that the archived r2SCAN and PBE0 beta '
                  'series coincide to 0.1 %, far below their block errors, so at least one of them is mislabelled.',
                  'The conductivity cells listed in to_predict are not served with either estimator.']}
    lab = {'conductivity': lab_cond, 'sublattice_dynamics': lab_sub, 'heat_capacity': lab_cp, 'li_msd': P['li_msd_750K'],
           'current_breakdown': P['current_breakdown'], 'li_pair_correlation': P['li_pair_correlation'],
           'order_parameters': P['order_parameters'], 'costs': COSTS}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    tg = hidden['targets']
    values = [tg['transition_temperature_K']['PBE0'], tg['transition_temperature_K']['r2SCAN'], *sigma_ref.values(), hr_ref,
              ea_below_ref, ea_above_ref, lis_ref, ratio_ref, *ne_hidden.values()]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': R}, indent=1) + '\n')
    for name, path in (('public', TASK/'environment/lab/public.json'), ('tests', TASK/'tests/public.json')):
        d = json.loads(path.read_text())
        d['solver_budget'] = BUDGET
        path.write_text(json.dumps(d, indent=1) + '\n')
    dex = lambda pred, k: math.log10(pred / sigma_ref[k])
    report = {'transitions': {f: {k: tr[f][k] for k in ('lo', 'hi', 'Tc', 'melting_onset', 'melting_peak')} for f in F3},
              'sigma_ref': sigma_ref, 'hr_ref': hr_ref, 'ea_below_ref': ea_below_ref, 'ea_below_window': hidden['reference_details']['ea_below_window_K'],
              'ea_above_ref': ea_above_ref, 'ratio_ref': ratio_ref, 'lis_ref': lis_ref, 'recoil_prediction': recoil,
              'abl': {k: ({kk: (round(vv, 4), round(dex(vv, kk), 3)) for kk, vv in v.items()} if isinstance(v, dict) and set(v) == set(CELLS) else v)
                      for k, v in A.items()},
              'rob': {k: ({kk: (round(vv, 4), round(dex(vv, kk), 3)) for kk, vv in v.items()} if isinstance(v, dict) and set(v) == set(CELLS) else v)
                      for k, v in R.items()},
              'catalog_kb': round((TASK/'environment/agent/catalog.json').stat().st_size / 1024, 1),
              'lab_kb': round((TASK/'environment/lab/data.json').stat().st_size / 1024, 1)}
    print(json.dumps(report, indent=1, default=float))


if __name__ == '__main__':
    main()
