"""Generate the staged, lab, private and ablation data for cnt-kohn-anomaly-nonadiabatic.

    bash run_in_image.sh make_task.py          (runtime image; reads /cache/pool from extract.py)

Steps:
  1. checks that the lab evaluator (selfenergy.py) reproduces the archived self-energy tables;
  2. computes the reference answers from the archived tables (all 37 temperatures, including the
     withheld ones) and the archived threshold model;
  3. computes the robustness variants and the ablation (wrong-path) values from the same data;
  4. writes the catalog, lab data, private targets, ablation values, public budgets and calibration.json.
"""
import json
import math
import shutil
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/cnt-kohn-anomaly-nonadiabatic'
LAB = TASK/'environment/lab'
RELEASED = [50.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1578.87]
TARGET_T = [100.0, 300.0]
TSTAR_T = [500.0, 510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0]
M_STAR = 2.0
BUDGET = 26
THRESHOLDS = {'omega_nonadiabatic_abs_meV': 1.2, 'linewidth_rel': 0.06, 'omega_adiabatic_abs_meV': 1.5,
              'instability_temperature_rel': 0.08, 'threshold_coupling_rel': 0.06}
EQ9 = (7.278, 5.186, 0.763, 0.518)          # alpha, beta, delta, zeta of the source's Eq. (9)
TB = M.Tables()


def table_name(t):
    return '1578.87' if abs(t - 1578.87) < 1e-6 else str(int(t))


def lab_pi(nu, t, omegas, nk=1000000, eta=None, processes=('k+q',), pairs=((1, 1), (1, 2), (2, 1), (2, 2))):
    eta = M.KB_EV * t / 100 * 1000 if eta is None else eta
    return M.SE.self_energy(M.couplings(nu, processes), t, omegas, nk, eta, pairs)


def check_evaluator():
    """Max deviation of the lab evaluator (nk = 10^6, eta = kT/100) from archived tables at several T."""
    dev = {}
    for nu in (1, 20, 23):
        for t in (50.0, 100.0, 300.0, 520.0, 1578.87):
            wq, w, re, im, st = TB.get(nu, t)
            idx = [int(np.argmin(abs(w - x))) for x in (0.0, wq * 0.9, wq)]
            r, i = lab_pi(nu, t, w[idx])
            dev[f'{nu}@{t:g}'] = {'re': float(np.max(np.abs(r - re[idx]))), 'im': float(np.max(np.abs(i - im[idx]))),
                                  'static': float(abs(r[0] - st))}
    worst = max(max(v.values()) for v in dev.values())
    assert worst < 0.2, dev
    return dev, worst


def solve_set(pdc_of, wq_of=None, scale=1.0, gamma_term=True, source='archive'):
    """Omega_NA / gamma at 100, 300 K and Omega_ad(100 K) for nu = 20 with a given subtraction rule."""
    out = {'omega_nonadiabatic_meV': {}, 'linewidth_meV': {}}
    for t in TARGET_T:
        wq, w, re, im, st = TB.get(20, t)
        wq = wq if wq_of is None else wq_of(20)
        o, g = M.pole(w, re, im, wq, pdc_of(20, t), gamma_term=gamma_term, scale=scale)
        out['omega_nonadiabatic_meV'][table_name(t)] = o
        out['linewidth_meV'][table_name(t)] = g
    wq, w, re, im, st = TB.get(20, 100.0)
    wq = wq if wq_of is None else wq_of(20)
    o2 = wq ** 2 + 2 * wq * (st - pdc_of(20, 100.0)) * scale
    out['omega_adiabatic_100K_meV'] = math.copysign(math.sqrt(abs(o2)), o2)
    return out


def tstar(pdc_of, wq_of=None, scale=1.0, gamma_term=True, temps=None):
    temps = TB.temps if temps is None else temps
    fs = []
    for t in temps:
        wq, w, re, im, st = TB.get(1, t)
        wq = wq if wq_of is None else wq_of(1)
        fs.append(M.instability_function(wq, st, M.im_slope(w, im), pdc_of(1, t), scale, gamma_term))
    if all(f < 0 for f in fs):
        return max(temps)                     # unstable at every archived temperature
    return M.crossing(temps, fs)


def main():
    dev, worst = check_evaluator()
    pdc = {nu: TB.get(nu, 1578.87)[4] for nu in (1, 20, 23)}
    assert abs(M.T_SIGMA - 1578.87) < 0.01
    std = lambda nu, t: pdc[nu]
    ref = solve_set(std)
    ref['rbm_instability_temperature_K'] = tstar(std)
    wq = {nu: TB.get(nu, 1578.87)[0] for nu in (1, 20, 23)}
    g0 = {nu: M.threshold_model(1.0 / M_STAR, wq[nu]) for nu in (1, 20)}
    ref['threshold_coupling_meV'] = {str(nu): g0[nu] for nu in (1, 20)}
    # archived Fig. 1 output reproduced by the rebuilt model (a = 0.5, omega = 100 meV -> 123.655 meV)
    assert abs(M.threshold_model(0.5, 100.0) - 123.65510) < 1e-3

    # ---------------- robustness variants (defensible; must pass) ----------------
    rob = {}
    auth = {'omega_nonadiabatic_meV': {}, 'linewidth_meV': {}}
    for t in TARGET_T:
        w_, w, re, im, st = TB.get(20, t)
        o, g = M.authors_iteration(w, w_, re, im, pdc[20])
        auth['omega_nonadiabatic_meV'][table_name(t)], auth['linewidth_meV'][table_name(t)] = o, g
    signed = []
    for t in TSTAR_T:
        w_, w, re, im, st = TB.get(1, t)
        o, _ = M.authors_iteration(w, w_, re, im, pdc[1])
        signed.append(math.copysign(o * o, o))
    auth['rbm_instability_temperature_K'] = M.crossing(TSTAR_T, signed)
    rob['authors_iteration'] = auth
    nog = solve_set(std, gamma_term=False)
    rob['no_linewidth_term'] = {'omega_nonadiabatic_meV': nog['omega_nonadiabatic_meV'],
                                'rbm_instability_temperature_K': tstar(std, gamma_term=False)}
    onshell = {}
    for t in TARGET_T:
        w_, w, re, im, st = TB.get(20, t)
        onshell[table_name(t)] = float(-np.interp(w_, w, im))
    rob['onshell_linewidth'] = {'linewidth_meV': onshell}
    dyn = {nu: float(np.interp(wq[nu], TB.get(nu, 1578.87)[1], TB.get(nu, 1578.87)[2])) for nu in (1, 20)}
    dref = solve_set(lambda nu, t: dyn[nu])
    dref['rbm_instability_temperature_K'] = tstar(lambda nu, t: dyn[nu])
    rob['dynamic_reference_at_sigma'] = dref
    # log interpolation / fit of the released static tables only
    rel_static = {nu: [TB.get(nu, t)[4] for t in RELEASED] for nu in (1, 20)}
    lt = np.log(RELEASED)
    s20 = np.interp(math.log(100.0), [math.log(50.0), math.log(700.0)], [rel_static[20][0], rel_static[20][1]])
    b, a = np.polyfit(lt, rel_static[1], 1)
    tlog = float(math.exp((pdc[1] - wq[1] / 2 - a) / b))
    rob['released_tables_log_fit'] = {'omega_adiabatic_100K_meV': math.sqrt(wq[20] ** 2 + 2 * wq[20] * (s20 - pdc[20])),
                                      'rbm_instability_temperature_K': tlog}
    a_, b_, d_, z_ = EQ9
    rob['eq9_threshold'] = {'threshold_coupling_meV': {str(nu): a_ + b_ * wq[nu] ** d_ / M_STAR ** z_ for nu in (1, 20)}}
    # interpolation of the archived Fig. 1 output: rows follow main()'s loop (a = 0.2, 0.5, 1, 2, 5; omega = 10 ... 200 meV)
    rows = [l for l in json.loads((M.POOL/'archive.json').read_text())['files']['Fig1/backup/output.csv'].split('\n') if l.strip()]
    last = [[float(x) for x in r.split(',') if x.strip() and 'nan' not in x][-1] for r in rows]
    a_row = [last[5 + j] for j in range(5)]                      # a = 0.5 (m* = 2): rows 5-9
    rob['fig1_output_interpolation'] = {'threshold_coupling_meV': {str(nu): float(np.interp(wq[nu], [10, 50, 100, 150, 200], a_row))
                                                                   for nu in (1, 20)}}

    # lab evaluator routes
    grid = np.linspace(110.0, 135.0, 41)

    def lab_route(nk, eta=None, processes=('k+q',), pairs=((1, 1), (1, 2), (2, 1), (2, 2)), with_tstar=True):
        sig_eta = None if eta is None else eta
        p20 = lab_pi(20, 1578.87, [0.0], nk, sig_eta, processes, pairs)[0][0]
        out = {'omega_nonadiabatic_meV': {}, 'linewidth_meV': {}}
        for t in TARGET_T:
            r, i = lab_pi(20, t, grid, nk, eta, processes, pairs)
            o, g = M.pole(grid, r, i, wq[20], p20)
            out['omega_nonadiabatic_meV'][table_name(t)], out['linewidth_meV'][table_name(t)] = o, g
        s = lab_pi(20, 100.0, [0.0], nk, eta, processes, pairs)[0][0]
        o2 = wq[20] ** 2 + 2 * wq[20] * (s - p20)
        out['omega_adiabatic_100K_meV'] = math.copysign(math.sqrt(abs(o2)), o2)
        if with_tstar:
            p1 = lab_pi(1, 1578.87, [0.0], min(nk, 1000000), sig_eta, processes, pairs)[0][0]
            ts = [50.0, 100.0, 200.0, 300.0, 400.0] + TSTAR_T + [700.0, 1000.0, 1578.87]
            fs = []
            for t in ts:
                r, i = lab_pi(1, t, [0.0, 0.5], min(nk, 1000000), eta, processes, pairs)
                fs.append(M.instability_function(wq[1], r[0], i[1] / 0.5, p1))
            out['rbm_instability_temperature_K'] = (max(ts) if all(f < 0 for f in fs) else M.crossing(ts, fs))
        return out

    rob['lab_evaluator_1e6'] = lab_route(1000000)
    rob['lab_evaluator_eta_2p5meV'] = {k: v for k, v in lab_route(1000000, eta=2.5, with_tstar=False).items()
                                       if k != 'omega_adiabatic_100K_meV'}

    # ---------------- ablations (realistic wrong paths; must fail) ----------------
    abl = {}
    abl['no_dfpt_subtraction'] = solve_set(lambda nu, t: 0.0)
    abl['no_dfpt_subtraction']['rbm_instability_temperature_K'] = tstar(lambda nu, t: 0.0)
    same = solve_set(lambda nu, t: TB.get(nu, t)[4])
    same['rbm_instability_temperature_K'] = tstar(lambda nu, t: TB.get(nu, t)[4])
    abl['same_temperature_reference'] = same
    t_ev = 0.01 / M.KB_EV

    def static_at(nu, t):
        lo, hi = max(x for x in TB.temps if x <= t), min(x for x in TB.temps if x >= t)
        return float(np.interp(math.log(t), [math.log(lo), math.log(hi)], [TB.get(nu, lo)[4], TB.get(nu, hi)[4]]))
    p_ev = {nu: static_at(nu, t_ev) for nu in (1, 20)}
    abl['smearing_read_as_eV'] = solve_set(lambda nu, t: p_ev[nu])
    abl['smearing_read_as_eV']['rbm_instability_temperature_K'] = tstar(lambda nu, t: p_ev[nu])
    ad300 = math.sqrt(wq[20] ** 2 + 2 * wq[20] * (TB.get(20, 300.0)[4] - pdc[20]))
    abl['adiabatic_for_nonadiabatic'] = {'omega_nonadiabatic_meV': {'100': ref['omega_adiabatic_100K_meV'], '300': ad300}}
    abl['nonadiabatic_for_adiabatic'] = {'omega_adiabatic_100K_meV': ref['omega_nonadiabatic_meV']['100']}
    hdr = {}
    for nu in (1, 20):
        w_, w, re, im, st = TB.get(nu, 1578.87)
        hdr[nu] = float(im[int(np.argmin(np.abs(w)))])        # header "Real Sigma(w=0)" = column 8 at w ~ 0
    abl['header_column_order'] = solve_set(lambda nu, t: hdr[nu])
    abl['header_column_order']['rbm_instability_temperature_K'] = tstar(lambda nu, t: hdr[nu])
    both = solve_set(std, scale=2.0)
    both['rbm_instability_temperature_K'] = tstar(std, scale=2.0)
    abl['both_processes_summed'] = both
    abl['sorted_intraband_only'] = lab_route(1000000, pairs=((1, 1), (2, 2)))
    coarse = lab_route(100000, with_tstar=False)
    abl['coarse_mesh_1e5'] = {'linewidth_meV': coarse['linewidth_meV'], 'omega_nonadiabatic_meV': coarse['omega_nonadiabatic_meV']}
    native = lab_route(2000, with_tstar=False)
    abl['native_archived_grid'] = {'omega_nonadiabatic_meV': native['omega_nonadiabatic_meV'],
                                   'linewidth_meV': native['linewidth_meV']}
    arch = json.loads((M.POOL/'archive.json').read_text())
    scan = arch['dfpt_scan']
    fd10 = {nu: scan['frequency_cm1'][scan['modes'].index(nu)][0] * M.CM1_MEV for nu in (1, 20)}
    fdq = solve_set(std, wq_of=lambda nu: fd10[nu])
    fdq['rbm_instability_temperature_K'] = tstar(std, wq_of=lambda nu: fd10[nu])
    abl['fd_scan_reference_frequency'] = fdq
    ts_scan = [s * 1e-3 * M.RY_EV / M.KB_EV for s in scan['sigma_mRy']][::-1]
    om_scan = [v * M.CM1_MEV for v in scan['frequency_cm1'][0]][::-1]
    abl['dfpt_scan_as_temperature_proxy'] = {'rbm_instability_temperature_K':
                                             M.crossing(ts_scan, [math.copysign(o * o, o) for o in om_scan])}
    abl['inverted_mass_parameter'] = {'threshold_coupling_meV': {str(nu): M.threshold_model(M_STAR, wq[nu]) for nu in (1, 20)}}
    # interpolation between the released tables (50 K and 700 K) instead of computing the withheld temperatures
    end = {}
    for t in (50.0, 700.0):
        w_, w, re, im, st = TB.get(20, t)
        end[t] = M.pole(w, re, im, w_, pdc[20])
    lin = lambda t, j: end[50.0][j] + (end[700.0][j] - end[50.0][j]) * (t - 50.0) / 650.0
    abl['interpolate_released_tables'] = {'omega_nonadiabatic_meV': {table_name(t): lin(t, 0) for t in TARGET_T},
                                          'linewidth_meV': {table_name(t): lin(t, 1) for t in TARGET_T},
                                          'rbm_instability_temperature_K': tstar(std, temps=[50.0, 700.0])}

    # ---------------- archived quantities quoted in the report / rubric ----------------
    z = np.load(M.POOL/'coupling.npz')
    gk = {}
    for nu in (1, 20):
        g = z[f'g_{nu}_p']
        gk[str(nu)] = {'|g_12| at k=1575 (k_F)': float(g[0, 1, 1575]), '|g_21| at k=1575': float(g[1, 0, 1575])}
    facts = {'pi_sigma_meV': pdc, 'omega_Q_meV': wq, 'fd_dfpt_10mRy_meV': fd10, 'T_sigma_K': M.T_SIGMA,
             'evaluator_max_deviation_meV': worst, 'g_at_fermi_points_meV': gk, 'header_path_pi_sigma_meV': hdr,
             'pi_sigma_if_eV_meV': p_ev}

    # ---------------- write task files ----------------
    (LAB/'data').mkdir(parents=True, exist_ok=True)
    shutil.copy2(HERE/'selfenergy.py', LAB/'selfenergy.py')
    shutil.copy2(M.POOL/'coupling.npz', LAB/'data/coupling.npz')
    text = json.loads((M.POOL/'pit_text.json').read_text())
    (LAB/'data/tables.json').write_text(json.dumps({table_name(t): text[table_name(t)] for t in RELEASED}))
    files = arch['files']
    (LAB/'data/archive.json').write_text(json.dumps({
        'bands': arch['bands'], 'phonons': arch['phonons'], 'dfpt_scan': scan, 'files': files,
        'fermi_shift_comment': arch['fermi_shift_comment'], 'k_plus_minus_q_comment': arch['k_plus_minus_q_comment']}))
    requests = {
        'archive_file': {'kind': 'archive_file', 'file': '<one of archive_files>', 'cost': 1},
        'bands': {'kind': 'bands', 'cost': 1},
        'phonons': {'kind': 'phonons', 'cost': 1},
        'coupling': {'kind': 'coupling', 'mode': '<1, 20 or 23>', 'cost': 1},
        'dfpt_scan': {'kind': 'dfpt_scan', 'cost': 1},
        'self_energy_table': {'kind': 'self_energy_table', 'mode': '<1, 20 or 23>', 'temperature_K': '<released temperature>', 'cost': 1},
        'self_energy': {'kind': 'self_energy', 'mode': '<1, 20 or 23>', 'temperature_K': '<30-2000>', 'omega_meV': '<list>',
                        'nk': '<2000-4000000>', 'eta_meV': '<(0, 50]>', 'processes': '<list>', 'band_pairs': '<list>',
                        'cost': '1 if nk <= 100000, 2 if nk <= 1000000, else 3'}}
    public = {'task_id': 'cnt-kohn-anomaly-nonadiabatic', 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
              'requests': requests, 'catalog': '/app/catalog.json'}
    (LAB/'public.json').write_text(json.dumps(public, indent=1) + '\n')
    (TASK/'tests/public.json').write_text(json.dumps({'task_id': 'cnt-kohn-anomaly-nonadiabatic', 'solver_budget': BUDGET,
                                                      'min_report_chars': 400}, indent=1) + '\n')
    file_desc = {
        'README.md': 'archive README',
        'dft_inputs/scf.in': 'ground-state input (Quantum ESPRESSO 7.3)',
        'dft_inputs/ph.in': 'DFPT input',
        'dft_inputs/epw.in': 'EPW input of the Wannier / electron-phonon run',
        'dft_inputs/specfun/epw.in': 'EPW input template of the phonon self-energy (spectral function) runs',
        'dft_inputs/specfun/input.py': 'script that generated the self-energy runs from the template',
        'dft_inputs/specfun/pathq.dat': 'q-point list of the self-energy runs',
        'Fig1/coupling_calculator.c': 'source of the analytic 1D model of Fig. 1 (threshold coupling)',
        'Fig1/fig1.py': 'plotting script of Fig. 1',
        'Fig1/backup/output.csv': 'archived output of coupling_calculator (Fig. 1 data)'}
    assert set(file_desc) == set(files)
    catalog = {
        'system': ('(3,3) single-wall carbon nanotube (diameter about 4 A): 12 C atoms per cell, period c = 2.4477 A along '
                   'the tube axis z, 13 A of vacuum; LDA, norm-conserving PseudoDojo pseudopotential, 100 Ry; Quantum '
                   'ESPRESSO 7.3 (DFPT) and EPW with two Wannier bands around the Fermi level.'),
        'source': ('E. Marazzi, S. Ponce, J.-C. Charlier, G.-M. Rignanese, Materials Cloud Archive record vs5fd-dbv98 '
                   '(DOI 10.24435/materialscloud:yx-x1, CC BY 4.0). All data are retrieved author results; the '
                   'self_energy evaluator recomputes self-energies from the archived coupling tables.'),
        'nesting_wavevector': ('All coupling and self-energy data refer to one phonon wavevector Q = 0.424 (units of 2 pi / c) '
                               'along the tube axis.'),
        'modes_at_Q': ('36 phonon branches, numbered at Q from the lowest frequency. Coupling data, archived self-energy '
                       'tables and the evaluator cover nu = 1 (radial breathing mode, RBM), nu = 20 and nu = 23 (in-plane '
                       'optical A1 modes).'),
        'units': {'energy': 'meV for phonon frequencies, self-energies and couplings; eV for electron energies',
                  'temperature': 'K', 'k_B': '8.617333262145179e-05 eV/K', 'Ry': '13.605693122994 eV',
                  'cm-1_to_meV': 0.12398419843320026},
        'budget': BUDGET,
        'purchase_costs': {k: v['cost'] for k, v in requests.items()},
        'requests': {
            'archive_file': {'format': {'kind': 'archive_file', 'file': '<name>'},
                             'returns': 'the verbatim text of one small archive file', 'archive_files': file_desc},
            'bands': {'format': {'kind': 'bands'},
                      'returns': ('DFT band structure from Fig2/e_band/band.eig: 12 bands at 5000 k along the tube axis '
                                  '(k in units of 2 pi / alat, alat = 13.0 A), eV, with the Fermi-level line of the '
                                  "authors' band-plot script (fig2_a_b.py) quoted verbatim")},
            'phonons': {'format': {'kind': 'phonons'},
                        'returns': ('DFPT phonon dispersion from Fig2/ph_band/matdyn.freq: 36 branches at 1001 q from '
                                    'Gamma along the tube axis (q in units of 2 pi / alat), cm-1')},
            'coupling': {'format': {'kind': 'coupling', 'mode': 20},
                         'returns': ('the archived EPW coupling tables Fig2/g/g+q.<n m nu> and g-q.<n m nu> for one mode: '
                                     'uniform 2000-point k grid (crystal units, -0.5 <= k < 0.5); energies (eV) of the two '
                                     'Wannier bands at k and at k+Q (or k-Q), ordered by energy at each k as in the '
                                     'archive; omega_Q (meV); |g_mn(k)| (meV) for the band pairs (n, m), n at k and m at '
                                     "k+-Q; plus a comment quoted from the authors' coupling-plot script")},
            'dfpt_scan': {'format': {'kind': 'dfpt_scan'},
                          'returns': ('DFPT frequencies at Q of nu = 1, 20, 23 computed with Fermi-Dirac smearing '
                                      'sigma = 10, 5, 3, 2, 1, 0.5 mRy, as tabulated in the authors\' Fig. 2 script '
                                      '(their "Eq. (1)" curves), cm-1')},
            'self_energy_table': {'format': {'kind': 'self_energy_table', 'mode': 20, 'temperature_K': 700},
                                  'released_temperatures_K': RELEASED,
                                  'returns': ('the archived EPW phonon self-energy file Fig2/PiT/<T> for one mode: its two '
                                              'header lines and its 801 data lines (frequencies w from -500 to 150 meV), '
                                              'verbatim. Tables at other temperatures exist in the archive but are '
                                              'withheld.')},
            'self_energy': {
                'format': {'kind': 'self_energy', 'mode': '1 | 20 | 23', 'temperature_K': 'number, 30-2000',
                           'omega_meV': 'list of 1-41 real frequencies in [-500, 500]', 'nk': 'integer, 2000-4000000',
                           'eta_meV': 'number in (0, 50]', 'processes': 'non-empty list from ["k+q", "k-q"]',
                           'band_pairs': 'non-empty list of distinct [n, m], n and m in {1, 2}'},
                'returns': 'Re and Im of Pi_nu(w, T) (meV) at the requested frequencies',
                'definition': ('For one process p (k -> k+Q or k -> k-Q): Pi(w, T) = (2/N_k) sum_k sum_(n,m) '
                               '|g_mn(k)|^2 [f(e_m(k+pQ)) - f(e_n(k))] / [e_m(k+pQ) - e_n(k) - w - i eta]; the result '
                               'is summed over the selected processes and band pairs. f is the Fermi-Dirac function at '
                               'temperature_K with the chemical potential fixed at -1.0380 eV (fermi_energy of the archived '
                               'EPW self-energy input). The archived 2000-point bands and |g|^2 are linearly interpolated '
                               'onto nk uniform k-points after the two bands have been followed through their crossings; '
                               'band_pairs use the archive\'s energy ordering (1 = lower band) at each k.'),
                'sign_convention': 'Re Pi < 0 softens the mode; Im Pi(w > 0) < 0.'}},
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    hidden = {'reference': ref, 'thresholds': THRESHOLDS}
    (LAB/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = (list(ref['omega_nonadiabatic_meV'].values()) + list(ref['linewidth_meV'].values()) +
              [ref['omega_adiabatic_100K_meV'], ref['rbm_instability_temperature_K']] +
              list(ref['threshold_coupling_meV'].values()))
    (LAB/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    (HERE/'calibration.json').write_text(json.dumps({'reference': ref, 'robustness': rob, 'ablations': abl, 'facts': facts,
                                                     'evaluator_check': dev}, indent=1, default=float) + '\n')
    print(json.dumps({'reference': ref, 'facts': facts}, indent=1, default=float))


if __name__ == '__main__':
    main()
