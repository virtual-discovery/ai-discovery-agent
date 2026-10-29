"""Generate the staged, lab, private and ablation data for mos2-doping-dome-cdw from pool.json.

    python3 -B make_task.py          (run in the runtime image from this folder; about 3 minutes)

Writes environment/agent/catalog.json, environment/lab/{data.json,eliashberg.py}, private/hidden.json,
private/hidden_values.json, solution/ablation_values.json, the budget in both public.json files, and
reference.json here (all reference, ablation and robustness numbers, for write_design.py).
"""
import json
import math
import shutil
from pathlib import Path
import numpy as np
import eliashberg as E

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mos2-doping-dome-cdw'
POOL = json.loads((HERE/'pool.json').read_text())
MU = 0.13
BUDGET = 34
COSTS = {'spectrum': 1, 'dos': 1, 'eliashberg': 2, 'eliashberg_custom': 2, 'moments_table': 1, 'pes': 1,
         'soft_mode': 1, 'mo_displacement': 1, 'archived_tc_table': 1, 'a1g_raman': 1, 'experiment_dome': 1}
HIDDEN_X = ['0.15', '0.17']
RELEASED_X = ['0.12', '0.13', '0.14', '0.16', '0.18', '0.19']
SOC_X = ['0.18', '0.19']
CDW = ['dop1', 'dop2', 'dop3']
FLAVORS = ['A', 'NA', 'FA']
A_BOHR = float(POOL['qe_inputs']['1x1H']['celldm(1)'])
BOHR_A = 0.529177210903
AREA_A2 = math.sqrt(3) / 2 * (A_BOHR * BOHR_A) ** 2
PER_1E14 = 1 / (AREA_A2 * 1e-16) / 1e14          # 10^14 cm^-2 per (electron per 1x1 cell)
RY_MEV = 13.605693122994e3
THRESHOLDS = {'tc_1x1_rel': 0.03, 'critical_doping_abs': 0.002, 'tc_soc_rel': 0.03, 'tc_cdw_rel': 0.04}
CAPTION_1E14 = {'dop1': 3.0, 'dop2': 4.3, 'dop3': 4.6}   # paper, Fig. 4 caption (n1, n2, n3)


def spectrum(key):
    rows = np.array(POOL['spectra'][key]['rows'])
    return E.Spectrum(rows[:, 0], rows[:, 1])


def dos(key, n=None):
    d = POOL['dos'][key]
    return E.DOS(d['energy_eV'], d['total'], float(key.split('/')[1]) if n is None else n)


def tc(key, mode, n_override=None, dos_key=None, spec=None):
    s = spec or spectrum(key)
    d = dos(dos_key or '/'.join(key.split('/')[:2]), n_override) if mode == 'full' else None
    return E.critical_temperature(s, MU, d)[0]


def archived(col, soc=False):
    t = POOL['tc_tables_fig2e']['tc_013']
    i = t['columns'].index(col)
    return {'%.2f' % r[1]: r[i] for r in t['rows'] if bool(r[0]) == soc}


def lerp(x, x0, y0, x1, y1):
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def pes_delta(label):
    e = POOL['pes'][label]['energies_Ry']
    return (e[-1] - e[0]) * RY_MEV


def pes_well_delta(label):
    """CDW-well minimum from a cubic fit to the last five path points (may lie slightly beyond the endpoint)."""
    e = np.array(POOL['pes'][label]['energies_Ry'])
    s = np.linspace(0, 1, len(e))
    p = np.polyfit(s[-5:], e[-5:], 3)
    grid = np.linspace(0.6, 1.3, 7001)
    return (np.polyval(p, grid).min() - e[0]) * RY_MEV


def crossing(fn, quadratic=False):
    xs = [0.74, 0.76, 0.77, 0.78]
    d = [fn('pes_%03d' % round(100 * x)) for x in xs]
    if quadratic:
        c = np.polyfit(xs[:3], d[:3], 2)
        roots = [r.real for r in np.roots(c) if abs(r.imag) < 1e-12 and 0.74 <= r.real <= 0.77]
        return roots[0] / 4
    for i in range(3):
        if d[i] > 0 >= d[i + 1]:
            return lerp(0, d[i], xs[i], d[i + 1], xs[i + 1]) / 4
    raise RuntimeError('no crossing')


def main():
    shutil.copy2(HERE/'eliashberg.py', TASK/'environment/lab/eliashberg.py')
    DA, AD, CA, MM = archived('DA'), archived('AD'), archived('CA'), archived('MM')
    DAs, ADs, CAs, MMs = archived('DA', True), archived('AD', True), archived('CA', True), archived('MM', True)
    # --- solver reproduces the archived hidden-doping values (asserted) ---
    for x in HIDDEN_X + ['0.14', '0.16', '0.18']:
        assert abs(tc(f'1x1H/{x}/FA', 'full') / DA[x] - 1) < 2e-3, x
    # --- critical doping from the PES (energy of the 2x2 CDW endpoint relative to 1x1 H) ---
    deltas = {k: pes_delta(k) for k in POOL['pes']}
    xc = crossing(pes_delta)
    xc_quad = crossing(pes_delta, quadratic=True)
    xc_well = crossing(pes_well_delta)
    soft = POOL['fig3']['1x1_w_Gamma.txt']['rows']
    i = next(k for k, r in enumerate(soft) if r[1] < 0)
    x_soft = (soft[i - 1][0] + soft[i][0]) / 2 / PER_1E14
    hyst = POOL['fig3']['Mo_hysteresis_cdw.txt']['rows']
    j = next(k for k, r in enumerate(hyst) if r[1] > 0.005)
    x_spinodal = hyst[j][0] / PER_1E14
    soc_at = lambda x, table: lerp(x, 0.18, table['0.18'], 0.19, table['0.19'])
    tc_soc = soc_at(xc, DAs)
    # --- high-doping CDW: Tc in the constant-DOS limit, mapping of the archive labels via the authors' table ---
    cdw = {}
    for d in CDW:
        s = spectrum(f'2x2CDW-high/{d}/FA')
        cdw[d] = {'CA': tc(f'2x2CDW-high/{d}/FA', 'constant'), 'AD': E.allen_dynes(s.lam, s.wlog, MU, s.w2),
                  'MM': E.allen_dynes(s.lam, s.wlog, MU, correct=False)}
    table = {round(r[0], 2): r[3] for r in POOL['fig5']['tcs_CDW_800K.txt']['rows']}
    mapping = {d: min(table, key=lambda q: abs(table[q] - cdw[d]['AD'])) for d in CDW}
    assert mapping == {'dop1': 1.4, 'dop2': 1.5, 'dop3': 1.6}, mapping
    assert all(abs(cdw[d]['AD'] / table[mapping[d]] - 1) < 2e-3 for d in CDW)
    x_cdw = {d: mapping[d] / 4 for d in CDW}
    cdw_targets = {'0.35': cdw['dop1']['CA'], '0.40': cdw['dop3']['CA']}
    hidden = {'targets': {'tc_1x1_K': {x: DA[x] for x in HIDDEN_X}, 'critical_doping_e_per_cell': xc,
                          'tc_1x1_soc_at_critical_doping_K': tc_soc, 'tc_cdw_high_K': cdw_targets},
              'thresholds': THRESHOLDS,
              'provenance': {'tc_1x1_K': 'archived fig2e/tc_013.dat, column DA (Eliashberg, full DOS, FA, mu*=0.13)',
                             'critical_doping': 'linear interpolation of E(2x2 CDW endpoint) - E(1x1 H) of fig3/pes_074 and pes_076, divided by 4',
                             'tc_soc': 'archived fig2e/tc_013.dat SOC rows (DA) at 0.18 and 0.19, interpolated linearly to the critical doping',
                             'tc_cdw_high': 'eliashberg.py constant-DOS Tc of the archived fig4 2x2cdw_dop1/dop3 FA spectra; dop1 = 1.4 and dop3 = 1.6 e per 2x2 cell from the archived fig5/tcs_CDW_800K.txt'}}
    values = list(hidden['targets']['tc_1x1_K'].values()) + [xc, tc_soc] + list(cdw_targets.values())
    # --- ablations: realistic wrong paths, computed from the archive ---
    A = {}
    A['tc_linear_interpolation'] = {'0.15': (DA['0.14'] + DA['0.16']) / 2, '0.17': (DA['0.16'] + DA['0.18']) / 2}
    A['tc_unscaled_average_spectrum'] = {}
    for x, (lo, hi) in {'0.15': ('0.14', '0.16'), '0.17': ('0.16', '0.18')}.items():
        r0, r1 = np.array(POOL['spectra'][f'1x1H/{lo}/FA']['rows']), np.array(POOL['spectra'][f'1x1H/{hi}/FA']['rows'])
        avg = E.Spectrum(r0[:, 0], (r0[:, 1] + r1[:, 1]) / 2)
        A['tc_unscaled_average_spectrum'][x] = tc(None, 'full', float(x), f'1x1H/{lo}', avg)
    A['tc_NA_spectra'] = {x: tc(f'1x1H/{x}/NA', 'full') for x in HIDDEN_X}
    A['tc_A_spectra'] = {x: tc(f'1x1H/{x}/A', 'full') for x in HIDDEN_X}
    na_soc = {x: tc(f'1x1H-SOC/{x}/NA', 'full') for x in SOC_X}
    A['soc_NA_spectra'] = soc_at(xc, na_soc)
    A['tc_constant_dos'] = {x: CA[x] for x in HIDDEN_X}
    A['soc_constant_dos'] = soc_at(xc, CAs)
    A['tc_allen_dynes'] = {x: AD[x] for x in HIDDEN_X}
    A['soc_allen_dynes'] = soc_at(xc, ADs)
    A['cdw_allen_dynes'] = {'0.35': cdw['dop1']['AD'], '0.40': cdw['dop3']['AD']}
    A['tc_mcmillan'] = {x: MM[x] for x in HIDDEN_X}
    A['soc_mcmillan'] = soc_at(xc, MMs)
    A['cdw_mcmillan'] = {'0.35': cdw['dop1']['MM'], '0.40': cdw['dop3']['MM']}
    A['xc_soft_mode_instability'] = x_soft
    A['soc_at_soft_mode_instability'] = soc_at(x_soft, DAs)
    A['xc_cdw_spinodal'] = x_spinodal
    A['soc_at_cdw_spinodal'] = soc_at(x_spinodal, DAs)
    A['xc_per_2x2_cell'] = 4 * xc
    A['xc_areal_density_1e14'] = xc * PER_1E14
    nosoc = {x: DA[x] for x in ('0.18', '0.19')}
    A['soc_omitted'] = soc_at(xc, nosoc)
    cap = {d: CAPTION_1E14[d] / PER_1E14 for d in CDW}
    A['cdw_caption_dopings'] = {'0.35': lerp(0.35, cap['dop1'], cdw['dop1']['CA'], cap['dop2'], cdw['dop2']['CA']),
                                '0.40': lerp(0.40, cap['dop2'], cdw['dop2']['CA'], cap['dop3'], cdw['dop3']['CA'])}
    A['cdw_table_order'] = {'0.35': cdw['dop2']['CA'],
                            '0.40': lerp(0.40, 0.35, cdw['dop2']['CA'], 0.375, cdw['dop3']['CA'])}
    t1 = {'%.2f' % r[0]: r[3] for r in POOL['fig5']['tcs_800K.txt']['rows']}
    A['tc_archived_table'] = {'0.15': (t1['0.14'] + t1['0.16']) / 2, '0.17': (t1['0.16'] + t1['0.18']) / 2}
    A['cdw_archived_table'] = {'0.35': table[1.4], '0.40': table[1.6]}
    # --- robustness: defensible alternative estimates ---
    R = {}
    lam_tab = {'%.2f' % r[0]: r[3] for r in POOL['moments']['lambda']['rows']}
    wlog_tab = {'%.2f' % r[0]: r[3] * 1e-3 for r in POOL['moments']['omega_log_meV']['rows']}
    R['tc_custom_rigid_band'] = {}
    R['tc_custom_rigid_band_upper_dos'] = {}
    R['tc_moments_no_shape'] = {}
    R['tc_moments_ratio'] = {}
    for x, (lo, hi) in {'0.15': ('0.14', '0.16'), '0.17': ('0.16', '0.18')}.items():
        r0, r1 = np.array(POOL['spectra'][f'1x1H/{lo}/FA']['rows']), np.array(POOL['spectra'][f'1x1H/{hi}/FA']['rows'])
        avg = E.Spectrum(r0[:, 0], (r0[:, 1] + r1[:, 1]) / 2)
        scaled = E.Spectrum(r0[:, 0], (r0[:, 1] + r1[:, 1]) / 2 * lam_tab[x] / avg.lam)
        R['tc_custom_rigid_band'][x] = tc(None, 'full', float(x), f'1x1H/{lo}', scaled)
        R['tc_custom_rigid_band_upper_dos'][x] = tc(None, 'full', float(x), f'1x1H/{hi}', scaled)
        s0, s1 = spectrum(f'1x1H/{lo}/FA'), spectrum(f'1x1H/{hi}/FA')
        ad = lambda lam, wl, w2: E.allen_dynes(lam, wl, MU, w2)
        ratio = (DA[lo] / ad(s0.lam, s0.wlog, s0.w2) + DA[hi] / ad(s1.lam, s1.wlog, s1.w2)) / 2
        R['tc_moments_ratio'][x] = ad(lam_tab[x], wlog_tab[x], (s0.w2 + s1.w2) / 2) * ratio
        ratio0 = (DA[lo] / ad(lam_tab[lo], wlog_tab[lo], None) + DA[hi] / ad(lam_tab[hi], wlog_tab[hi], None)) / 2
        R['tc_moments_no_shape'][x] = ad(lam_tab[x], wlog_tab[x], None) * ratio0
    # lambda-weighted mix of the neighbours' spectra (matches the tabulated lambda without rescaling the shape)
    R['tc_lambda_weighted_mix'] = {}
    G = {'single_neighbour_shape_lower': {}, 'single_neighbour_shape_upper': {}}
    for x, (lo, hi) in {'0.15': ('0.14', '0.16'), '0.17': ('0.16', '0.18')}.items():
        r0, r1 = np.array(POOL['spectra'][f'1x1H/{lo}/FA']['rows']), np.array(POOL['spectra'][f'1x1H/{hi}/FA']['rows'])
        s0, s1 = E.Spectrum(r0[:, 0], r0[:, 1]), E.Spectrum(r1[:, 0], r1[:, 1])
        t = (lam_tab[x] - s0.lam) / (s1.lam - s0.lam)
        R['tc_lambda_weighted_mix'][x] = tc(None, 'full', float(x), f'1x1H/{lo}', E.Spectrum(r0[:, 0], (1 - t) * r0[:, 1] + t * r1[:, 1]))
        for name, r, d in (('single_neighbour_shape_lower', r0, lo), ('single_neighbour_shape_upper', r1, hi)):
            s = E.Spectrum(r[:, 0], r[:, 1])
            G[name][x] = tc(None, 'full', float(x), f'1x1H/{d}', E.Spectrum(r[:, 0], r[:, 1] * lam_tab[x] / s.lam))
    R['xc_quadratic'] = xc_quad
    R['soc_at_xc_quadratic'] = soc_at(xc_quad, DAs)
    R['xc_cdw_well_minimum'] = xc_well
    R['soc_at_xc_well_minimum'] = soc_at(xc_well, DAs)
    R['soc_ratio_method'] = lerp(xc, 0.18, DAs['0.18'] / DA['0.18'], 0.19, DAs['0.19'] / DA['0.19']) * \
        float(np.polyval(np.polyfit([0.17, 0.18, 0.19], [DA['0.17'], DA['0.18'], DA['0.19']], 2), xc))
    # --- lab data (released items only) and catalog ---
    released = {}
    for key in POOL['spectra']:
        phase, label, fl = key.split('/')
        if phase == '1x1H' and label in HIDDEN_X:
            continue
        released[key] = POOL['spectra'][key]
    dos_rel = {k: v for k, v in POOL['dos'].items() if k.split('/')[1] not in HIDDEN_X or k.startswith('1x1H-SOC')}
    tcs1 = POOL['fig5']['tcs_800K.txt']
    tables = {'1x1H': {'source_file': 'fig5/tcs_800K.txt', 'header': tcs1['header'],
                       'rows': [r for r in tcs1['rows'] if '%.2f' % r[0] not in HIDDEN_X],
                       'note': 'rows at 0.15 and 0.17 withheld'},
              '2x2CDW': {'source_file': 'fig5/tcs_CDW_800K.txt', 'header': POOL['fig5']['tcs_CDW_800K.txt']['header'],
                         'rows': POOL['fig5']['tcs_CDW_800K.txt']['rows']},
              '1Tprime': {'source_file': 'fig5/tcs_1tprime.txt', 'header': POOL['fig5']['tcs_1tprime.txt']['header'],
                          'rows': POOL['fig5']['tcs_1tprime.txt']['rows']}}
    data = {'costs': COSTS,
            'spectra': {k: {'source_file': v['file'], 'header': v['header'], 'rows': v['rows']} for k, v in released.items()},
            'dos': {k: {'source_files': v['files'], 'energy_eV': v['energy_eV'], 'total': v['total'],
                        'K_dz2_weighted': v['K_dz2'], 'remainder_Q': v['Q_rest']} for k, v in dos_rel.items()},
            'moments': POOL['moments'],
            'pes': {k: {'source_file': f'fig3/{k}', 'energies_Ry_per_cell': v['energies_Ry']} for k, v in POOL['pes'].items()},
            'soft_mode': {'1x1H': POOL['fig3']['1x1_w_Gamma.txt'], '2x2CDW': POOL['fig3']['2x2_w_Gamma.txt']},
            'mo_displacement': {'1x1H': POOL['fig3']['Mo_hysteresis_1x1.txt'], '2x2CDW': POOL['fig3']['Mo_hysteresis_cdw.txt']},
            'archived_tc_tables': tables,
            'a1g': {'theory': POOL['fig5']['a1g.txt'], 'experiment': POOL['fig5']['exp_a1g.txt']},
            'experiment_dome': POOL['fig5']['exp_dome_1']}
    catalog = make_catalog(released, dos_rel)
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': R}, indent=1) + '\n')
    for pub in (TASK/'environment/lab/public.json', TASK/'tests/public.json'):
        p = json.loads(pub.read_text())
        p['solver_budget'] = BUDGET
        pub.write_text(json.dumps(p, indent=1) + '\n')
    ref = {'targets': hidden['targets'], 'thresholds': THRESHOLDS, 'ablations': A, 'robustness': R, 'gray_zone': G,
           'pes_delta_meV_per_2x2_cell': deltas, 'xc_variants': {'linear': xc, 'quadratic': xc_quad, 'well_minimum': xc_well,
                                                                  'soft_mode': x_soft, 'cdw_spinodal': x_spinodal},
           'cdw_high': cdw, 'cdw_mapping_e_per_2x2': mapping, 'per_1e14': PER_1E14,
           'archived_tc013': {'DA': DA, 'AD': AD, 'CA': CA, 'MM': MM, 'DA_soc': DAs, 'AD_soc': ADs, 'CA_soc': CAs}}
    (HERE/'reference.json').write_text(json.dumps(ref, indent=1) + '\n')
    print(json.dumps({'targets': hidden['targets'], 'xc_variants': ref['xc_variants'], 'cdw': cdw}, indent=1))


def make_catalog(released, dos_rel):
    qe = POOL['qe_inputs']
    spectra = []
    for key, v in released.items():
        phase, label, fl = key.split('/')
        spectra.append({'phase': phase, 'label': label, 'flavor': fl, 'archive_file': v['file']})
    return {
        'material': 'Electron-doped MoS2 monolayer (1x1 H, 2x2 CDW and 1T\' phases). Released DFT/DFPT/EPW data of '
                    'Girotto Erhardt, Berges, Ponce and Novko (2025), Materials Cloud record 51yyy-ksf72.',
        'cells_and_doping': {
            'convention': 'Doping is added as a uniform background charge: Quantum ESPRESSO tot_charge = -q means q extra '
                          'electrons per simulation cell. Different archive items use different simulation cells.',
            'qe_system_blocks': {'1x1H (input/1x1H/sc.in)': qe['1x1H'], '2x2 CDW, low doping (input/2x2CDW_lowdop/sc.in)': qe['2x2CDW_lowdop'],
                                 '2x2 CDW, high doping (input/2x2CDW_highdop/sc.in)': qe['2x2CDW_highdop'],
                                 '1T\' (input/1T_prime/sc.in)': qe['1T_prime']},
            'lattice': 'ibrav = 4 (hexagonal); celldm(1) is the in-plane lattice constant in bohr (1 bohr = 0.529177 A). '
                       'The 2x2 cells repeat the 1x1 cell twice along a1 and a2; the 1T\' cell holds two formula units.'},
        'requests': {
            'spectrum': {'example': {'kind': 'spectrum', 'phase': '1x1H', 'label': '0.14', 'flavor': 'FA'}, 'cost': COSTS['spectrum'],
                         'returns': 'Eliashberg spectral function: omega (meV), alpha^2F, cumulative lambda(omega).'},
            'dos': {'example': {'kind': 'dos', 'phase': '1x1H', 'label': '0.14'}, 'cost': COSTS['dos'],
                    'returns': 'Electronic DOS of the two lowest conduction bands (states per eV, per spin and 1x1 cell) on a 1 meV grid: '
                               'total, the dz2-weighted part (K valleys) and the remainder (Q valleys).'},
            'eliashberg': {'example': {'kind': 'eliashberg', 'phase': '1x1H', 'label': '0.14', 'flavor': 'FA', 'dos': 'full', 'mu_star': 0.13},
                           'cost': COSTS['eliashberg'],
                           'returns': 'Tc (K) from the linearized isotropic Migdal-Eliashberg equations for an archived spectrum, '
                                      'with lambda, omega_log and omega_2 of the spectrum and solver diagnostics. "dos": "full" keeps '
                                      'the energy-dependent DOS (1x1H and 1x1H-SOC only), "constant" uses the constant-DOS limit.'},
            'eliashberg_custom': {'example': {'kind': 'eliashberg_custom', 'spectrum_meV': [[0.15, 0.0003], [0.3, 0.0007], ['...', '...']],
                                              'dos': {'phase': '1x1H', 'label': '0.14'}, 'electrons_per_cell': 0.14, 'mu_star': 0.13},
                                  'cost': COSTS['eliashberg_custom'],
                                  'returns': 'The same solver for a spectrum you supply (3-600 rows of [omega_meV, alpha^2F], omega increasing '
                                             'in (0, 100] meV). "dos" is "constant" or an archived DOS; with an archived DOS, '
                                             'electrons_per_cell sets the electron count in those bands (rigid band).'},
            'moments_table': {'example': {'kind': 'moments_table'}, 'cost': COSTS['moments_table'],
                              'returns': 'The authors\' table of lambda and omega_log (A, NA, FA) for the 1x1 H phase at 0.11-0.19 '
                                         '(archive fig2a-d/lambda_2d.txt and omegalog_2d.txt; paper Fig. 2d).'},
            'pes': {'example': {'kind': 'pes', 'label': 'pes_074'}, 'cost': COSTS['pes'],
                    'returns': 'Total energies (Ry per simulation cell) at 11 equally spaced points of the paper\'s Fig. 3(c) '
                               'path, from the 1x1 H structure (first point) to the 2x2 CDW structure (last point), in the '
                               '2x2 supercell. Labels are the archive file names.'},
            'soft_mode': {'example': {'kind': 'soft_mode', 'branch': '1x1H'}, 'cost': COSTS['soft_mode'],
                          'returns': 'Frequency (meV) of the Gamma* mode of the 2x2 supercell (folded from M of the 1x1 cell) along '
                                     'the 1x1 H or the 2x2 CDW branch, versus charge concentration n in 10^14 cm^-2 (paper Fig. 3e; '
                                     'negative values denote imaginary frequencies).'},
            'mo_displacement': {'example': {'kind': 'mo_displacement', 'branch': '2x2CDW'}, 'cost': COSTS['mo_displacement'],
                                'returns': 'Mo displacement (A) after relaxing a 2x2 supercell started from the 1x1 H or the CDW structure, '
                                           'versus n in 10^14 cm^-2 (paper Fig. 3d).'},
            'archived_tc_table': {'example': {'kind': 'archived_tc_table', 'phase': '2x2CDW'}, 'cost': COSTS['archived_tc_table'],
                                  'returns': 'The authors\' Tc table of a phase (archive fig5/), as archived: its header names the doping '
                                             'column and three Tc columns (adiabatic, NA, fully NA). The archive does not state how these '
                                             'Tc were computed. Phases: 1x1H, 2x2CDW, 1Tprime.'},
            'a1g_raman': {'example': {'kind': 'a1g_raman'}, 'cost': COSTS['a1g_raman'],
                          'returns': 'A1g frequency versus doping across phases (adiabatic, nonadiabatic; meV) and the Raman data on '
                                     'LixMoS2 (Li per MoS2, cm^-1) shown in the paper\'s Fig. 5(a).'},
            'experiment_dome': {'example': {'kind': 'experiment_dome'}, 'cost': COSTS['experiment_dome'],
                                'returns': 'Measured Tc (K) of gated MoS2 versus carrier density (10^14 cm^-2) as plotted in Fig. 5(b).'}},
        'spectra_available': spectra,
        'dos_available': [{'phase': k.split('/')[0], 'label': k.split('/')[1]} for k in dos_rel],
        'pes_available': sorted(POOL['pes']),
        'labels': {
            '1x1H': 'label = tot_charge per 1x1 cell (archive directories fig2e/data/dop_0XX)',
            '1x1H-SOC': 'fully relativistic spectra and DOS at 0.18 and 0.19 (archive dop_018_soc, dop_019_soc); both use one '
                        'Wannier Hamiltonian, with the electron count set by the doping',
            '2x2CDW-high': 'archive labels dop1, dop2, dop3 (fig4/2x2cdw_dopN_*); the archive gives no doping for them. The '
                           'paper\'s Fig. 4 caption: "the largest doping values considered in this work; n1 = 3.0 x 10^14 cm^-2, '
                           'n2 = 4.3 x 10^14 cm^-2, n3 = 4.6 x 10^14 cm^-2".',
            '1Tprime': 'label "archived" (fig4/1tprime_*); the paper\'s Fig. 4(a) quotes n = 2.62 x 10^14 cm^-2',
            'pes': 'archive files fig3/pes_074, pes_076, pes_077, pes_078; the colour scale of Fig. 3(c) spans n = 2.1-2.2 x 10^14 cm^-2'},
        'withheld': 'The 1x1 H spectra, DOS and archived Tc of the 0.15 and 0.17 calculations are withheld. The SOC calculations exist only at 0.18 and 0.19. No DOS is archived for the CDW and 1T\' phases.',
        'methods': {
            'dft': 'Quantum ESPRESSO, PBE, fully relativistic norm-conserving PseudoDojo potentials, 100 Ry, Coulomb truncation '
                   '(assume_isolated = 2D), Fermi-Dirac smearing 0.005 Ry (Tscf = 800 K), 24x24 k (1x1), 12x12 q.',
            'epc': 'EPW, 11 Wannier functions, fine grids 360x360 k and 120x120 q (1x1); alpha^2F at the 800 K electronic temperature.',
            'flavors': {'A': 'adiabatic (static DFPT phonons)', 'NA': 'nonadiabatic phonon frequency shifts',
                        'FA': 'fully nonadiabatic: frequency shifts plus phonon linewidths (Girotto and Novko, PRB 107, 064310 (2023))'},
            'dos': 'From the authors\' Wannier Hamiltonians on a 313x313 grid (fig2e/dos.py); energies relative to the Fermi level '
                   'of the doped cell at the 800 K smearing.',
            'eliashberg_solver': 'Linearized isotropic Migdal-Eliashberg equations on the Matsubara axis (paper Eqs. 1-9). Full DOS: '
                                 'energy-dependent N(eps) over the archived energy window, with the energy shift chi and conservation of '
                                 'the electron count, mu* referenced to omega_2 and rescaled to the Matsubara cutoff (20 omega_2). '
                                 'Constant DOS: mu*_C = mu* / (1 + mu* ln(omega_2/omega_C)). Tc below 1 K is reported as such. '
                                 'The solver reproduces the authors\' archived Eliashberg tables.'},
        'purchase_costs': COSTS, 'budget': BUDGET}


if __name__ == '__main__':
    main()
