"""Write the staged, lab, private and ablation files of mnte-altermagnet-chiral-magnons.

Run in the runtime image after fetch.py and calib.py:
    MNTE_CACHE=/cache python3 -B make_task.py
"""
import glob
import json
import os
import numpy as np
import forward as F
import reduce as R
import sym

HERE = os.path.dirname(os.path.abspath(__file__))
TASK = os.path.join(HERE, '..', '..', 'tasks', 'mnte-altermagnet-chiral-magnons')
X = R.X
CAL = json.load(open(os.path.join(HERE, 'calibration.json')))
BUDGET = 8
COSTS = {'hyspec_tile': 1, 'hyspec_cut': 1, 'pnd': 1}
TILES = {'u_low_v_low': (slice(0, 22), slice(0, 32)), 'u_low_v_high': (slice(0, 22), slice(32, 63)),
         'u_high_v_low': (slice(22, 44), slice(0, 32)), 'u_high_v_high': (slice(22, 44), slice(32, 63))}
THRESHOLDS = {'dJ_window_meV': [0.055, 0.12], 'in20_mean_tol_qh': 0.012, 'in20_separation_rel_tol': 0.35,
              'q_score_min': 0.97}


def r6(x):
    return None if (x is None or not np.isfinite(x)) else float('%.6g' % x)


def load_dat(path):
    a = np.loadtxt(path)
    return a


def grid_map(path):
    a = load_dat(path)
    u = np.round(a[:, 2].reshape(44, 63)[:, 0], 4)
    v = np.round(a[:, 3].reshape(44, 63)[0], 4)
    I = a[:, 0].reshape(44, 63)
    E = a[:, 1].reshape(44, 63)
    return u, v, [[r6(x) for x in row] for row in I], [[r6(x) for x in row] for row in E]


def cut(path, axis):
    a = load_dat(path)
    return {'axis': axis, 'x': [r6(x) for x in a[:, 2]], 'intensity': [r6(x) for x in a[:, 0]], 'error': [r6(x) for x in a[:, 1]]}


def fig2b(path):
    rows = [l.split(',') for l in open(path).read().splitlines()[2:] if l.strip()]
    a = np.array(rows, float)
    return a


def header_only(path):
    lines = open(path, errors='replace').read().splitlines()
    out = []
    for l in lines:
        if l.startswith('DATA_:'):
            break
        out.append(l)
    return '\n'.join(out)


def main():
    u, v, sym_I, sym_E = grid_map(os.path.join(X, 'Figure2/Fig2c.dat'))
    _, _, chi_I, chi_E = grid_map(os.path.join(X, 'Figure2/Fig2e.dat'))
    U, V = np.meshgrid(u, v, indexing='ij')
    Egrid = F.hyspec_energy(U, V)
    gray, green = fig2b(os.path.join(X, 'Figure2/Fig2b_gray.dat')), fig2b(os.path.join(X, 'Figure2/Fig2b_green.dat'))
    assert np.allclose(gray[:, 0] - gray[0, 0], green[:, 0] - green[0, 0], atol=1e-3)
    centre = gray[np.argmax(gray[:, 1] + green[:, 1]), 0]
    lab = {}
    for qname, I, E, text in (('symmetric', sym_I, sym_E, 'symmetric S = (I_x^{-0} + I_x^{+0}) / 2'),
                              ('chiral', chi_I, chi_E, 'chiral M = (I_x^{-0} - I_x^{+0}) / (2 |P_i|)')):
        for tname, (su, sv) in TILES.items():
            lab['hyspec_tile:%s:%s' % (qname, tname)] = {
                'quantity': text, 'tile': tname, 'u': u[su].tolist(), 'v': v[sv].tolist(),
                'intensity': [row[sv] for row in I[su]], 'error': [row[sv] for row in E[su]],
                'layout': 'intensity[i][j] at (u[i], v[j]); null = no coverage'}
    lab.update({
        'hyspec_cut:g': dict(cut(os.path.join(X, 'Figure2/Fig2g.dat'), 'u'), quantity='symmetric', integrated='v from 0.40 to 0.50'),
        'hyspec_cut:h': dict(cut(os.path.join(X, 'Figure2/Fig2h.dat'), 'v'), quantity='symmetric', integrated='u from 0.95 to 1.05'),
        'hyspec_cut:i': dict(cut(os.path.join(X, 'Figure2/Fig2i.dat'), 'u'), quantity='chiral', integrated='v from 0.40 to 0.50'),
        'hyspec_cut:j': dict(cut(os.path.join(X, 'Figure2/Fig2j.dat'), 'v'), quantity='chiral', integrated='u from 0.95 to 1.05'),
        'pnd:HYSPEC': {'reflection_hkl': [1, 0, -1], 'frame': 'HYSPEC slice frame (same as the maps and cuts)',
                       'polarization_axis': 'y (Blume-Maleyev)', 'temperature_K': 20,
                       'delta_omega_deg': [r6(x - centre) for x in gray[:, 0]],
                       'I_y_plus0': [r6(x) for x in gray[:, 1]], 'I_y_plus0_error': [r6(x) for x in gray[:, 2]],
                       'I_y_minus0': [r6(x) for x in np.interp(gray[:, 0] - gray[0, 0], green[:, 0] - green[0, 0], green[:, 1])],
                       'I_y_minus0_error': [r6(x) for x in np.interp(gray[:, 0] - gray[0, 0], green[:, 0] - green[0, 0], green[:, 2])]},
        'pnd:IN20:A': {'reflection_hkl': [1, 0, 1], 'frame': 'IN20 instrument frame (same as the inelastic scan files)',
                       'file': '043731', 'raw_file': open(os.path.join(X, 'Figure3/Fig3b_red'), errors='replace').read()},
        'pnd:IN20:B': {'reflection_hkl': [1, 0, 1], 'frame': 'IN20 instrument frame (same as the inelastic scan files)',
                       'file': '043806', 'raw_file': open(os.path.join(X, 'Figure3/Fig3b_blue'), errors='replace').read()},
    })
    withheld = []
    for st, folder in (('A', 'Fig3d'), ('B', 'Fig3c')):
        for f in sorted(glob.glob(os.path.join(X, 'Figure3', folder, '*'))):
            withheld.append({'state': st, 'file': os.path.basename(f), 'header': header_only(f)})
    catalog = {
        'material': 'MnTe, hexagonal NiAs type, P6_3/mmc (No. 194); altermagnet with T_N = 307 K',
        'structure': {'a_A': sym.A_LAT, 'c_A': sym.C_LAT, 'Mn_fractional': {'A': [0, 0, 0], 'B': [0, 0, 0.5]},
                      'Te_fractional': [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]],
                      'cartesian_axes': 'x || a, z || c; a* lies in the ab plane, 30 deg from a (perpendicular to b)',
                      'reciprocal_lattice_units': 'r.l.u. of this hexagonal cell (a* = 4 pi / (sqrt(3) a), c* = 2 pi / c)'},
        'magnetic_structure': ('Collinear: all Mn at z = 0 (sublattice A) parallel, all Mn at z = 1/2 (sublattice B) antiparallel to A; '
                               'moments in the ab plane along an a*-type axis. The crystal has three 120-degree orientational domains '
                               'and, for each, two antiphase domains with reversed Neel vector n = S_A - S_B. Field cooling through T_N '
                               'changes the antiphase-domain balance.'),
        'prior_hamiltonian': {
            'source': 'Liu et al., Phys. Rev. Lett. 133, 156702 (2024): unpolarized inelastic neutron scattering, the Hamiltonian used for the model curves of the present data set',
            'convention': 'H = sum over Mn pairs (each pair once) J S_i.S_j + D sum_i (S_i^z)^2, S = 5/2, J > 0 antiferromagnetic, z || c, energies in meV',
            'shells': [{'distance_A': 3.3555, 'pairs': 'A-B (along c)', 'J_meV': 3.99},
                       {'distance_A': 4.148, 'pairs': 'same sublattice (in plane)', 'J_meV': -0.120},
                       {'distance_A': 5.3353, 'pairs': 'A-B', 'J_meV': 0.472},
                       {'distance_A': 9.8313, 'pairs': 'same sublattice, 12 bonds per Mn',
                        'J_meV': 'average over the 12 bonds = 0.0230; how this exchange is distributed over the bonds is not fixed by the unpolarized data'}],
            'D_meV': 0.0482},
        'hyspec': {
            'instrument': 'HYSPEC (SNS), time of flight, half-polarized: incident polarization P_i approximately parallel to Q (x), scattered polarization not analysed',
            'sample': 'S1, 20 K, after field cooling; crystal orientation fixed (omega = 32 deg), E_i = 37 meV',
            'coordinates': 'A slice pixel (u, v) is the wave vector Q = u (1, 0, -2.23) + v (1, 0, 1.6) in r.l.u. (scattering plane (h 0 l)). '
                           'At fixed crystal orientation each pixel has a single energy transfer, given by energy_grid_meV.',
            'energy_grid_meV': {'u': u.tolist(), 'v': v.tolist(), 'E': [[round(float(x), 2) for x in row] for row in Egrid],
                                'formula': 'E = E_i - 2.0721 |k_i - Q|^2 (meV, 1/A) with E_i = 37 meV, k_i in the (h 0 l) plane at 32 deg from a*, and the reduction lattice a = 4.100 A, c = 6.843 A'},
            'quantities': {'symmetric': '(I_x^{-0} + I_x^{+0}) / 2', 'chiral': '(I_x^{-0} - I_x^{+0}) / (2 |P_i|)',
                           'notation': 'I_alpha^{if}: polarization axis alpha (Blume-Maleyev: x || Q, y in the scattering plane, z = x cross y || b), incident i and scattered f spin states +, -, or 0 (not analysed)'},
            'cuts': {'g': 'symmetric, along u, integrated over v = 0.40-0.50', 'h': 'symmetric, along v, integrated over u = 0.95-1.05',
                     'i': 'chiral, along u, integrated over v = 0.40-0.50', 'j': 'chiral, along v, integrated over u = 0.95-1.05'},
            'pnd': 'Polarized diffraction (rocking scan) at the (1 0 -1) Bragg peak in the slice frame: I_y^{+0} and I_y^{-0}. (1 0 -1) has nuclear (Te) and magnetic (Mn) contributions.'},
        'in20': {
            'instrument': 'IN20 (ILL) triple-axis spectrometer with Cryopad, full polarization analysis',
            'sample': 'S2, a different crystal; states A and B were prepared by field cooling in opposite fields',
            'pnd': 'Raw ILL files of rocking scans at the (1 0 1) Bragg peak in the IN20 instrument frame, the same frame as the inelastic scan files, '
                   'with polarization along y and along z (non-spin-flip channels ++ and --).',
            'flipper_convention': 'In the raw files, I9 and I10 are the incident and scattered flipper currents: I9 > 0 selects incident spin + along P1, '
                                  'I10 > 0 selects scattered spin + along P2 (PX1..PZ2 give P1 and P2 in the Blume-Maleyev frame). '
                                  'So (I9, I10) = (+2, -2) is the channel I^{+-}.',
            'withheld_inelastic_scans': {'note': 'Constant-energy scans at 32 meV with P1 = P2 = x, measured in channels I_x^{+-} and I_x^{-+} for states A and B. '
                                                 'The lab does not serve them; their headers (everything above DATA_:) are given here.',
                                         'scans': withheld}},
        'requests': {'hyspec_tile': {'example': {'kind': 'hyspec_tile', 'quantity': 'symmetric', 'tile': 'u_low_v_low'},
                                     'quantity': ['symmetric', 'chiral'],
                                     'tile': {k: {'u_range': [float(u[su][0]), float(u[su][-1])], 'v_range': [float(v[sv][0]), float(v[sv][-1])]}
                                              for k, (su, sv) in TILES.items()},
                                     'cost': COSTS['hyspec_tile']},
                     'hyspec_cut': {'example': {'kind': 'hyspec_cut', 'cut': 'g'}, 'cut': ['g', 'h', 'i', 'j'], 'cost': COSTS['hyspec_cut']},
                     'pnd': {'example': {'kind': 'pnd', 'instrument': 'IN20', 'state': 'A'},
                             'options': [{'instrument': 'HYSPEC'}, {'instrument': 'IN20', 'state': 'A'}, {'instrument': 'IN20', 'state': 'B'}],
                             'cost': COSTS['pnd']}},
        'budget': BUDGET,
        'units': 'HYSPEC intensities are in the archive units (arbitrary, common to maps and cuts); IN20 counts with monitor M1.'}
    truth = CAL['in20_truth']
    comb = truth['combined']
    hidden = {'in20': {'low_qh': comb['centres'][0], 'high_qh': comb['centres'][1], 'mean_qh': comb['mean'],
                       'separation_qh': comb['separation'],
                       'chirality': {st: {'low': truth[st]['signs'][0], 'high': truth[st]['signs'][1]} for st in ('A', 'B')},
                       'significance_sigma': {st: truth[st]['significance'] for st in ('A', 'B')}},
              'reference': {'dJ_meV': CAL['dj_ref'], 'q_max_rlu': [1 / 3, 0, 0.25]},
              'thresholds': THRESHOLDS}
    values = [comb['centres'][0], comb['centres'][1], comb['mean'], comb['separation'], CAL['dj_ref']]
    P = CAL['in20_predictions']
    pk = lambda name: {'low': P[name]['centres'][0], 'high': P[name]['centres'][1]}
    abl = {'no_split': {'dJ': 0.0, 'peaks': pk('no_split')},
           'chiral_dipole_full_chirality': {'dJ': CAL['dipole']['dJ'], 'peaks': {'low': CAL['dipole']['peaks_nominal_l'][0], 'high': CAL['dipole']['peaks_nominal_l'][1]},
                                            'peaks_true_l': {'low': CAL['dipole']['peaks_true_l'][0], 'high': CAL['dipole']['peaks_true_l'][1]}},
           'nominal_l': {'peaks': pk('nominal_l')}, 'nominal_l_h_scaled': {'peaks': pk('nominal_l_h_scaled')},
           'half_dj': {'dJ': P['half_dj']['dJ'], 'peaks': pk('half_dj')},
           'q_gamma': [0.0, 0.0, 0.0], 'q_measured_line': [0.375, 0.0, 0.22]}
    rob = {'h_scaled': {'peaks': pk('h_scaled')}, 'dj_g_only': {'dJ': P['dj_g_only']['dJ'], 'peaks': pk('dj_g_only')},
           'dj_h_only': {'dJ': P['dj_h_only']['dJ'], 'peaks': pk('dj_h_only')},
           'dj_published': {'dJ': 0.0902, 'peaks': pk('dj_published')}, 'q_equivalent': [2 / 3, 0.0, 0.75]}
    w = lambda p, obj: open(os.path.join(TASK, p), 'w').write(json.dumps(obj, indent=1) + '\n')
    w('environment/agent/catalog.json', catalog)
    open(os.path.join(TASK, 'environment/lab/data.json'), 'w').write(json.dumps(lab) + '\n')
    w('environment/lab/private/hidden.json', hidden)
    open(os.path.join(TASK, 'environment/lab/private/hidden_values.json'), 'w').write(json.dumps(values) + '\n')
    w('solution/ablation_values.json', {'ablations': abl, 'robustness': rob})
    pub = {'task_id': 'mnte-altermagnet-chiral-magnons', 'solver_budget': BUDGET, 'reserved_validation_charge': 0}
    w('environment/lab/public.json', pub)
    w('tests/public.json', {'task_id': 'mnte-altermagnet-chiral-magnons', 'solver_budget': BUDGET, 'min_report_chars': 400})
    print(json.dumps({'hidden': hidden, 'values': values}, indent=1))


if __name__ == '__main__':
    main()
