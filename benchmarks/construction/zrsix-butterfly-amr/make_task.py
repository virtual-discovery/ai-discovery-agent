"""Write the catalog, the private targets, the ablation and robustness values and calibration.json for
zrsix-butterfly-amr. Every number is computed through the lab's own task_lab.run on the assembled lab data.

    bash run_in_image.sh make_task.py
"""
import copy, json, sys
from pathlib import Path
import numpy as np

TASK = Path('/repo/benchmarks/tasks/zrsix-butterfly-amr')
LAB = TASK/'environment/lab'
sys.path.insert(0, str(LAB))
import task_lab as TL  # noqa: E402

HERE = Path(__file__).resolve().parent
T = 30.0
INPUT_EF = {'ZrSiS': 8.9117, 'ZrSiSe': 8.9117, 'ZrSiTe': 8.8965}      # E_FERMI of the archived input.dat files
CELL = {c: TL.model(c).A for c in TL.COMPOUNDS}
THRESH = {'mu_abs_eV': 0.003, 'density_rel': 0.06, 'angle_abs_deg': 3.0, 'exponent_abs': 0.3,
          'mr_b_over_c_rel': 0.15, 'zrsite_mr_c_over_b_rel': 0.25, 'mr_ratio_se_over_te_rel': 0.2}


def count(c, mesh, mus):
    return TL.run({'kind': 'electron_count', 'compound': c, 'mesh': mesh, 'mu_eV': list(map(float, mus)),
                   'temperature_K': T}, None)


def mu_at(c, mesh, nfill):
    lo, hi = TL.SPEC['count_range'][c]
    for _ in range(3):
        g = np.linspace(lo, hi, 200)
        n = np.array(count(c, mesh, g)['electrons_per_cell'])
        i = int(np.searchsorted(n, nfill))
        lo, hi = g[max(i - 1, 0)], g[min(i, 199)]
    g = np.linspace(lo, hi, 200)
    r = count(c, mesh, g)
    n = np.array(r['electrons_per_cell']); i = int(np.searchsorted(n, nfill))
    w = (nfill - n[i - 1]) / (n[i] - n[i - 1])
    occ = np.array(r['band_occupation']); occ = occ[i - 1] + w * (occ[i] - occ[i - 1])
    return float(g[i - 1] + w * (g[i] - g[i - 1])), occ


def density(c, occ, spin=2):
    A = CELL[c]
    return float(spin * occ[12:].sum() / (abs(np.linalg.det(A)) * 1e-24))


def mr(c, plane, ang, mesh, mu, bands=None):
    bands = bands or TL.SPEC['transport_bands'][c]
    r = TL.run({'kind': 'magnetotransport', 'compound': c, 'field': {'plane': plane, 'angle_deg': ang}, 'mesh': mesh,
                'mu_eV': [float(mu)], 'temperature_K': T, 'bands': bands}, None)
    rho = np.array(r['results'][0]['rho_tau_Ohm_m_ps'])
    X = list(r['Btau_Tps'])
    return X, rho[:, 0, 0] / rho[0, 0, 0] - 1, r['results'][0]['open_orbit_fraction_by_band'], rho


def peak(th, y, how='parabola'):
    th, y = np.asarray(th, float), np.asarray(y, float)
    i = int(np.argmax(y)); i = min(max(i, 1), len(y) - 2)
    if how == 'parabola':
        p = np.polyfit(th[i - 1:i + 2], y[i - 1:i + 2], 2)
        return float(-p[1] / (2 * p[0]))
    # cubic polynomial through 5 points around the maximum, maximised on a fine grid
    j = slice(max(i - 2, 0), min(i + 3, len(y)))
    p = np.polyfit(th[j], y[j], 3)
    f = np.linspace(th[j][0], th[j][-1], 2001)
    return float(f[np.argmax(np.polyval(p, f))])


def fields(mu_se, mu_te, mesh=160, angles=None, how='parabola', se_te_btau=(5.0, 2.0)):
    angles = angles or [float(a) for a in TL.SPEC['field_menu']['ZrSiSe']['bc']]
    X, m0, _, _ = mr('ZrSiSe', 'bc', 0.0, mesh, mu_se)
    _, m90, _, _ = mr('ZrSiSe', 'bc', 90.0, mesh, mu_se)
    i5, i10, i20 = X.index(5.0), X.index(10.0), X.index(20.0)
    scan = {a: mr('ZrSiSe', 'bc', a, mesh, mu_se)[1][i5] for a in angles}
    _, t0, _, _ = mr('ZrSiTe', 'bc', 0.0, mesh, mu_te)
    _, t90, _, _ = mr('ZrSiTe', 'bc', 90.0, mesh, mu_te)
    xs, xt = X.index(se_te_btau[0]), X.index(se_te_btau[1])
    return {'butterfly_peak_angle_deg': peak(list(scan), list(scan.values()), how),
            'mr_b_over_c': float(m90[i5] / m0[i5]),
            'mr_exponent_b': float(np.log(m90[i20] / m90[i10]) / np.log(2)),
            'zrsite_mr_c_over_b': float(t0[i5] / t90[i5]),
            'mr_ratio_se_over_te': float(m0[xs] / t0[xt])}, scan


def catalog(spec):
    return {
        'material': 'ZrSiX (X = S, Se, Te), P4/nmm; GGA Wannier models without spin-orbit coupling released by Zhang and '
                    'Yazyev (2025), Materials Cloud record tyan3-q3j67',
        'hamiltonians': spec['files'],
        'orbitals': 'per cell: Zr1 d(5), Zr2 d(5), Si1 s p(4), Si2 s p(4), X1 p(3), X2 p(3); no spin index (24 bands)',
        'cells_A': {c: {'a': round(float(CELL[c][0, 0]), 7), 'c': round(float(CELL[c][2, 2]), 7),
                        'source': 'POSCAR / wannier90.win of the compound'} for c in TL.COMPOUNDS},
        'budget': json.loads((LAB/'public.json').read_text())['solver_budget'],
        'requests': {
            'bands': {'cost': 1, 'format': {'kind': 'bands', 'compound': 'ZrSiSe', 'k_frac': [[0, 0, 0]], 'velocities': False},
                      'limits': 'k_frac: 1-500 points; velocities true returns dE/dk (eV Angstrom) per band'},
            'fermi_slice': {'cost': 1, 'format': {'kind': 'fermi_slice', 'compound': 'ZrSiSe', 'origin': [-0.5, -0.5, 0.0],
                                                  'u': [1, 0, 0], 'v': [0, 1, 0], 'n': [61, 61]},
                            'limits': 'k = origin + i/(nu-1) u + j/(nv-1) v (fractional); 2 <= nu, nv <= 121; returns bands 10-15'},
            'electron_count': {'cost': 1, 'format': {'kind': 'electron_count', 'compound': 'ZrSiSe', 'mesh': 80,
                                                     'mu_eV': [8.9], 'temperature_K': 30},
                               'limits': 'mesh in ' + str(spec['count_meshes']) + ' (Gamma-centred N x N x N/2); 1-200 '
                                         'chemical potentials within mu_range_eV; temperature 0-300 K'},
            'orbit': {'cost': 1, 'format': {'kind': 'orbit', 'compound': 'ZrSiSe', 'band': 12, 'k_frac': [[0.2, 0.2, 0.0]],
                                            'field_direction': [0, 0, 1], 'max_length_invA': 20},
                      'limits': '1-8 start points; any field direction (Cartesian); arc length 0.1-40 1/Angstrom'},
            'magnetotransport': {'cost': {'40': 1, '80': 1, '160': 3},
                                 'format': {'kind': 'magnetotransport', 'compound': 'ZrSiSe',
                                            'field': {'plane': 'bc', 'angle_deg': 45}, 'mesh': 80, 'mu_eV': [8.91],
                                            'temperature_K': 30, 'bands': [12, 13]},
                                 'limits': '1-4 chemical potentials within mu_range_eV; temperature 10-40 K; bands from '
                                           'the list below; field from field_menu'}},
        'electron_count': {'mu_range_eV': spec['count_range'], 'meshes': spec['count_meshes']},
        'magnetotransport': {
            'mu_range_eV': spec['transport_mu_range'], 'bands': spec['transport_bands'], 'meshes': [40, 80, 160],
            'Btau_Tps': [float(x) for x in np.load(LAB/'data/mr_ZrSiSe.npz')['X']],
            'field_menu': spec['field_menu'],
            'field_planes': {'bc': 'B = (0, sin a, cos a): from c (a = 0) towards b (a = 90)',
                             'ac': 'B = (sin a, 0, cos a): from c (a = 0) towards a (a = 90)',
                             'ab': 'B = (cos a, sin a, 0): in the ab plane, a from the a axis'},
            'method': 'Chambers formula, one relaxation time for all bands; orbits of every start point of the '
                      'Gamma-centred N x N x N/2 mesh traced on the start point\'s constant-energy surface (the 40 and 80 '
                      'meshes are subsets of the 160 mesh); results scale as rho = rho_tau / tau.'},
    }


def main():
    spec = TL.SPEC
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog(spec), indent=1) + '\n')
    ref_mu, occ = {}, {}
    for c in TL.COMPOUNDS:
        ref_mu[c], occ[c] = mu_at(c, 160, 12)
    ref = {'neutral_mu_eV': ref_mu, 'carrier_density_cm3': {c: density(c, occ[c]) for c in ('ZrSiSe', 'ZrSiTe')}}
    f160, scan160 = fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'])
    ref.update(f160)
    # robustness: mesh-120 potentials, cubic peak finder, the 5 deg scan restricted to the oracle's angles
    mu120 = {c: mu_at(c, 120, 12)[0] for c in TL.COMPOUNDS}
    occ120 = {c: mu_at(c, 120, 12)[1] for c in ('ZrSiSe', 'ZrSiTe')}
    rob_mu120 = dict(neutral_mu_eV=mu120, carrier_density_cm3={c: density(c, occ120[c]) for c in occ120},
                     **fields(mu120['ZrSiSe'], mu120['ZrSiTe'])[0])
    rob_cubic = fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'], how='cubic')[0]
    rob_oracle_angles = fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'], angles=[0.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 90.0])[0]
    # ablations
    ab = {}
    ab['input_e_fermi'] = dict(neutral_mu_eV=dict(INPUT_EF),
                               **fields(INPUT_EF['ZrSiSe'], INPUT_EF['ZrSiTe'])[0])
    ab['zrsis_copied_input'] = {'neutral_mu_eV': dict(ref_mu, ZrSiS=INPUT_EF['ZrSiS'])}
    m13 = {c: mu_at(c, 160, 13)[0] for c in TL.COMPOUNDS}
    ab['thirteen_filled_bands'] = {'neutral_mu_eV': dict(ref_mu, ZrSiTe=m13['ZrSiTe']), 'mu13': m13}
    ab['carrier_density_one_spin'] = {'carrier_density_cm3': {c: density(c, occ[c], spin=1) for c in ('ZrSiSe', 'ZrSiTe')}}
    ab['mesh_80'] = fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'], mesh=80)[0]
    ab['mesh_40'] = fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'], mesh=40)[0]
    end = max((0.0, 90.0), key=lambda a: scan160[a])
    ab['closed_orbit_mass_model'] = {'butterfly_peak_angle_deg': end}
    ab['samples_at_equal_btau'] = {'mr_ratio_se_over_te': fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'], se_te_btau=(5.0, 5.0))[0]['mr_ratio_se_over_te']}
    ab['samples_at_equal_btau_swapped'] = {'mr_ratio_se_over_te': fields(ref_mu['ZrSiSe'], ref_mu['ZrSiTe'], se_te_btau=(2.0, 5.0))[0]['mr_ratio_se_over_te']}
    hidden = {'targets': ref, 'thresholds': THRESH,
              'scan_160_Btau5': {str(k): float(v) for k, v in scan160.items()}}
    (LAB/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(ref['neutral_mu_eV'].values()) + list(ref['carrier_density_cm3'].values()) + \
        [ref[k] for k in f160]
    (LAB/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    rob = {'mesh120_potentials': rob_mu120, 'cubic_peak': rob_cubic, 'oracle_angles': rob_oracle_angles}
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': ab, 'robustness': rob}, indent=1) + '\n')
    (HERE/'calibration.json').write_text(json.dumps({'reference': ref, 'robustness': rob, 'ablations': ab,
                                                      'thresholds': THRESH}, indent=1) + '\n')
    print(json.dumps({'reference': ref, 'robustness': rob, 'ablations': ab}, indent=1))


if __name__ == '__main__':
    main()
