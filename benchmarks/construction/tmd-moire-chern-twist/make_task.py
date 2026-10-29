"""Write the staged, lab, private and ablation data of tmd-moire-chern-twist.

    bash benchmarks/construction/tmd-moire-chern-twist/run_in_image.sh make_task.py

Reads the extracted archive (fetch.py, mounted at /data), c3_labels.json (transcribed Fig. 2 labels) and the
reference calculations in reference/ (reference.py). Writes into benchmarks/tasks/tmd-moire-chern-twist.
"""
import cmath
import hashlib
import json
import math
import shutil
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/tmd-moire-chern-twist'
SRC = Path('/data/x/source-data')
AGENT, LAB = TASK/'environment/agent', TASK/'environment/lab'
BUDGET = 27
BAND_RUNS = {'tWSe2': {'3.15': 'fig2c/3.15', '1.70': 'fig2c/1.70', '1.54': 'fig2c/1.54', '1.25': 'fig2c/1.25'},
             'tMoTe2': {'3.89': 'fig2d/3.89', '3.15': 'fig2d/3.15', '2.14': 'fig2d/2.14', '1.25': 'fig2d/1.25'}}
GRIDS = {'3.15': 320, '1.70': 300, '1.54': 320, '1.25': 360}
C3VAL = {'xi': cmath.exp(1j * math.pi / 3), 'xis': cmath.exp(-1j * math.pi / 3), 'm1': -1 + 0j}
CONJ = {'xi': 'xis', 'xis': 'xi', 'm1': 'm1'}
SET_FOLDER = {'P': 'fig1b', 'Q': 'fig1a'}
SET_PARAMS = {'P': {'effective_mass_me': 0.60, 'monolayer_lattice_constant_A': 3.52},
              'Q': {'effective_mass_me': 0.62, 'monolayer_lattice_constant_A': 3.472}}
THRESHOLDS = {'hartree_xm_abs_meV': 1.0, 'first_star_abs_meV': 0.25, 'first_star_rel': 0.05}


def bands(rel):
    raw = (SRC/'fig2'/rel/'Es-band.dat').read_text().strip().split('\n\n')
    B = [np.array([[float(x) for x in l.split()] for l in b.strip().split('\n')]) for b in raw]
    return B[0][:, 0], np.array([b[:, 1] for b in B])[::-1] * 1000.0      # descending energy, meV


def chern_from_labels(g, k, kp):
    z = -C3VAL[g] * C3VAL[k] * C3VAL[kp]
    c = round(cmath.phase(z) / (2 * math.pi / 3))
    return int(c)


def complex_pair(z):
    return [round(z.real, 6) + 0.0, round(z.imag, 6) + 0.0]


def flip(*key):
    return int(hashlib.sha256('|'.join(key).encode()).hexdigest(), 16) % 2 == 1


def c3_records(labels):
    out = {}
    for mat, runs in BAND_RUNS.items():
        out[mat] = {}
        for angle, rel in runs.items():
            x, E = bands(rel)
            out[mat][angle] = {}
            for point, col, lab_key in (('gamma', 0, 'gamma'), ('kappa', 40, 'kappa'), ('kappa_prime', None, 'kappa')):
                states = []
                for n in range(3):
                    lab_up = labels[mat][angle][lab_key][n]
                    pair = [{'spin': 'up', 'C3z_eigenvalue': complex_pair(C3VAL[lab_up])},
                            {'spin': 'down', 'C3z_eigenvalue': complex_pair(C3VAL[CONJ[lab_up]])}]
                    if flip(mat, angle, point, str(n)):
                        pair = pair[::-1]
                    for s in pair:
                        if col is not None:
                            e = 0.5 * (E[2 * n, col] + E[2 * n + 1, col])
                            s = {'energy_meV': round(float(e), 2), **s}
                        states.append(s)
                assert col is None or all(abs(E[2 * n, col] - E[2 * n + 1, col]) < 0.2 for n in range(3))
                rec = {'states': [{'state': i + 1, **s} for i, s in enumerate(states)],
                       'note': 'the six highest moire valence states at this point, highest first (pairs of states '
                               'that are degenerate within 0.2 meV share the pair-averaged energy)'}
                if col is None:
                    rec['note'] += '; energies are not archived at kappa_prime (it is not on the band path)'
                out[mat][angle][point] = rec
    return out


def main():
    labels = json.loads((HERE/'c3_labels.json').read_text())['labels']
    cont = json.loads((HERE/'reference/continuum.json').read_text())
    hart = json.loads((HERE/'reference/hartree.json').read_text())
    for d in (AGENT/'archive', LAB/'data'):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    # ---- targets ---------------------------------------------------------------------------------------------
    dft_chern = {}
    for mat, runs in labels.items():
        for angle, lab in runs.items():
            dft_chern[f'{mat}_{angle}'] = [chern_from_labels(lab['gamma'][n], lab['kappa'][n], lab['kappa'][n]) for n in range(2)]
    rows = cont['rows']
    def row(label, theta, ns, nk, variant=None):
        return next(r for r in rows if r['set'] == label and r['theta'] == theta and r['nshell'] == ns and r['nk'] == nk
                    and r['variant'] == variant)
    cont_chern = {s: row(s, 3.89, 5, 18)['chern'] for s in 'PQ'}
    for s in 'PQ':
        assert row(s, 3.89, 7, 24)['chern'] == cont_chern[s] == row(s, 3.89, 4, 12)['chern']
    p214 = row('P', 2.14, 7, 24)['chern']
    assert row('P', 2.14, 5, 18)['chern'] == p214
    match = min(cont['rms_vs_dft_3p89_meV'], key=cont['rms_vs_dft_3p89_meV'].get)
    assert cont_chern[match] == dft_chern['tMoTe2_3.89']
    xm = {a: round(hart[a]['xm_fourier'], 4) for a in GRIDS}
    fs = {a: round(hart[a]['first_star_xm'], 4) for a in ('1.70', '1.25')}
    targets = {'dft_chern_K_up': dft_chern, 'continuum_chern_3p89': cont_chern, 'continuum_set_matching_dft': match,
               'continuum_P_chern_2p14': p214, 'hartree_xm_meV': xm, 'hartree_first_star_xm_meV': fs}
    (LAB/'private').mkdir(exist_ok=True)
    (LAB/'private/hidden.json').write_text(json.dumps({'targets': targets, 'thresholds': THRESHOLDS,
                                                       'reference_details': {'hartree': hart}}, indent=1) + '\n')
    (LAB/'private/hidden_values.json').write_text(json.dumps(list(xm.values()) + list(fs.values())) + '\n')

    # ---- lab data --------------------------------------------------------------------------------------------
    archive = {'bands': {}, 'c3': c3_records(labels), 'linecut': {}, 'fourier': {}}
    for mat, runs in BAND_RUNS.items():
        archive['bands'][mat] = {}
        for angle, rel in runs.items():
            x, E = bands(rel)
            archive['bands'][mat][angle] = {
                'k_distance_invA': [round(float(v), 6) for v in x],
                'path': 'gamma (index 0) - m (index 20) - kappa (index 40) - gamma (index 60) of the moire Brillouin zone',
                'energies_meV': [[round(float(v), 3) for v in e] for e in E],
                'note': f'archived SIESTA band energies (Es-band.dat), {len(E)} states, highest first; energies relative '
                        'to the valence-band maximum as archived'}
    for angle, n in GRIDS.items():
        bot = 'pos_2d_bot_sm.dat' if angle == '3.15' else 'pot_2d_bot_sm.dat'
        top = np.loadtxt(SRC/'fig3/fig3a-d'/angle/'pot_2d_top_sm.dat').reshape(n, n)
        b = np.loadtxt(SRC/'fig3/fig3a-d'/angle/bot).reshape(n, n)
        np.savez_compressed(LAB/'data'/f'hartree_{angle}.npz', top=top, bot=b, n=n)
        cut = np.loadtxt(SRC/'fig3/fig3e'/angle/'pot_diag-top.dat')
        archive['linecut'][angle] = {'columns': ['distance along the long diagonal a1 + a2 from the origin, in units of |a1|',
                                                 'archived value (eV)'],
                                     'data': [[round(float(u), 6), float('%.10g' % v)] for u, v in cut],
                                     'note': 'archived file fig3e/<angle>/pot_diag-top.dat, written by plot-pot.py'}
        g = np.loadtxt(SRC/'fig3/fig3f'/angle/'pot_top_bot_G.dat')
        keep = g[g[:, 0] <= 7.0001]
        archive['fourier'][angle] = {'columns': ['|g| / g1', 'archived Fourier amplitude (eV)'],
                                     'data': [[round(float(u), 6), float('%.10g' % v)] for u, v in keep],
                                     'note': 'archived file fig3f/<angle>/pot_top_bot_G.dat (rows with |g|/g1 <= 7); '
                                             'read by plot-fft-pot.py'}
        shutil.copy(SRC/'fig3/fig3a-d'/angle/'POSCAR', LAB/'data'/f'POSCAR_{angle}')
    archive['hartree_map_note'] = {
        'grid': 'v[i][j] is the value at fractional coordinates (i/(n-1), j/(n-1)) of the lattice vectors a1, a2 of the '
                'relaxed moire cell of this angle (CELL in the structure request); row i = n-1 repeats row 0',
        'quantity': 'smoothed (coarse-grained) DFT Hartree potential in the vacuum above the top-layer surface and below '
                    'the bottom-layer surface (eV), archived files pot_2d_top_sm.dat / pot_2d_bot_sm.dat; top = larger z'}
    (LAB/'data/archive.json').write_text(json.dumps(archive) + '\n')

    # ---- staged files ----------------------------------------------------------------------------------------
    A = AGENT/'archive'
    for s, folder in SET_FOLDER.items():
        (A/f'fig1/set_{s}').mkdir(parents=True)
        for f in ('Delta-top.dat', 'Delta-bot.dat', 'Delta-T.dat', 'plot-sky.py'):
            shutil.copy(SRC/'fig1'/folder/f, A/f'fig1/set_{s}'/f)
    for mat, angle, rel, fdf in (('tWSe2', '3.15', 'fig2c/3.15', 'wse2.fdf'), ('tMoTe2', '3.89', 'fig2d/3.89', 'mote2.fdf')):
        dst = A/f'fig2/{mat}_{angle}'
        dst.mkdir(parents=True)
        shutil.copy(SRC/'fig2'/rel/'Es-band.dat', dst/'Es-band.dat')
        shutil.copy(SRC/'fig2'/rel/fdf, dst/fdf)
    for angle in ('3.15', '1.25'):
        dst = A/f'fig2/tWSe2_relaxation_{angle}'
        dst.mkdir(parents=True)
        for f in ('CONTCAR', 'ILD-int.dat', 'displace-int.dat'):
            shutil.copy(SRC/'fig2/fig2a-b'/angle/f, dst/f)
    dst = A/'fig3/tWSe2_3.15'
    dst.mkdir(parents=True)
    for f in ('pot_2d_top_sm.dat', 'pos_2d_bot_sm.dat', 'plot-pot.input', 'plot-pot.py'):
        shutil.copy(SRC/'fig3/fig3a-d/3.15'/f, dst/f)
    shutil.copy(SRC/'fig3/fig3e/3.15/pot_diag-top.dat', dst/'pot_diag-top.dat')
    shutil.copy(SRC/'fig3/fig3e/plot-diag-pot.py', A/'fig3/plot-diag-pot.py')
    shutil.copy(SRC/'fig3/fig3f/plot-fft-pot.py', A/'fig3/plot-fft-pot.py')

    catalog = {
        'study': 'X.-W. Zhang, C. Wang, X. Liu, Y. Fan, T. Cao, D. Xiao, Polarization-driven band topology evolution in '
                 'twisted MoTe2 and WSe2, Nat. Commun. 15, 4223 (2024); data: Zenodo 10888769 (CC BY 4.0)',
        'systems': {'tWSe2': {'angles_deg': list(BAND_RUNS['tWSe2'])}, 'tMoTe2': {'angles_deg': list(BAND_RUNS['tMoTe2'])}},
        'staged_files': {
            'archive/fig1/set_P, archive/fig1/set_Q': 'the two continuum-model parameter maps of the paper\'s Fig. 1 '
                '(one from the local-stacking approximation, one fitted to a 3.89 deg tMoTe2 DFT calculation; the labels '
                'P and Q are arbitrary): Delta-bot.dat, Delta-top.dat (columns Re, Im, |.| in meV, 100 x 100 grid) and '
                'Delta-T.dat (Re, Im, |.|, 10 x 10 grid), with the authors\' plotting script',
            'archive/fig2/tWSe2_3.15, archive/fig2/tMoTe2_3.89': 'archived SIESTA moire valence bands (Es-band.dat: blocks of '
                '(k distance in 1/A, energy in eV) per state along gamma-m-kappa-gamma, 61 points) and the SIESTA input',
            'archive/fig2/tWSe2_relaxation_*': 'relaxed tWSe2 structures (CONTCAR) and interpolated interlayer distance / '
                'in-plane displacement maps',
            'archive/fig3/tWSe2_3.15': 'surface Hartree potentials of 3.15 deg tWSe2 (top: pot_2d_top_sm.dat; bottom: '
                'pos_2d_bot_sm.dat, archived name), the line cut and the authors\' scripts (the Fourier-amplitude files are sold by hartree_fourier)'},
        'continuum_parameter_sets': {s: {**SET_PARAMS[s], 'files': f'archive/fig1/set_{s}'} for s in 'PQ'},
        'requests': {
            'dft_bands': {'format': {'kind': 'dft_bands', 'material': 'tWSe2 | tMoTe2', 'angle': '<angle string>'},
                          'cost': 1, 'returns': 'archived SIESTA moire valence bands along gamma-m-kappa-gamma'},
            'c3': {'format': {'kind': 'c3', 'material': 'tWSe2 | tMoTe2', 'angle': '<angle string>',
                              'point': 'gamma | kappa | kappa_prime'}, 'cost': 1,
                   'returns': 'spin and C3z eigenvalue (complex, [Re, Im]) of the six highest moire valence states at '
                              'the point, as determined by the authors'},
            'hartree_map': {'format': {'kind': 'hartree_map', 'angle': '3.15 | 1.70 | 1.54 | 1.25'}, 'cost': 3,
                            'returns': 'tWSe2 top- and bottom-surface Hartree potential maps (n x n; a few MB of JSON)'},
            'hartree_linecut': {'format': {'kind': 'hartree_linecut', 'angle': '...'}, 'cost': 1,
                                'returns': 'the archived line cut pot_diag-top.dat of that angle'},
            'hartree_fourier': {'format': {'kind': 'hartree_fourier', 'angle': '...'}, 'cost': 1,
                                'returns': 'the archived Fourier-amplitude file pot_top_bot_G.dat of that angle (|g|/g1 <= 7)'},
            'structure': {'format': {'kind': 'structure', 'angle': '...'}, 'cost': 1,
                          'returns': 'relaxed tWSe2 structure (POSCAR) of that angle'}},
        'budget': BUDGET}
    (AGENT/'catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    public = {'task_id': 'tmd-moire-chern-twist', 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
              'requests': {k: {**v['format'], 'cost': v['cost']} for k, v in catalog['requests'].items()},
              'catalog': '/app/catalog.json'}
    (LAB/'public.json').write_text(json.dumps(public, indent=1) + '\n')
    (TASK/'tests/public.json').write_text(json.dumps({'task_id': 'tmd-moire-chern-twist', 'solver_budget': BUDGET,
                                                      'min_report_chars': 400}, indent=1) + '\n')

    # ---- ablation values -------------------------------------------------------------------------------------
    def chern_all(rule):
        out = {}
        for mat, runs in labels.items():
            for angle, lab in runs.items():
                out[f'{mat}_{angle}'] = rule(mat, angle, lab)
        return out
    def phase_c(z):
        return int(round(cmath.phase(z) / (2 * math.pi / 3)))
    abl = {
        'dft_spin_down_sector': chern_all(lambda m, a, l: [-c for c in dft_chern[f'{m}_{a}']]),
        'dft_no_minus_sign': chern_all(lambda m, a, l: [phase_c(C3VAL[l['gamma'][n]] * C3VAL[l['kappa'][n]] ** 2)
                                                         for n in range(2)]),
        'dft_kappa_prime_conjugate': chern_all(lambda m, a, l: [phase_c(-C3VAL[l['gamma'][n]]) for n in range(2)]),
    }
    # bands indexed by state order (states 1, 2 of the flat list = "bands 1, 2"), spin labels ignored
    c3 = archive['c3']
    def state_order(m, a, l):
        out = []
        for n in range(2):
            zs = {}
            for point in ('gamma', 'kappa', 'kappa_prime'):
                v = c3[m][a][point]['states'][n]['C3z_eigenvalue']
                zs[point] = complex(*v)
            out.append(phase_c(-zs['gamma'] * zs['kappa'] * zs['kappa_prime']))
        return out
    abl['dft_state_index_as_band'] = chern_all(state_order)
    abl['dft_from_fixed_continuum_2p14'] = {**dft_chern, 'tMoTe2_2.14': p214}
    for v in ('layers_swapped', 'plot_script_basis', 'left_handed_basis'):
        abl[f'continuum_{v}'] = {s: row(s, 3.89, 5, 18, v)['chern'] for s in 'PQ'}
    abl['continuum_set_literature_prior'] = 'Q' if match == 'P' else 'P'
    abl['continuum_2p14_equals_dft'] = dft_chern['tMoTe2_2.14']
    abl['hartree_xm_mx_swapped'] = {a: round(hart[a]['mx_fourier'], 4) for a in GRIDS}
    abl['first_star_xm_mx_swapped'] = {a: round(hart[a]['first_star_mx'], 4) for a in fs}
    abl['hartree_xm_first_star_only'] = {a: round(hart[a]['first_star_xm'], 4) for a in GRIDS}
    abl['first_star_equals_full_xm'] = {a: xm[a] for a in fs}
    abl['first_star_from_magnitudes_3p15_sign'] = {a: round(hart[a]['first_star_xm_from_magnitudes_with_3p15_sign'], 4) for a in fs}
    abl['hartree_xm_no_sign_flip'] = {a: -abs(xm[a]) for a in GRIDS}
    rob = {'hartree_xm_nearest_grid': {a: round(hart[a]['xm_nearest'], 4) for a in GRIDS},
           'hartree_xm_bilinear': {a: round(hart[a]['xm_bilinear'], 4) for a in GRIDS},
           'hartree_xm_archived_linecut': {a: round(hart[a]['xm_linecut'], 4) for a in GRIDS},
           'first_star_full_grid_fft': {a: round(hart[a]['first_star_xm_fullgrid'], 4) for a in fs}}
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')

    calib = {'targets': targets, 'ablations': abl, 'robustness': rob, 'continuum_rms_vs_dft_3p89_meV': cont['rms_vs_dft_3p89_meV']}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    print(json.dumps(targets, indent=1))
    for k, v in abl.items():
        print(k, v)


if __name__ == '__main__':
    main()
