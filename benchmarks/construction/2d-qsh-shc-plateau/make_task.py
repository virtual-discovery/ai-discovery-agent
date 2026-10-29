"""Write the staged, lab, private and ablation files of 2d-qsh-shc-plateau (runs in the runtime image).

    bash run_in_image.sh make_task.py

Inputs: the verified archive members under /data/x (fetch.py), the lab models in environment/lab/data/*.npz
(build_models.py) and the reference calculations in reference/*.json (reference.py). Asserts the archive facts the
task relies on and prints the calibration numbers.
"""
import json
import math
import re
import shutil
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/2d-qsh-shc-plateau'
SRC = Path('/data/x/Data')
MATERIALS = {'Bi2_P-3m1': 'Bi2_P-3m1', 'Hf2Br2_P-3m1': 'Hf2Br2_P-3m1', 'CuLi2As_P-6m2': 'CuLi2As_P-6m2',
             'Ti2I2N2_Pmm2': 'Ti2I2N2_Pmm2', 'Li2Tl2_P3m1': 'Li2Tl2_P3m1', 'Bi2In2_P4_nmm': 'Bi2In2_P4|nmm',
             'Hg4O2_P1': 'Hg4O2_P1', 'KSnAs_P3m1': 'KSnAs_P3m1'}
STAGED = ['aiida.in', 'aiida.out', 'aiida.win', 'cut_aiida.wout']
TOP = {'Bi2': 'Bi2_P-3m1', 'Hf2Br2': 'Hf2Br2_P-3m1', 'CuLi2As': 'CuLi2As_P-6m2', 'Ti2N2I2': 'Ti2I2N2_Pmm2'}
SPIN_DIRECTION = ['Bi2_P-3m1', 'Hg4O2_P1']
BUDGET = 40
THRESHOLDS = {'plateau_rel': 0.03, 'plateau_abs_floor': 0.02, 'theta_deg': 2.5, 'phi_deg': 5.0}
E2H = 3.874045864e-5
BOHR = 0.529177210903


def cell_of(folder):
    t = (folder/'aiida.in').read_text()
    return np.array([[float(x) for x in l.split()] for l in t.split('CELL_PARAMETERS')[1].split('\n')[1:4]])


def angles(v):
    x, y, z = v
    return {'theta_deg': math.degrees(math.atan2(math.hypot(x, y), z)), 'phi_deg': math.degrees(math.atan2(y, x))}


def plateau_ok(value, ref):
    return abs(value - ref) <= max(THRESHOLDS['plateau_rel'] * abs(ref), THRESHOLDS['plateau_abs_floor'])


def read_top(name, comp):
    folder = SRC/'Top_materials'/name/f'SHC-{comp}'
    win = (folder/'aiida.win').read_text()
    rows = np.loadtxt(folder/'aiida-shc-fermiscan.dat')
    cell = re.search(r'Begin Unit_Cell_Cart\s*\n\s*bohr\s*\n(.*?)End Unit_Cell_Cart', win, re.S | re.I).group(1)
    c_bohr = float(cell.strip().split('\n')[2].split()[2])
    return win, rows, c_bohr


def main():
    ref = {n: json.loads((HERE/'reference'/f'{n}.json').read_text()) for n in MATERIALS}
    lab, agent = TASK/'environment/lab', TASK/'environment/agent'
    # ---------------------------------------------------------------- staged archive files
    if (agent/'archive').exists():
        shutil.rmtree(agent/'archive')
    for n, folder in MATERIALS.items():
        (agent/'archive'/n).mkdir(parents=True)
        for f in STAGED:
            src = SRC/'Wannier_Hamiltonians/insulators'/folder/f
            if f == 'cut_aiida.wout':      # header and final state; the iteration logs are omitted
                lines = src.read_text().split('\n')
                head = next(i for i, l in enumerate(lines) if 'Extraction of optimally-connected subspace' in l
                            or 'Initial State' in l)
                fin = max(i for i, l in enumerate(lines) if 'Final State' in l)
                text = '\n'.join(lines[:head] + ['', ' [... disentanglement and minimisation iterations omitted ...]', ''] + lines[fin:])
                (agent/'archive'/n/'cut_aiida.wout.excerpt').write_text(text)
            elif f == 'aiida.out':          # the band energies of the final scf iteration are omitted
                lines = src.read_text().split('\n')
                end = max(i for i, l in enumerate(lines) if 'End of self-consistent calculation' in l)
                ef = next(i for i in range(end, len(lines)) if 'the Fermi energy is' in lines[i])
                nk = sum(1 for l in lines[end:ef] if l.strip().startswith('k ='))
                text = '\n'.join(lines[:end + 1] + ['', f'     [... band energies and occupations at the {nk} k-points omitted ...]', ''] + lines[ef:])
                (agent/'archive'/n/'aiida.out.excerpt').write_text(text)
            else:
                shutil.copy2(src, agent/'archive'/n/f)
    shutil.copy2(SRC/'README.txt', agent/'archive/README.txt')
    # ---------------------------------------------------------------- lab data
    meta = {'materials': {}}
    for n, folder in MATERIALS.items():
        A = cell_of(SRC/'Wannier_Hamiltonians/insulators'/folder)
        z = np.load(lab/'data'/f'{n}.npz')
        meta['materials'][n] = {'archive_folder': f'Wannier_Hamiltonians/insulators/{folder}', 'cell_A': A.tolist(),
                                'num_wann': int(z['H'].shape[1])}
    (lab/'data/meta.json').write_text(json.dumps(meta, indent=1) + '\n')
    archived = {}
    for top in TOP:
        for comp in ('xyx', 'xyy', 'xyz'):
            win, rows, _ = read_top(top, comp)
            archived[f'Top_materials/{top}/SHC-{comp}'] = {
                'files': {'aiida.win': win},
                'fermiscan_columns': ['No.', 'Fermi energy (eV)', 'SHC ((hbar/e)*S/cm)'],
                'fermiscan_rows': [[int(r[0]), float(r[1]), float(r[2])] for r in rows],
                'note': 'aiida-shc-fermiscan.dat and aiida.win of this archive folder, values exactly as in the files'}
    (lab/'data/archived_shc.json').write_text(json.dumps(archived) + '\n')
    for f in ('qsh_tb.py', 'shc_worker.py'):
        shutil.copy2(HERE/f, lab/f)
    # ---------------------------------------------------------------- archive facts the task relies on
    facts = {}
    for n, r in ref.items():
        vbm, cbm = r['band_edges']['vbm'][0], r['band_edges']['cbm'][0]
        facts[n] = {'nocc': r['nocc'], 'wannier_vbm': vbm, 'wannier_cbm': cbm, 'wannier_gap': cbm - vbm,
                    'scf_fermi': r['scf_fermi_eV'], 'win_vbm': r['win_VBM_eV'], 'win_cbm': r['win_CBM_eV'],
                    'scf_fermi_in_valence_band': r['scf_fermi_eV'] < vbm,
                    'win_mid_gap_in_valence_band': r['win_mid_gap_eV'] < vbm}
        facts[n]['win_fermi_energy_is_vbm'] = abs(r['win_fermi_energy_eV'] - r['win_VBM_eV']) < 1e-9
    assert sum(f['win_fermi_energy_is_vbm'] for f in facts.values()) == 7      # all but Hg4O2 (set to the DFT mid-gap)
    assert facts['Bi2In2_P4_nmm']['win_mid_gap_in_valence_band']
    assert sum(f['scf_fermi_in_valence_band'] for f in facts.values()) >= 5
    # ---------------------------------------------------------------- references
    sheet = {n: {lab_: np.array(r['shc'][lab_]['sheet_xyz']) for lab_ in r['shc']} for n, r in ref.items()}
    MU_MID, MU_SCF, MU_WINMID, MU_WINEF = 0, 1, 2, 3
    z2 = {}
    for n, r in ref.items():
        best = [r['z2'][f'201x300_axis{a}'][0] for a in (0, 1)]
        assert best[0] == best[1], n
        z2[n] = best[0]
    plateau = {n: float(sheet[n]['reference'][2, MU_MID]) for n in MATERIALS}
    spin_dir = {n: angles(sheet[n]['reference'][:, MU_MID]) for n in SPIN_DIRECTION}
    hidden = {'targets': {'z2': z2, 'plateau_sigma_z': plateau, 'spin_direction': spin_dir},
              'thresholds': THRESHOLDS,
              'reference_definition': 'archived Wannier models, wannier_centres position convention, eta = 1e-4 eV, '
                                      '180 x 180 mesh + 8 x 8 refinement of points with |Omega| > 10 A^2, '
                                      'mu at the mid-gap of the 300 x 300 band edges; Z2 from 201 x 300 Wilson loops',
              'facts': facts}
    (lab/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(plateau.values()) + [v for d in spin_dir.values() for v in d.values()]
    (lab/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    needles = {t for v in values for t in {f'{v:.4g}', f'{v:.5g}', f'{v:.6g}', repr(v)} if len(t.replace('-', '').replace('.', '')) >= 4}
    hits = [f'{p.relative_to(TASK)}: {n}' for p in agent.rglob('*') if p.is_file() for n in needles if n in p.read_text(errors='ignore')]
    print('leak scan of staged files:', hits or 'clean')
    # ---------------------------------------------------------------- ablations and robustness
    def pl(label, mu=MU_MID, comp=2):
        return {n: float(sheet[n][label][comp, mu]) for n in MATERIALS}

    def sd(label, mu=MU_MID):
        return {n: angles(sheet[n][label][:, mu]) for n in SPIN_DIRECTION}

    top_plateau, top_vec = {}, {}
    for top, n in TOP.items():
        vec = []
        for comp in ('xyx', 'xyy', 'xyz'):
            win, rows, c_bohr = read_top(top, comp)
            vbm = float(re.search(r'# VBM = (\S+)', win).group(1))
            cbm = float(re.search(r'# CBM = (\S+)', win).group(1))
            i = int(np.argmin(abs(rows[:, 1] - 0.5 * (vbm + cbm))))
            vec.append(rows[i, 2] * c_bohr * BOHR * 1e-8 / E2H)
        top_vec[n] = vec
        top_plateau[n] = vec[2]
    thick = {}
    for n, folder in MATERIALS.items():
        t = (SRC/'Wannier_Hamiltonians/insulators'/folder/'aiida.in').read_text()
        zs = [float(l.split()[3]) for l in t.split('ATOMIC_POSITIONS')[1].split('K_POINTS')[0].strip().split('\n')[1:]]
        thick[n] = max(zs) - min(zs)
    A = {
        'plateau_at_scf_fermi': {'plateau_sigma_z': pl('reference', MU_SCF), 'spin_direction': sd('reference', MU_SCF)},
        'plateau_at_win_fermi_energy': {'plateau_sigma_z': pl('reference', MU_WINEF)},
        'plateau_at_dft_mid_gap': {'plateau_sigma_z': pl('reference', MU_WINMID)},
        'mesh_120_unrefined': {'plateau_sigma_z': pl('mesh120_eta1meV')},
        'mesh_60_unrefined': {'plateau_sigma_z': pl('mesh60_eta1meV'), 'spin_direction': sd('mesh60_eta1meV')},
        'lattice_position_convention': {'plateau_sigma_z': pl('lattice'), 'spin_direction': sd('lattice')},
        'relaxed_archive_values': {'plateau_sigma_z': {**plateau, **top_plateau},
                                   'spin_direction': {**spin_dir, 'Bi2_P-3m1': angles(top_vec['Bi2_P-3m1'])}},
        'paper_table_values': {'plateau_sigma_z': {**plateau, 'Bi2_P-3m1': 0.5, 'Hf2Br2_P-3m1': -0.9,
                                                   'CuLi2As_P-6m2': 1.0, 'Ti2I2N2_Pmm2': -1.0}},
        'z2_from_plateau_magnitude': {'z2': {n: int(abs(v) >= 0.5) for n, v in plateau.items()}},
        'z2_from_coarse_wcc': {'z2': {n: ref[n]['z2']['51x80_axis1'][0] for n in MATERIALS}},
        'atomic_thickness_normalisation': {'plateau_sigma_z': {n: plateau[n] * thick[n] / ref[n]['c_A'] for n in MATERIALS}},
        'e_over_4pi_units': {'plateau_sigma_z': {n: 2 * v for n, v in plateau.items()}},
        'sigma_z_only_spin_direction': {'spin_direction': {n: {'theta_deg': 0.0 if plateau[n] > 0 else 180.0, 'phi_deg': 0.0}
                                                           for n in SPIN_DIRECTION}},
    }
    R = {
        'authors_broadening_5meV': {'plateau_sigma_z': pl('eta5meV'), 'spin_direction': sd('eta5meV')},
        'mesh_240_unrefined_eta_1meV': {'plateau_sigma_z': pl('mesh240_eta1meV'), 'spin_direction': sd('mesh240_eta1meV')},
        'z2_wcc_201x200_axis0': {'z2': {n: ref[n]['z2']['201x200_axis0'][0] for n in MATERIALS}},
    }
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': R}, indent=1) + '\n')
    for p in (lab/'public.json', TASK/'tests/public.json'):
        d = json.loads(p.read_text())
        d['solver_budget'] = BUDGET
        p.write_text(json.dumps(d, indent=1) + '\n')
    # ---------------------------------------------------------------- catalog
    catalog = {
        'study': 'J. Zhou, S. Ponce, J.-C. Charlier, High-throughput calculations of spin Hall conductivity in non-magnetic '
                 '2D materials, npj 2D Mater. Appl. 9, 39 (2025); Materials Cloud Archive record rfcdy-wcy32 (CC BY 4.0)',
        'budget': BUDGET,
        'materials': {n: {'archive_folder': meta['materials'][n]['archive_folder'],
                          'staged_files': sorted(f'/app/archive/{n}/{p.name}' for p in (agent/'archive'/n).iterdir()),
                          'num_wann': meta['materials'][n]['num_wann']} for n in MATERIALS},
        'archived_shc_folders': sorted(archived),
        'requests': {
            'bands': {'cost': 1, 'format': {'kind': 'bands', 'material': '<material>', 'k_frac': [[0.0, 0.0], [0.5, 0.0]], 'spin': False},
                      'limits': 'up to 1000 k-points; spin = true adds <sigma_x>, <sigma_y>, <sigma_z> of every band'},
            'band_edges': {'cost': 1, 'format': {'kind': 'band_edges', 'material': '<material>', 'band': 10, 'mesh': 120},
                           'limits': 'mesh 12-300 (Gamma-centred N x N); returns max of band n, min of band n+1 and the smallest direct gap between them, with their k'},
            'wcc': {'cost': '1 if len(k) * n_loop <= 6000, 2 if <= 24000, else 3',
                    'format': {'kind': 'wcc', 'material': '<material>', 'bands': 10, 'loop_axis': 1, 'k': [0.0, 0.25, 0.5], 'n_loop': 60},
                    'limits': 'up to 201 k values, n_loop 20-300; Wilson loop of bands 1..n along reciprocal axis loop_axis at each reduced k of the other axis'},
            'shc': {'cost': 'ceil((N1*N2 + factor^2 * max_points) / 40000) + (number of spin components - 1)',
                    'format': {'kind': 'shc', 'material': '<material>', 'spin': ['z'], 'mesh': [60, 60], 'eta_eV': 0.005,
                               'position': '<lattice or wannier_centres>', 'mu_eV': [-3.5],
                               'refine': {'factor': 8, 'threshold_A2': 10.0, 'max_points': 500}},
                    'limits': 'mesh N1, N2 in 12-300 (Gamma-centred); N1*N2 + factor^2*max_points <= 160000; eta_eV 1e-4 to 0.05; up to 41 mu values; refine optional',
                    'parameters': {
                        'spin': 'spin axes gamma among x, y, z (Cartesian axes of the archived cell)',
                        'eta_eV': 'constant broadening eta in the Kubo denominator (e_n - e_m)^2 + eta^2',
                        'position': 'how the velocity operator places the Wannier functions: "lattice" = at the lattice point of '
                                    'their home cell (the phase convention of aiida_hr.dat); "wannier_centres" = at their final '
                                    'centres in cut_aiida.wout (diagonal position operator)',
                        'mu_eV': 'chemical potentials; T = 0 step occupation',
                        'refine': 'optional adaptive refinement as in Wannier90 berry_curv_adpt_kmesh: coarse points where the '
                                  'occupied spin Berry curvature exceeds threshold_A2 (Angstrom^2) for any requested spin and mu '
                                  '(at most max_points, largest first) are replaced by the average over a factor x factor sub-mesh of their mesh cell'},
                    'output': 'sigma^gamma_xy(mu) in Wannier90 units, (hbar/e) S/cm, computed with the volume of the archived cell '
                              '(area x c, vacuum included), as in the archived aiida-shc-fermiscan.dat files'},
            'ribbon': {'cost': 2, 'format': {'kind': 'ribbon', 'material': '<material>', 'finite_axis': 1, 'width': 40,
                                             'k': [0.0, 0.1, 0.2], 'window_eV': [-3.7, -3.4]},
                       'limits': 'width * num_wann <= 1300; up to 31 k values; window width <= 0.6 eV; open ribbon finite along lattice vector finite_axis'},
            'archived_shc': {'cost': 1, 'format': {'kind': 'archived_shc', 'folder': 'Top_materials/Bi2/SHC-xyz'},
                             'limits': 'returns the archived aiida-shc-fermiscan.dat rows and the aiida.win of that folder'}},
        'conventions': {
            'k': 'reduced coordinates of the reciprocal basis of CELL_PARAMETERS in aiida.in; k_z = 0',
            'energies': 'eV on the energy scale of each archived calculation',
            'bands': 'Wannier bands numbered from 1 (lowest Wannier band) within each model',
            'hamiltonian': 'archived aiida_hr.dat with the use_ws_distance image selection of Wannier90 (as in the archived runs); '
                           'it reproduces the archived aiida_band.dat',
            'spin_operator': 'Pauli matrices on the spinor index of the Wannier functions (ordered 1 up, 1 down, 2 up, ...); '
                             'the archive has no .spn file',
            'spin_berry_curvature': 'Omega^gamma_n = sum_{m != n} -2 Im[<n|j^gamma|m><m|v_y|n>] / ((e_n - e_m)^2 + eta^2), '
                                    'j^gamma = {sigma^gamma, v_x}/2 (Qiao et al., PRB 98, 214402 (2018); Wannier90 berry_task = shc)'},
    }
    (agent/'catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    # ---------------------------------------------------------------- calibration printout
    def fails(abl):
        bad = []
        for n, v in abl.get('plateau_sigma_z', {}).items():
            if not plateau_ok(v, plateau[n]):
                bad.append(f'plateau:{n} {v:.4f} vs {plateau[n]:.4f}')
        for n, v in abl.get('z2', {}).items():
            if v != z2[n]:
                bad.append(f'z2:{n}')
        for n, d in abl.get('spin_direction', {}).items():
            dth = d['theta_deg'] - spin_dir[n]['theta_deg']
            dph = (d['phi_deg'] - spin_dir[n]['phi_deg'] + 180) % 360 - 180
            if abs(dth) > THRESHOLDS['theta_deg'] or abs(dph) > THRESHOLDS['phi_deg']:
                bad.append(f'dir:{n} dth={dth:.1f} dph={dph:.1f}')
        return bad
    report = {'z2': z2, 'plateau': plateau, 'spin_direction': spin_dir, 'facts': facts,
              'ablations_failing_checks': {k: fails(v) for k, v in A.items()},
              'robustness_failing_checks': {k: fails(v) for k, v in R.items()}}
    (HERE/'calibration.json').write_text(json.dumps(report, indent=1, default=float) + '\n')
    print(json.dumps(report, indent=1, default=float))


if __name__ == '__main__':
    main()
