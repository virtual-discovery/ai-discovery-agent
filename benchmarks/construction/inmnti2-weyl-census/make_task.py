"""Generate the staged, lab, private and ablation data of inmnti2-weyl-census from the extracted archive.

    bash benchmarks/construction/inmnti2-weyl-census/run_in_image.sh make_task.py [--reuse-census]

Reads DATA/x (fetch.py), writes into benchmarks/tasks/inmnti2-weyl-census:
  environment/lab/data/{hubbard,standard}.npz   Wannier hr.dat as npz (complex64)
  environment/lab/data/archive.json             archived DFT bands, Kubo outputs, WannierTools outputs
  environment/lab/tbmodel.py                    copy of ./tbmodel.py (the lab's evaluator)
  environment/lab/private/hidden.json, hidden_values.json
  environment/agent/catalog.json, environment/agent/archive/...   staged archive files
  solution/ablation_values.json, the budget in both public.json files
and census.json + calibration.json here. All references use the lab's own evaluator on the lab's npz files.
"""
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import model as M
import tbmodel as T

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/inmnti2-weyl-census'
LAB, AGENT = TASK/'environment/lab', TASK/'environment/agent'
BUDGET = 28
COSTS = {'wannier_bands': '1 (<= 500 k-points)', 'gap_scan': '1 (mesh <= 24) or 2 (mesh 25-36)', 'find_nodes': '1 (<= 8 seeds)',
         'berry_flux': '1 (<= 8 spheres)',
         'dft_bands': 3, 'kubo': 2, 'wannier_tools': 2}
MODELS = {'hubbard': 'hubbard/wannier/InMnTi2_hr.dat', 'standard': 'standard/wannier/InMnTi2_hr.dat'}
DFT_RUNS = ['hubbard/nscf', 'hubbard/nscfXYZ', 'standard/nscf', 'standard/nscfXYZ']
KUBO_RUNS = ['10', '30', '58', '120', '210', '300']
WT_RUNS = {'hubbard/find': ['hubbard/wtool/wt_find.in', 'hubbard/wtool/find/Nodes.dat', 'hubbard/wtool/find/WT.out'],
           'hubbard/weyl': ['hubbard/wtool/wt_weyl.in', 'hubbard/wtool/weyl/WT.out'],
           'standard/find': ['standard/wtool/find/wt_find.in', 'standard/wtool/find/Nodes.dat'],
           'standard/weyl': ['standard/wtool/weyl/wt_weyl.in', 'standard/wtool/weyl/Cherns.dat']}
STAGED = ['README.txt',
          'hubbard/scf/scf.in', 'hubbard/scf/scf.out', 'hubbard/relax/relax.in', 'hubbard/nscf/nscf.in',
          'hubbard/nscfXYZ/nscf.in', 'hubbard/wannier/InMnTi2.win', 'hubbard/wannier/InMnTi2.wout.excerpt',
          'hubbard/wannier/pw2wan.in', 'hubbard/wannier/pw2wan.out', 'hubbard/wannier/kubo/InMnTi2.win',
          'hubbard/wberri/run_wberri.py',
          'standard/scf/scf.in', 'standard/scf/scf.out', 'standard/nscf/nscf.in', 'standard/nscfXYZ/nscf.in',
          'standard/wannier/InMnTi2.win', 'standard/wannier/InMnTi2.wout.excerpt', 'standard/wannier/pw2wan.in',
          'standard/wannier/pw2wan.out', 'standard/wberri/run_wberri.py']
THRESHOLDS = {'node_energy_abs_meV': 3.0, 'pair_separation_rel': 0.02, 'kz_offset_abs_invA': 0.002,
              'pbe_node_energy_abs_meV': 4.0, 'ahc_xy_rel': 0.035, 'ahc_forbidden_abs_S_cm': 1.5}


def relax_excerpt(text):
    """Input echo, every BFGS step's energy/pressure/cell, the final coordinates and the end of the file."""
    lines = text.splitlines()
    head = lines[:lines.index(next(l for l in lines if 'number of k points' in l)) + 1]
    keys = ('!    total energy', 'P=', 'total magnetization', 'absolute magnetization', 'new unit-cell volume',
            'the Fermi energy', 'number of scf cycles', 'number of bfgs steps', 'enthalpy new', 'enthalpy old',
            'bfgs converged', 'End of BFGS', 'Final enthalpy', 'Energy error', 'Cell gradient error', 'convergence has been')
    body = []
    for i, line in enumerate(lines[len(head):], len(head)):
        if any(k in line for k in keys):
            body.append(f'[line {i + 1}] ' + line.rstrip())
        if line.startswith('CELL_PARAMETERS'):
            body += [f'[line {i + 1 + j}] ' + l.rstrip() for j, l in enumerate(lines[i:i + 4])]
    start = next(i for i, l in enumerate(lines) if 'Begin final coordinates' in l)
    end = next(i for i, l in enumerate(lines) if 'End final coordinates' in l)
    return ('\n'.join(head) + '\n\n[... selected lines of the vc-relax trajectory ...]\n' + '\n'.join(body) +
            '\n\n' + '\n'.join(lines[start:end + 1]) + '\n\n[... last 30 lines of the file ...]\n' +
            '\n'.join(lines[-30:]) + '\n')


def wt_text(member):
    text = (M.X/member).read_text()
    if member.endswith('find/WT.out'):
        lines = text.splitlines()
        return '\n'.join(lines[:140]) + '\n[... lines omitted ...]\n' + '\n'.join(lines[-25:]) + '\n'
    if member.endswith('weyl/WT.out'):
        lines = text.splitlines()
        i = next(k for k, l in enumerate(lines) if 'Start of chirality' in l)
        keep = [l for l in lines[i:] if 'We are calculation' not in l and 'Wilson loop ik' not in l
                and 'k (Cartesian' not in l and l.strip()]
        return '\n'.join(lines[:140]) + '\n[... lines omitted ...]\n' + '\n'.join(keep[:40]) + '\n'
    return text


def convert_models():
    info = {}
    for name, member in MODELS.items():
        R, deg, H = M.read_hr(M.X/member)
        np.savez(LAB/'data'/f'{name}.npz', R=R.astype(np.int16), deg=deg.astype(np.int16), H=H.astype(np.complex64))
        info[name] = {'source': member, 'num_wann': int(H.shape[1]), 'nrpts': int(len(R))}
    return info


def dft_payload(run):
    model, kind = run.split('/')
    inp = (M.X/run/'nscf.in').read_text()
    k, e = M.qe_bands(M.X/run/'nscf.out')
    return {'run': run, 'files': [f'{run}/nscf.in', f'{run}/nscf.out'], 'nscf.in': inp,
            'nscf.out (selected lines)': M.qe_header(M.X/run/'nscf.out'),
            'k_points_2pi_over_alat': k.tolist(), 'eigenvalues_eV': e.tolist(),
            'note': 'k-points and band energies (eV) as printed in nscf.out, all bands, one row per k-point'}


def kubo_payload():
    out = {}
    for run in KUBO_RUNS:
        folder = M.X/'hubbard/wannier/kubo'/run
        out[run] = {}
        for f in sorted(folder.iterdir()):
            rows = np.loadtxt(f)
            rows = rows[rows[:, 0] <= 1.5 + 1e-9]
            comp = f.name.replace('InMnTi2-', '').replace('.dat', '').replace('kubo_', '')
            out[run][comp] = {'file': f'hubbard/wannier/kubo/{run}/{f.name}', 'rows': [[float('%.8g' % x) for x in r] for r in rows]}
    return out


SHIFTS = np.array([(i, j, k) for i in (-2, -1, 0, 1, 2) for j in (-2, -1, 0, 1, 2) for k in (-2, -1, 0, 1, 2)]) @ T.B


def ws_fold(k):
    """Representative of Cartesian k in the first Brillouin zone (smallest |k + G|)."""
    cand = np.asarray(k) + SHIFTS
    return cand[np.argmin(np.linalg.norm(cand, axis=1))]


def census(mdl, n, mesh, ef, flux=True):
    rows, nmin, gmin = mdl.gap_scan(n, mesh, keep=400)
    found = []
    for r in rows:
        x, g, e, _ = mdl.minimize_gap(np.array(r['k_cart_invA']), n)
        kf = x @ np.linalg.inv(T.B)
        kf -= np.floor(kf + 0.5)
        if any(np.linalg.norm(((kf - f['kf']) - np.round(kf - f['kf'])) @ T.B) < 2e-3 for f in found):
            continue
        found.append({'kf': kf, 'k': ws_fold(kf @ T.B), 'gap': g, 'E': e - ef})
    nodes = [f for f in found if f['gap'] < 1e-6]
    near = [f for f in found if 1e-6 <= f['gap'] < 0.01]
    if flux:
        for f in nodes:
            f['chi'] = round(mdl.berry_flux(f['k'], 0.01, n, 12, 24)[0])
    clean = lambda f: {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in f.items()}
    return {'grid_minima': nmin, 'nodes': [clean(f) for f in nodes], 'avoided': [clean(f) for f in near]}


def dft_crossings():
    """Crossings of DFT bands 52/53 on the three archived DFT+U lines (hubbard/nscfXYZ)."""
    kc, pts = M.tpiba_b_path((M.X/'hubbard/nscfXYZ/nscf.in').read_text())
    kp, e = M.qe_bands(M.X/'hubbard/nscfXYZ/nscf.out')
    assert np.abs(kc / M.TPIBA - kp).max() < 1e-4
    d = (e - M.EF['hubbard']) * 1000
    out = {}
    for name, seg in M.xyz_segments(pts):
        ks = kc[seg]; s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(ks, axis=0), axis=1))]
        lo, hi = d[seg, M.NOCC_DFT - 1], d[seg, M.NOCC_DFT]
        g = hi - lo
        mins = [i for i in range(1, len(g) - 1) if g[i] < g[i - 1] and g[i] <= g[i + 1] and g[i] < 10]
        for i0 in mins:
            for w in (0.02, 0.03, 0.05):
                sc, ec, sl = M.crossing_on_line(s, lo / 1000, hi / 1000, i0, w)
                k = ks[0] + (ks[-1] - ks[0]) * sc / s[-1]
                out[f'{name}{i0}_w{w}'] = {'line': name, 'index': i0, 'window': w, 's': sc, 'E_meV': ec * 1000,
                                           'k': k.tolist(), 'gap_min_meV': float(g[i0]), 'slopes': sl}
    return out, kc, d


def main():
    reuse = '--reuse-census' in sys.argv
    (LAB/'data').mkdir(parents=True, exist_ok=True)
    (LAB/'private').mkdir(parents=True, exist_ok=True)
    for f in ('hubbard/scf/scf.out', 'standard/scf/scf.out'):
        assert f'the Fermi energy is    {M.EF[f.split("/")[0]]:.4f} ev' in (M.X/f).read_text()
    assert 'mp_grid = 9 9 9' in (M.X/'hubbard/wannier/InMnTi2.win').read_text()
    assert 'mp_grid = 12 12 12' in (M.X/'standard/wannier/InMnTi2.win').read_text()
    assert 'starting_magnetization' not in (M.X/'standard/scf/scf.in').read_text()
    model_info = convert_models()
    shutil.copy2(HERE/'tbmodel.py', LAB/'tbmodel.py')
    hub, std = T.Model(LAB/'data/hubbard.npz'), T.Model(LAB/'data/standard.npz')

    # --- band-index offsets: Wannier bands equal the DFT bands at the scf k-points (fitting grid) ---
    for name, mdl in (('hubbard', hub), ('standard', std)):
        k, e = M.qe_bands(M.X/name/'scf/scf.out')
        w = mdl.bands(k * M.TPIBA)
        off = M.OFFSET[name]
        dev = np.abs(w[:, 10:20] - e[:, off + 10:off + 20]).max()
        print(name, 'max |Wannier - DFT| at the scf k-points (bands 11-20):', round(dev * 1000, 2), 'meV')
        assert dev < (2e-3 if name == 'hubbard' else 1e-2), (name, dev)   # hubbard scf grid = its 9^3 fitting grid

    # --- archive payloads for the lab ---
    archive = {'dft': {run: dft_payload(run) for run in DFT_RUNS}, 'kubo': kubo_payload(),
               'wannier_tools': {run: {m: wt_text(m) for m in members} for run, members in WT_RUNS.items()},
               'models': model_info}
    (LAB/'data/archive.json').write_text(json.dumps(archive))

    # --- staged files ---
    if (AGENT/'archive').exists():
        shutil.rmtree(AGENT/'archive')
    for member in STAGED:
        target = AGENT/'archive'/member
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(M.X/member, target)
    (AGENT/'archive/hubbard/relax').mkdir(parents=True, exist_ok=True)
    (AGENT/'archive/hubbard/relax/relax.out.excerpt').write_text(relax_excerpt((M.X/'hubbard/relax/relax.out').read_text()))
    (AGENT/'archive/files_description.md').write_text((M.DATA/'files_description.md').read_text())

    # --- censuses with the lab evaluator ---
    cpath = HERE/'census.json'
    if reuse and cpath.exists():
        C = json.loads(cpath.read_text())
    else:
        C = {'hubbard_16': census(hub, 16, 40, M.EF['hubbard']),
             'hubbard_18': census(hub, 18, 40, M.EF['hubbard'], flux=False),
             'standard_18': census(std, 18, 36, M.EF['standard']),
             'standard_16': census(std, 16, 36, M.EF['standard'], flux=False)}
        coarse = std.gap_scan(18, 20, keep=400)[0]
        C['standard_18_coarse20_minima_below_50meV'] = sum(r['gap_eV'] < 0.05 for r in coarse)
        cpath.write_text(json.dumps(C, indent=1))
    for key in ('hubbard_16', 'hubbard_18', 'standard_18', 'standard_16'):
        for f in C[key]['nodes'] + C[key]['avoided']:
            f['k'] = ws_fold(f['k']).tolist()
    H16, S18 = C['hubbard_16']['nodes'], C['standard_18']['nodes']
    assert len(H16) == 24 and len(S18) == 24 and not C['hubbard_16']['avoided']
    assert sum(f['chi'] for f in H16) == 0 and sum(f['chi'] for f in S18) == 0

    # --- magnetic symmetry orbits ---
    ops_h = M.symmetry_ops((M.X/'hubbard/scf/scf.out').read_text())
    ops_s = M.symmetry_ops((M.X/'standard/scf/scf.out').read_text())
    assert len(ops_h) == 8 and sum(t for _, t in ops_h) == 4 and len(ops_s) == 24
    kh = np.array([f['k'] for f in H16]); ch = np.array([f['chi'] for f in H16])
    orb = M.orbits(kh, ch, ops_h, T.B, tol=0.05)
    orbit_sizes = sorted(len(g) for g in orb)
    assert orbit_sizes == [4, 4, 4, 4, 8]
    orb_unitary = sorted(len(g) for g in M.orbits(kh, ch, [o for o in ops_h if o[1] == 0], T.B, tol=0.05))
    orb_std = sorted(len(g) for g in M.orbits(np.array([f['k'] for f in S18]), np.array([f['chi'] for f in S18]),
                                                ops_s, T.B, tol=0.03, nonmagnetic=True))
    fam_a = [H16[i] for i in next(g for g in orb if len(g) == 8)]
    assert all(abs(f['k'][2]) < 0.03 for f in fam_a)
    planar = sorted([sum(abs(f['k'][a]) < 0.03 for f in H16) for a in range(3)])

    # --- DFT truth on the archived lines ---
    X, kline, dline = dft_crossings()
    xs = sorted([v for v in X.values() if v['line'] == 'x' and v['window'] == 0.03], key=lambda v: v['s'])
    yv = next(v for v in X.values() if v['line'] == 'y' and v['window'] == 0.03)
    zv = next(v for v in X.values() if v['line'] == 'z' and v['window'] == 0.03)
    assert len(xs) == 2 and all(v['gap_min_meV'] < 5 for v in (xs[0], xs[1], yv, zv))
    E_ref = float(np.mean([xs[0]['E_meV'], yv['E_meV'], zv['E_meV']]))
    sep_ref = xs[1]['s'] - xs[0]['s']
    node1, node2 = np.array(xs[0]['k']), np.array(xs[1]['k'])          # node1 at kz < 0
    chi1 = round(hub.berry_flux(node1, 0.03, 16, 16, 32)[0])
    chi2 = round(hub.berry_flux(node2, 0.03, 16, 16, 32)[0])
    assert chi1 == -chi2 and abs(chi1) == 1
    plus_node_kz_sign = np.sign(node1[2]) if chi1 == 1 else np.sign(node2[2])
    kz_ref = float(plus_node_kz_sign * abs(zv['k'][2]))
    assert node1[2] < 0 and chi1 == 1 and kz_ref < 0
    pbe_ref = float(np.mean([f['E'] for f in S18]) * 1000)

    # --- Kubo ---
    K = archive['kubo']
    a0 = {run: {c: K[run][f'A_{c}']['rows'][0][1] for c in ('xy', 'yz', 'zx')} for run in KUBO_RUNS}
    assert all(K[run][f'A_{c}']['rows'][0][0] == 0 for run in KUBO_RUNS for c in ('xy', 'yz', 'zx'))
    ahc_ref = a0['300']['xy']
    assert abs(a0['210']['xy'] / ahc_ref - 1) < 0.01

    # --- Wannier-model paths (wrong) ---
    def nearest_opposite(nodes, f):
        return min(np.linalg.norm(np.array(g['k']) - np.array(f['k'])) for g in nodes if g['chi'] == -f['chi'])
    wan_E = float(np.mean([f['E'] for f in fam_a]) * 1000)
    wan_sep = float(np.mean([nearest_opposite(fam_a, f) for f in fam_a]))
    wan_kz = float(np.mean([f['k'][2] for f in fam_a if f['chi'] == 1]))
    near1 = min(H16, key=lambda f: np.linalg.norm(np.array(f['k']) - node1))
    near2 = min(H16, key=lambda f: np.linalg.norm(np.array(f['k']) - node2))
    grid_E = float((near1['E'] + near2['E']) / 2 * 1000)
    grid_sep = float(np.linalg.norm(np.array(near1['k']) - np.array(near2['k'])))
    std_sep = float(np.mean([nearest_opposite(S18, f) for f in S18]))
    ef_line = float([l for l in archive['dft']['standard/nscfXYZ']['nscf.out (selected lines)'] if 'Fermi energy' in l][0].split()[-2])

    hidden = {'targets': {'node_count': {'dftu': 24, 'pbe': 24}, 'orbit_sizes': orbit_sizes,
                          'node_energy_meV': round(E_ref, 2), 'pair_separation_invA': round(sep_ref, 5),
                          'kz_offset_plus_invA': round(kz_ref, 5), 'pbe_node_energy_meV': round(pbe_ref, 2),
                          'ahc_S_per_cm': {'xy': round(ahc_ref, 3), 'yz': 0.0, 'zx': 0.0}},
              'thresholds': THRESHOLDS}
    (LAB/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (LAB/'private/hidden_values.json').write_text(json.dumps(
        [round(E_ref, 2), round(sep_ref, 5), round(kz_ref, 5), round(pbe_ref, 2), round(ahc_ref, 3), sep_ref, E_ref, pbe_ref,
         ahc_ref, kz_ref]) + '\n')

    abl = {'standard_model_as_magnetic': {'node_energy_meV': pbe_ref, 'pair_separation_invA': std_sep,
                                          'kz_offset_plus_invA': 0.0, 'orbit_sizes': orb_std},
           'band_index_from_exclude_bands': {'node_count': {'dftu': len(C['hubbard_18']['nodes']), 'pbe': 24}},
           'same_band_index_both_models': {'node_count': {'dftu': 24, 'pbe': len(C['standard_16']['nodes'])},
                                           'pbe_node_energy_meV': float(np.mean([f['E'] for f in C['standard_16']['nodes']]) * 1000)},
           'wannier_nodes_trusted': {'node_energy_meV': wan_E, 'pair_separation_invA': wan_sep, 'kz_offset_plus_invA': wan_kz},
           'grid_validated_nearest_pair': {'node_energy_meV': grid_E, 'pair_separation_invA': grid_sep},
           'cubic_single_orbit': {'orbit_sizes': orb_std},
           'unitary_operations_only': {'orbit_sizes': orb_unitary},
           'planar_families': {'orbit_sizes': planar},
           'magnetization_axis_split': {'orbit_sizes': [8, 16]},
           'chirality_of_upper_band': {'kz_offset_plus_invA': -kz_ref},
           'cubic_symmetric_node': {'kz_offset_plus_invA': 0.0},
           'pbe_count_with_avoided_crossings': {'node_count': {'dftu': 24, 'pbe': 24 + len(C['standard_18']['avoided'])}},
           'pbe_count_coarse_grid': {'node_count': {'dftu': 24, 'pbe': C['standard_18_coarse20_minima_below_50meV']}},
           'kubo_folder_as_broadening': {'ahc_S_per_cm': {'xy': a0['10']['xy'], 'yz': 0.0, 'zx': 0.0}},
           'kubo_mesh_58': {'ahc_S_per_cm': {'xy': a0['58']['xy'], 'yz': 0.0, 'zx': 0.0}},
           'raw_hall_components': {'ahc_S_per_cm': {'xy': a0['300']['xy'], 'yz': a0['300']['yz'], 'zx': a0['300']['zx']}},
           'nscf_path_fermi_level': {'pbe_node_energy_meV': pbe_ref - (ef_line - M.EF['standard']) * 1000}}
    rob = {'node_energy_y_line': {'node_energy_meV': yv['E_meV']},
           'node_energy_x_line': {'node_energy_meV': xs[0]['E_meV']},
           'node_energy_window_extremes': {'node_energy_meV': min(v['E_meV'] for v in X.values())},
           'separation_from_x_and_z_crossings': {'pair_separation_invA': float(np.linalg.norm(np.array(zv['k']) - node1))},
           'kz_from_x_line_passage': {'kz_offset_plus_invA': float(node1[2])},
           'ahc_mesh_210': {'ahc_S_per_cm': {'xy': a0['210']['xy'], 'yz': 0.0, 'zx': 0.0}},
           'pbe_single_node_extremes': {'pbe_node_energy_meV': float(min(f['E'] for f in S18) * 1000)}}
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    for p in (LAB/'public.json', TASK/'tests/public.json'):
        d = json.loads(p.read_text()); d['solver_budget'] = BUDGET; p.write_text(json.dumps(d, indent=1) + '\n')

    catalog = {
        'material': 'InMnTi2, inverse full-Heusler (space group 216), noncollinear DFT with spin-orbit coupling',
        'source': 'Materials Cloud Archive record hkprp-xw234 (DOI 10.24435/materialscloud:ph-3c, CC BY 4.0): '
                  'D. Grassano, L. Binci, N. Marzari, input/output files of the InMnTi2 Weyl-semimetal study',
        'staged_files': {'root': '/app/archive', 'files': sorted(str(p.relative_to(AGENT/'archive')) for p in (AGENT/'archive').rglob('*') if p.is_file()),
                         'note': 'Archive files copied unchanged, except *.excerpt (selected lines of large outputs).'},
        'wannier_models': {name: {**info, 'fermi_energy_eV_from_scf': M.EF[name]} for name, info in model_info.items()},
        'dft_runs_purchasable': DFT_RUNS,
        'kubo_runs_purchasable': {'folders': KUBO_RUNS, 'location': 'hubbard/wannier/kubo/<folder>/',
                                  'components': {run: sorted(K[run]) for run in KUBO_RUNS},
                                  'served_range': 'rows with omega <= 1.5 eV of each archived file (columns as in the file)'},
        'wannier_tools_runs_purchasable': {run: members for run, members in WT_RUNS.items()},
        'costs': COSTS, 'budget': BUDGET,
        'conventions': {'k_cart_invA': 'Cartesian k in 1/Angstrom including 2 pi, axes of CELL_PARAMETERS',
                        'k_frac': 'coordinates in the reciprocal basis b_i of CELL_PARAMETERS (k_cart = k_frac @ B)',
                        'reciprocal_basis_B_invA': T.B.round(8).tolist(),
                        'wannier_band_numbering': 'Wannier bands numbered from 1 (lowest Wannier band) in each model',
                        'energies': 'all lab energies are absolute eV on the scale of the model / run that produced them'}}
    (AGENT/'catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')

    cal = {'node_energy_meV': {'reference': E_ref, 'definition': 'mean of the DFT 52/53 crossings on the x, y, z lines (fit window 0.03 1/A)',
                               'robustness': {k: v['E_meV'] for k, v in X.items()},
                               'nearest_ablation': {'grid_validated_nearest_pair': grid_E, 'wannier_nodes_trusted': wan_E,
                                                    'standard_model': pbe_ref},
                               'wannier_family_energies_meV': sorted(f['E'] * 1000 for f in fam_a), 'tolerance': '+-3 meV'},
           'pair_separation_invA': {'reference': sep_ref, 'robustness': {'x_and_z_crossings': rob['separation_from_x_and_z_crossings']['pair_separation_invA']},
                                    'nearest_ablation': {'grid_validated_nearest_pair': grid_sep, 'wannier_nodes_trusted': wan_sep,
                                                         'standard_model': std_sep}, 'tolerance': '+-2%'},
           'kz_offset_plus_invA': {'reference': kz_ref, 'robustness': {'x_line_passage': float(node1[2])},
                                   'nearest_ablation': {'cubic_symmetric_node': 0.0, 'chirality_of_upper_band': -kz_ref},
                                   'wannier_plus_family_kz': sorted(f['k'][2] for f in fam_a if f['chi'] == 1),
                                   'not_ablation': {'paper_0.0047_read_as_invA': -0.0047}, 'tolerance': '+-0.002 1/A'},
           'pbe_node_energy_meV': {'reference': pbe_ref, 'spread': [min(f['E'] for f in S18) * 1000, max(f['E'] for f in S18) * 1000],
                                   'nearest_ablation': {'nscf_path_fermi_level': abl['nscf_path_fermi_level']['pbe_node_energy_meV'],
                                                        'same_band_index_both_models': abl['same_band_index_both_models']['pbe_node_energy_meV']},
                                   'tolerance': '+-4 meV'},
           'ahc_S_per_cm': {'reference_xy_300': ahc_ref, 'by_folder': a0, 'tolerance': 'xy +-3.5%; |yz|, |zx| <= 1.5 S/cm (symmetry: 0)'},
           'node_count': {'reference': {'dftu': 24, 'pbe': 24}, 'nearest_ablation': {k: abl[k]['node_count'] for k in abl if 'node_count' in abl[k]}},
           'orbit_sizes': {'reference': orbit_sizes, 'nearest_ablation': {'cubic': orb_std, 'unitary_only': orb_unitary, 'planar': planar,
                                                                          'm_axis_split': [8, 16]}},
           'dft_crossings': X, 'wannier_family_A': fam_a}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1, default=float) + '\n')
    print(json.dumps({'E_ref': E_ref, 'sep_ref': sep_ref, 'kz_ref': kz_ref, 'pbe_ref': pbe_ref, 'ahc_ref': ahc_ref,
                      'orbit_sizes': orbit_sizes, 'orb_unitary': orb_unitary, 'orb_std': orb_std, 'planar': planar,
                      'wan_E': wan_E, 'wan_sep': wan_sep, 'wan_kz': wan_kz, 'grid_E': grid_E, 'grid_sep': grid_sep,
                      'std_sep': std_sep, 'ef_line': ef_line, 'a0': a0,
                      'counts': {k: len(v['nodes']) for k, v in C.items() if isinstance(v, dict)},
                      'avoided_std': len(C['standard_18']['avoided']), 'coarse': C['standard_18_coarse20_minima_below_50meV']},
                     indent=1, default=float))


if __name__ == '__main__':
    main()
