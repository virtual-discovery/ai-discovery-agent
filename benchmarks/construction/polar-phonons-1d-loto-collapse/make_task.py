"""Generate the staged, lab, private and ablation data of polar-phonons-1d-loto-collapse.

    bash run_in_image.sh make_task.py          # archive unpacked by fetch.py, mounted at /data

Writes (task folder benchmarks/tasks/polar-phonons-1d-loto-collapse):
  environment/agent/catalog.json             what can be bought, archive inputs, README
  environment/lab/data/index.json + blobs/   served archive files (gzip), parsed interpolated dispersions
  environment/lab/private/hidden.json        withheld explicit small-q LO frequencies and the 3D q->0 LO
  environment/lab/private/hidden_values.json
  solution/ablation_values.json              wrong paths and defensible variants, from the archived data
  environment/lab/public.json, tests/public.json (budget)
and calibration.json here.
"""
import gzip
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/polar-phonons-1d-loto-collapse'
ROOT = Path(os.environ.get('PP1D_ARCHIVE', '/data/materialscloud_data_Rivano'))
BUDGET = 30
COSTS = {'dfpt': 1, 'interpolated': 2, 'archived_file': 1}
MAX_TEXT = 400_000
TOL = {'0.01': 2.0, '0.025': 3.0, '0.05': 4.5, '0.075': 4.5, '0.1': 4.5, '0.125': 4.5}
TOL_3D = 3.0

SYSTEMS = {
    'BN_atomic_chain': 'BN atomic chain (polyyne-like B-N chain), 2 atoms per cell',
    'BN_nanotubes/44': 'armchair BN nanotube (4,4), 16 atoms per cell',
    'BN_nanotubes/55': 'armchair BN nanotube (5,5), 20 atoms per cell [hidden-target system]',
    'BN_nanotubes/66': 'armchair BN nanotube (6,6), 24 atoms per cell [hidden-target system]',
    'GaAs_nanowire': 'H-passivated wurtzite GaAs nanowire, 24 atoms per cell (the archived titles say 84)',
    'GaAs_bulk/ZB_phase': 'bulk zincblende GaAs (LDA)',
    'GaAs_bulk/WZ_phase': 'bulk wurtzite GaAs (LDA)',
}
HIDDEN_TUBES = {'BN55': 'BN_nanotubes/55', 'BN66': 'BN_nanotubes/66'}
TARGETS = {'BN55': {'0.01': 'dyn100', '0.025': 'dyn25', '0.05': 'dyn50', '0.075': 'dyn75', '0.125': 'dyn125'},
           'BN66': {'0.01': 'dyn100', '0.025': 'dyn25', '0.05': 'dyn50', '0.075': 'dyn75', '0.1': 'dyn10',
                    '0.125': 'dyn125'}}
WITHHELD_DYN = ({f'BN_nanotubes/55/1D_code/longwave/{n}' for n in ('dyn10', 'dyn100', 'dyn25', 'dyn50', 'dyn75', 'dyn125')}
                | {f'BN_nanotubes/66/1D_code/longwave/{n}' for n in ('dyn10', 'dyn100', 'dyn25', 'dyn50', 'dyn75', 'dyn125')}
                | {f'BN_nanotubes/66/1D_code/{n}' for n in ('dyn10', 'dyn100', 'dyn125', 'dyn15', 'dyn25', 'dyn50', 'dyn75')})
WITHHELD_FOLDERS_INTERP = {'BN_nanotubes/55/1D_code', 'BN_nanotubes/55/3D_code', 'BN_nanotubes/66/1D_code',
                           'BN_nanotubes/66/3D_code'}
SRC = 'retrieved author-archived file (Materials Cloud 08gsg-gpd20); no new calculation'


def folders():
    out = []
    for sysname in SYSTEMS:
        base = ROOT/sysname
        subs = [p for p in sorted(base.iterdir()) if p.is_dir() and p.name in ('1D_code', '3D_code')]
        out += [str(p.relative_to(ROOT)) for p in subs] if subs else [sysname]
    return out


def bc_of(folder):
    return '1D' if folder.endswith('1D_code') else '3D' if folder.endswith('3D_code') else 'bulk'


def first_q(path):
    text = path.read_text()
    m = re.search(r'q = \(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\)', text)
    return [float(m.group(k)) for k in (1, 2, 3)], text.count('Dynamical  Matrix in cartesian axes')


def strip_positions(text):
    """Drop the ATOMIC_POSITIONS block of a pw.x input (the positions are in every dyn-file header)."""
    out, skip = [], False
    for line in text.split('\n'):
        if line.strip().upper().startswith('ATOMIC_POSITIONS'):
            out.append(line.split()[0] + ' ... [atomic positions omitted in the catalog; see the dyn-file headers]')
            skip = True
            continue
        if skip and (not line.strip() or line.strip().upper().startswith(('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES'))):
            skip = False
        if not skip:
            out.append(line)
    return '\n'.join(out)


def store(name, payload_bytes):
    blob = TASK/'environment/lab/data/blobs'/(hashlib.sha256(name.encode()).hexdigest()[:20] + '.gz')
    blob.parent.mkdir(parents=True, exist_ok=True)
    blob.write_bytes(gzip.compress(payload_bytes, mtime=0))
    return blob.name


def interpolated_payload(folder):
    s_g = M.parse_dyn(ROOT/folder/'dyn1')
    nat = s_g['nat']
    modes = M.parse_modes(ROOT/folder/'matdyn.modes', nat)
    alat = s_g['celldm'][0]
    qz = [round(float(m[0][2]), 6) for m in modes]
    freqs = [[round(float(f), 3) for f in m[1]] for m in modes]
    zw = []
    for q, fr, ev in modes:
        w = np.sum(np.abs(ev[:, :, 2]) ** 2, axis=1) / np.sum(np.abs(ev) ** 2, axis=(1, 2))
        zw.append([round(float(x), 4) for x in w])
    if len(qz) > 400:          # the chain path has 2144 points: keep q <= 0.3 |Gamma-Z| densely, then every 10th
        zmax = max(abs(x) for x in qz)
        keep = [i for i, x in enumerate(qz) if abs(x) <= 0.3 * zmax or i % 10 == 0 or i == len(qz) - 1]
        qz, freqs, zw = [qz[i] for i in keep], [freqs[i] for i in keep], [zw[i] for i in keep]
    return {'folder': folder, 'archive_path': folder + '/matdyn.modes',
            'matdyn_in': (ROOT/folder/'matdyn.in').read_text(),
            'q_z_2pi_over_alat': qz, 'alat_bohr': alat,
            'frequencies_cm-1': freqs, 'z_displacement_weight': zw,
            'note': 'Frequencies (cm^-1) and, for each mode, the z (axial) share of the normalized displacement pattern, '
                    'parsed from the archived matdyn.modes (Gamma to Z path, q_in_band_form). Full eigenvectors are not served.',
            'source': SRC}


def build_catalog_and_lab():
    shutil.rmtree(TASK/'environment/lab/data', ignore_errors=True)
    index, items, interp_items, file_items, inputs = {}, [], [], [], {}
    digests = {}
    for folder in folders():
        fpath = ROOT/folder
        inputs[folder] = {p.name: strip_positions(p.read_text()) for p in sorted(fpath.iterdir())
                          if p.is_file() and p.suffix == '.in'}
        dyns = sorted([p for p in fpath.rglob('dyn*') if p.name != 'dyn0'],
                      key=lambda p: (len(p.relative_to(fpath).parts), p.name))
        for p in dyns:
            rel = str(p.relative_to(ROOT))
            if rel in WITHHELD_DYN:
                continue
            data = p.read_bytes()
            h = hashlib.sha256(data).hexdigest()
            if h in digests:          # identical archived copy: list it under the first path (both names are accepted)
                items[digests[h]]['also_archived_as'].append(rel)
                index[rel] = index[items[digests[h]]['id']]
                continue
            q, nstar = first_q(p)
            digests[h] = len(items)
            items.append({'id': rel, 'system': folder.split('/1D_code')[0].split('/3D_code')[0],
                          'boundary_conditions': bc_of(folder), 'q_recorded_2pi_over_alat': q,
                          'q_points_in_star': nstar, 'also_archived_as': [], 'cost': COSTS['dfpt']})
            index[rel] = {'kind': 'dfpt', 'blob': store(rel, data)}
        if (fpath/'matdyn.modes').exists() and folder not in WITHHELD_FOLDERS_INTERP:
            payload = interpolated_payload(folder)
            index['interpolated:' + folder] = {'kind': 'interpolated', 'blob': store('interp:' + folder, json.dumps(payload).encode())}
            interp_items.append({'folder': folder, 'archive_path': folder + '/matdyn.modes', 'cost': COSTS['interpolated']})
        for p in sorted(fpath.iterdir()):
            if p.is_file() and p.suffix == '.out' and not (p.name == 'matdyn.out' and folder in WITHHELD_FOLDERS_INTERP):
                rel = str(p.relative_to(ROOT))
                text = p.read_text(errors='replace')
                truncated = len(text) > MAX_TEXT
                body = {'path': rel, 'text': text[:MAX_TEXT], 'truncated': truncated, 'bytes': len(text), 'source': SRC}
                index['file:' + rel] = {'kind': 'archived_file', 'blob': store('file:' + rel, json.dumps(body).encode())}
                file_items.append({'path': rel, 'bytes': len(text), 'served_bytes': min(len(text), MAX_TEXT), 'cost': COSTS['archived_file']})
    for it in items:
        if not it['also_archived_as']:
            del it['also_archived_as']
    return index, items, interp_items, file_items, inputs


def lw_truth(folder, name, s):
    d = M.parse_dyn(ROOT/folder/'longwave'/name)
    qz = d['blocks'][0]['q'][2] * 2 * np.pi / s['alat']
    return qz / (np.pi / s['c']), M.lo(s, d['blocks'][0]['D'])


def archived_lo_curve(folder, s):
    """LO (max longitudinal mode charge, Born charges of the folder) along an archived matdyn.modes path."""
    modes = M.parse_modes(ROOT/folder/'matdyn.modes', s['nat'])
    qmax = max(abs(m[0][2]) for m in modes)
    qf = np.array([abs(m[0][2]) / qmax for m in modes])
    lo = []
    for q, fr, ev in modes:
        Pz = np.abs(np.einsum('aj,naj->n', s['Z'][:, 2, :], ev))
        lo.append(fr[np.argmax(Pz)])
    return qf, np.array(lo)


def calibrate_t_ratio(s44, pts44):
    """Oracle protocol: t/R of the released (4,4) tube from its explicit 1D runs, t constrained to [R, t(eps1d=1)]."""
    R = M.atomic_radius(s44)
    tmax = np.sqrt(s44['area'] * (s44['eps'][2, 2] - 1) / np.pi)
    best = None
    for t in np.arange(R, tmax, 0.05):
        ip = M.Interpolator(s44, lr={'t': t})
        r = np.sqrt(np.mean([(M.lo(s44, ip.D(qf * np.pi / s44['c'])) - f) ** 2 for qf, f in pts44]))
        if best is None or r < best[1]:
            best = (t, r)
    return best[0] / R, best


def predict(s, qfs, lr=None, gamma_only=False, qscale=1.0):
    GZ = np.pi / s['c']
    if gamma_only:
        return {k: M.lo(s, s['Dq'][0] + M.na_1d(s, float(k) * qscale * GZ, nG=0, ewald=False, **(lr or {}))) for k in qfs}
    ip = M.Interpolator(s, lr=lr)
    return {k: M.lo(s, ip.D(float(k) * qscale * GZ)) for k in qfs}


def main():
    # ---- verify the 1D kernel against the authors' archived interpolations -------------------------------------
    checks = {}
    for folder in ('BN_nanotubes/55/1D_code', 'BN_nanotubes/66/1D_code', 'GaAs_nanowire/1D_code', 'BN_atomic_chain/1D_code'):
        s = M.load_system(ROOT/folder)
        ip = M.Interpolator(s, lr={})
        modes = M.parse_modes(ROOT/folder/'matdyn.modes', s['nat'])
        qmax = max(abs(m[0][2]) for m in modes)
        step = 1 if len(modes) < 200 else 43
        dev = []
        for q, fr, ev in modes[1:60:step] if step == 1 else modes[1:1100:step]:
            f, Pz, _ = M.modes(s, ip.D(q[2] * 2 * np.pi / s['alat']))
            flo = f[np.argmax(Pz)]
            dev.append(abs(fr[np.argmin(np.abs(fr - flo))] - flo))
        checks[folder] = max(dev)
        assert max(dev) < 0.6, (folder, max(dev))
    # The archived (4,4) "1D" interpolation equals plain Fourier interpolation (no long-range term): record it.
    s44 = M.load_system(ROOT/'BN_nanotubes/44/1D_code')
    a44 = M.parse_modes(ROOT/'BN_nanotubes/44/1D_code/matdyn.modes', s44['nat'])
    ip0 = M.Interpolator(s44)
    dev0 = []
    for q, fr, ev in a44[1:16]:
        flo = M.lo(s44, ip0.D(q[2] * 2 * np.pi / s44['alat']))
        dev0.append(abs(fr[np.argmin(np.abs(fr - flo))] - flo))
    checks['BN_nanotubes/44/1D_code archived interpolation vs plain Fourier interpolation (max dev, q<=0.15)'] = max(dev0)
    assert max(dev0) < 1.5
    # ---- catalog and lab data -----------------------------------------------------------------------------------
    index, items, interp_items, file_items, inputs = build_catalog_and_lab()
    (TASK/'environment/lab/data').mkdir(parents=True, exist_ok=True)
    (TASK/'environment/lab/data/index.json').write_text(json.dumps(index, indent=0) + '\n')
    # ---- hidden truth -------------------------------------------------------------------------------------------
    sysd = {k: M.load_system(ROOT/v/'1D_code') for k, v in HIDDEN_TUBES.items()}
    sys3 = {k: M.load_system(ROOT/v/'3D_code') for k, v in HIDDEN_TUBES.items()}
    truth, purch015 = {}, {}
    for tube, tg in TARGETS.items():
        folder = HIDDEN_TUBES[tube] + '/1D_code'
        truth[tube] = {}
        for key, name in tg.items():
            qf, f = lw_truth(folder, name, sysd[tube])
            assert abs(qf - float(key)) < 1e-4, (tube, key, qf)
            truth[tube][key] = f
        purch015[tube] = lw_truth(folder, 'dyn15', sysd[tube])[1]
    qf3, lo3 = archived_lo_curve('BN_nanotubes/55/3D_code', sys3['BN55'])
    truth_3d = float(lo3[0])
    # ---- released (4,4) explicit runs used by the oracle protocol -------------------------------------------------
    pts44_all = {}
    for name in ('dyn25', 'dyn_05', 'dyn75', 'dyn10', 'dyn125', 'dyn15'):
        qf, f = lw_truth('BN_nanotubes/44/1D_code', name, s44)
        pts44_all[round(qf, 4)] = f
    oracle_pts = [(q, pts44_all[q]) for q in (0.025, 0.05, 0.1)]
    ratio, best44 = calibrate_t_ratio(s44, oracle_pts)
    G44 = M.lo(s44, s44['Dq'][0])
    # ---- predictions --------------------------------------------------------------------------------------------
    def table(fn):
        return {tube: fn(tube, sysd[tube], list(TARGETS[tube])) for tube in TARGETS}
    R = {k: M.atomic_radius(v) for k, v in sysd.items()}
    t_or = {k: ratio * R[k] for k in sysd}
    oracle = table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb]}))
    robust = {
        'authors_radius_d_over_4': table(lambda tb, s, q: predict(s, q, lr={})),
        'no_ewald_damping': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb], 'ewald': False})),
        'ewald_alpha_2': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb], 'alph': 2.0})),
        'radius_1p5_atomic': table(lambda tb, s, q: predict(s, q, lr={'t': 1.5 * R[tb]})),
    }
    gam = {tb: M.lo(s, s['Dq'][0]) for tb, s in sysd.items()}
    abl = {
        'standard_fourier_interpolation': table(lambda tb, s, q: predict(s, q)),
        'atomic_shell_radius': table(lambda tb, s, q: predict(s, q, lr={'t': R[tb]})),
        'cell_dielectric_as_cylinder_dielectric': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb], 'eps1d': s['eps'][2, 2]})),
        'asymptotic_kernel_in_interpolation': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb], 'kernel': M.delta_expansion})),
        'gamma_only_long_range_model': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb]}, gamma_only=True)),
        'q_fraction_of_2pi_over_L': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb]}, qscale=2.0)),
        'e2_hartree_in_long_range': table(lambda tb, s, q: predict(s, q, lr={'t': t_or[tb], 'scale': 0.5})),
    }
    # 3D boundary conditions: the archived 3D-periodic interpolations of the hidden tubes
    abl['three_d_boundary_condition_dispersion'] = {}
    for tb in TARGETS:
        qf, lo = archived_lo_curve(HIDDEN_TUBES[tb] + '/3D_code', sys3[tb])
        abl['three_d_boundary_condition_dispersion'][tb] = {k: float(np.interp(float(k), qf, lo)) for k in TARGETS[tb]}
    # 2D-like linear collapse fitted to the purchasable 0.15 run of each hidden tube
    abl['linear_2d_collapse'] = {}
    for tb in TARGETS:
        B = (purch015[tb] ** 2 - gam[tb] ** 2) / 0.15
        abl['linear_2d_collapse'][tb] = {k: float(np.sqrt(gam[tb] ** 2 + B * float(k))) for k in TARGETS[tb]}
    # family transfer: the released (4,4) shift curve added to the hidden tube's Gamma LO
    qs44 = [0.0] + sorted(pts44_all)
    sh44 = [0.0] + [pts44_all[q] - G44 for q in sorted(pts44_all)]
    abl['family_transfer_from_44'] = {tb: {k: float(gam[tb] + np.interp(float(k), qs44, sh44)) for k in TARGETS[tb]} for tb in TARGETS}
    # ---- 3D q->0 LO of BN(5,5) ----------------------------------------------------------------------------------
    s3 = sys3['BN55']
    s1 = sysd['BN55']
    t4 = s1['alat'] / 4
    eps1d = s1['area'] * (s1['eps'][2, 2] - 1) / (np.pi * t4 ** 2)
    v3 = {'oracle_3D_code_gamma': M.lo(s3, s3['Dq'][0] + M.na_3d(s3)),
          'robust_1D_code_gamma': M.lo(s1, s1['Dq'][0] + M.na_3d(s1)),
          'abl_cylinder_dielectric': M.lo(s3, s3['Dq'][0] + M.na_3d(s3, eps=np.eye(3) * eps1d)),
          'abl_cylinder_dielectric_and_volume': M.lo(s3, s3['Dq'][0] + M.na_3d(s3, eps=np.eye(3) * eps1d, omega=np.pi * t4 ** 2 * s3['c'])),
          'abl_perpendicular_dielectric': M.lo(s3, s3['Dq'][0] + M.na_3d(s3, eps=np.eye(3) * s3['eps'][0, 0])),
          'abl_unscreened': M.lo(s3, s3['Dq'][0] + M.na_3d(s3, eps=np.eye(3))),
          'abl_analytic_gamma_frequency': M.lo(s3, s3['Dq'][0])}
    # ---- write private data and ablation values ------------------------------------------------------------------
    thresholds = {'tolerance_cm-1': TOL, 'tolerance_3d_cm-1': TOL_3D}
    hidden = {'targets_1D_LO_cm-1': {tb: {k: round(v, 4) for k, v in truth[tb].items()} for tb in truth},
              'target_3D_LO_q0_cm-1': round(truth_3d, 4), 'thresholds': thresholds,
              'definition': 'explicit archived 1D-boundary-condition DFPT runs (longwave folders); LO = mode with the largest '
                            'longitudinal mode effective charge, Born charges of the Gamma run; 3D: q->0 point of the archived '
                            '3D-periodic matdyn.modes of BN(5,5)'}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [round(v, 4) for tb in truth for v in truth[tb].values()] + [round(truth_3d, 4)]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')

    def r4(d):
        return {tb: {k: round(float(v), 4) for k, v in d[tb].items()} for tb in d}
    ablation_values = {'ablations': {k: r4(v) for k, v in abl.items()},
                       'ablations_3d': {k: round(v, 4) for k, v in v3.items() if k.startswith('abl_')},
                       'robustness': {k: r4(v) for k, v in robust.items()},
                       'robustness_3d': {'1D_code_gamma': round(v3['robust_1D_code_gamma'], 4)}}
    (TASK/'solution/ablation_values.json').write_text(json.dumps(ablation_values, indent=1) + '\n')
    # ---- catalog ------------------------------------------------------------------------------------------------
    catalog = {
        'archive': {'citation': 'N. Rivano, N. Marzari, T. Sohier, npj Comput. Mater. 9, 194 (2023), DOI 10.1038/s41524-023-01140-2',
                    'record': 'Materials Cloud Archive 08gsg-gpd20, DOI 10.24435/materialscloud:46-wj (CC BY 4.0)',
                    'README.txt': (ROOT/'README.txt').read_text(),
                    'files_description.md': (ROOT.parent/'files_description.md').read_text() if (ROOT.parent/'files_description.md').exists() else ''},
        'systems': SYSTEMS,
        'folders': {f: {'boundary_conditions': bc_of(f), 'input_files': inputs[f]} for f in inputs},
        'conventions': {
            'dyn files': 'Quantum ESPRESSO ph.x fildyn format, as archived: q in cartesian units of 2pi/alat with alat = celldm(1) '
                         '(the lateral cell size for ibrav 6, given in the header for ibrav 0); force-constant matrices in Ry/bohr^2 '
                         'per atom pair (not mass-scaled); masses in Rydberg units (amu x 911.444); Gamma files of the 1D systems also '
                         'carry the dielectric tensor of the simulation cell and the Born effective charges (E-U: [field][displacement]).',
            'boundary conditions': "'1D' = the authors' modified Quantum ESPRESSO with a 1D Coulomb cutoff (folders 1D_code); "
                                   "'3D' = standard 3D-periodic Quantum ESPRESSO on the same supercell (folders 3D_code).",
            'interpolated': 'archived q2r/matdyn interpolations along Gamma-Z (matdyn.in shows the settings).'},
        'requests': {
            'dfpt': {'request': {'kind': 'dfpt', 'file': '<id from dfpt_runs>'}, 'cost': COSTS['dfpt'],
                     'returns': 'the archived dynamical-matrix file as text'},
            'interpolated': {'request': {'kind': 'interpolated', 'folder': '<folder from interpolated_dispersions>'},
                             'cost': COSTS['interpolated'],
                             'returns': 'frequencies along Gamma-Z and the axial share of each mode, parsed from the archived matdyn.modes'},
            'archived_file': {'request': {'kind': 'archived_file', 'file': '<path from archived_outputs>'},
                              'cost': COSTS['archived_file'],
                              'returns': f'an archived pw.x/ph.x/q2r.x/matdyn.x output as text (first {MAX_TEXT} characters)'}},
        'budget': BUDGET,
        'dfpt_runs': items,
        'interpolated_dispersions': interp_items,
        'archived_outputs': file_items,
        'withheld': ['BN(5,5) and BN(6,6), 1D boundary conditions: the archived explicit small-q runs below 0.15 |Gamma-Z| '
                     '(longwave folders and their copies)',
                     'BN(5,5) and BN(6,6), both boundary conditions: the archived interpolated dispersions (matdyn.modes) and matdyn.out'],
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    for f in ('environment/lab/public.json', 'tests/public.json'):
        p = json.loads((TASK/f).read_text())
        p['solver_budget'] = BUDGET
        if f.startswith('environment'):
            p['requests'] = {k: {'cost': v} for k, v in COSTS.items()}
            p['catalog'] = '/app/catalog.json'
        (TASK/f).write_text(json.dumps(p, indent=1) + '\n')
    # ---- calibration ---------------------------------------------------------------------------------------------
    def err(pred):
        return {tb: {k: round(pred[tb][k] - truth[tb][k], 2) for k in truth[tb]} for tb in truth}
    calib = {'kernel_checks_max_dev_cm-1': checks,
             'truth': {tb: {k: round(v, 3) for k, v in truth[tb].items()} for tb in truth},
             'truth_3d': round(truth_3d, 3),
             'purchasable_0.15_runs': {k: round(v, 3) for k, v in purch015.items()},
             'gamma_LO': {k: round(v, 3) for k, v in gam.items()},
             'oracle_protocol': {'t_over_R_from_44': round(ratio, 4), 't44_bohr': round(best44[0], 3), 'rms44': round(best44[1], 3),
                                 't_bohr': {k: round(v, 3) for k, v in t_or.items()}, 'R_atomic_bohr': {k: round(v, 3) for k, v in R.items()}},
             'oracle_errors': err(oracle),
             'robustness_errors': {k: err(v) for k, v in robust.items()},
             'ablation_errors': {k: err(v) for k, v in abl.items()},
             'three_d': {k: round(v, 3) for k, v in v3.items()},
             'tolerances': thresholds,
             'catalog_units': {'dfpt': len(items), 'interpolated': 2 * len(interp_items), 'archived_file': len(file_items)},
             'budget': BUDGET}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    print(json.dumps(calib, indent=1))


if __name__ == '__main__':
    main()
