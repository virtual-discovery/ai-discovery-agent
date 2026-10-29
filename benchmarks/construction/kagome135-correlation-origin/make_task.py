"""Write the staged files, lab data, catalog, private targets and ablation values (run in the runtime image).

    bash run_in_image.sh make_task.py

Inputs: /data/x (fetch.py), /data/models (make_model.py), /data/refs_*.json (make_refs.py, lab code path).
Writes into benchmarks/tasks/kagome135-correlation-origin and calibration.json here.
"""
import json
import shutil
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '.')
import kag_ssmf as K

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/kagome135-correlation-origin'
X = Path('/data/x')
FOLDERS = {'TiSb': ('CsTi3Sb5', 'CsTiSb'), 'VSb': ('CsV3Sb5', 'CsVSb'), 'CrSb': ('CsCr3Sb5', 'CsCrSb'),
           'TiBi': ('CsTi3Bi5', 'CsTiBi'), 'CrBi': ('CsCr3Bi5', 'CsCrBi'), 'CrBi_9.2': ('CsCr3Bi5_9.2', 'CsCrBi')}
ORB = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']                 # archive (.win) order, site-local frames
W90_DEFAULT = ['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']         # Wannier90 order of a plain "d" projection


def refs():
    out = {}
    for f in sorted(Path('/data').glob('refs_*.json')):
        out.update(json.loads(f.read_text()))
    return out


def zbar(r):
    return float(np.mean(r['Z']))


def stage():
    agent, lab = TASK/'environment/agent', TASK/'environment/lab'
    for d in (agent/'archive', agent/'models', lab/'data'):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    for folder, (key, stem) in FOLDERS.items():
        dst = agent/'archive'/folder
        dst.mkdir()
        for name in ['Fermi_level.dat', f'turn_{stem}.win', f'turn_{stem}.wout', f'turn_{stem}_bands.dat.gnu']:
            shutil.copy(X/folder/name, dst/name)
        m = np.load(f'/data/models/{key}.npz')
        arrays = {'R': m['R'].astype(np.int8), 'degen': m['degen'].astype(np.int16), 'H': m['H'].astype(np.complex128)}
        np.savez_compressed(agent/'models'/f'{folder}.npz', **arrays)
        np.savez_compressed(lab/'data'/f'{folder}.npz', **arrays, efermi=m['efermi'], cell=m['cell'])
    (agent/'archive/DFT').mkdir()
    shutil.copy(X/'DFT/scf.in', agent/'archive/DFT/scf.in')
    # scf.out without the per-k eigenvalue listing (between 'End of self-consistent calculation' and the
    # Fermi-energy line), which is redundant with the models and collides numerically with graded values.
    lines = (X/'DFT/scf.out').read_text().splitlines(keepends=True)
    a = next(i for i, l in enumerate(lines) if 'End of self-consistent calculation' in l)
    b = next(i for i, l in enumerate(lines) if 'the Fermi energy is' in l)
    (agent/'archive/DFT/scf.out').write_text(''.join(lines[:a + 1]) + '\n     [per-k eigenvalue listing removed]\n\n'
                                              + ''.join(lines[b:]))
    shutil.copy(HERE/'kag_ssmf.py', lab/'kag_ssmf.py')


def sigma_variants(model, mesh=(18, 18, 12)):
    """Windowed (|E - EF| < 4 eV) PDOS moment: a plausible wrong estimator of the kinetic scale."""
    ks, w = K.kmesh_tr(*mesh)
    e, v = np.linalg.eigh(K.hk(model, ks))
    wgt = (np.abs(v[:, :5, :]) ** 2) * w[:, None, None]
    E = e[:, None, :] - model['efermi']
    win = np.abs(E) < 4
    m1 = (wgt * E * win).sum((0, 2)) / (wgt * win).sum((0, 2))
    m2 = (wgt * (E - m1[None, :, None]) ** 2 * win).sum((0, 2)) / (wgt * win).sum((0, 2))
    H = model['H'].copy()
    H[:, :15, 15:] = 0
    H[:, 15:, :15] = 0
    return np.sqrt(m2), K.sigma_exact(dict(model, H=H))[:5]


CATALOG = {
    'material_family': 'CsM3Pn5 kagome metals (M = Ti, V, Cr; Pn = Sb, Bi), nonmagnetic PBE (Quantum ESPRESSO 7.1, no spin-orbit coupling), Wannier90 models of 33 orbitals per cell',
    'source': 'Zenodo record 17951035 (DOI 10.5281/zenodo.17951035, CC BY 4.0), S. Enzner and G. Sangiovanni, "135 Kagome Wannier Hamiltonians"',
    'staged_files': {
        'archive': '/app/archive/<model>/: Fermi_level.dat, turn_*.win, turn_*.wout, turn_*_bands.dat.gnu (DFT bands), unchanged from the record; /app/archive/DFT/: the record\'s example QE input scf.in and output scf.out (per-k eigenvalue listing removed; pseudopotentials, rel.in/rel.out and the Wannier band files turn_*_band.dat not staged)',
        'models': '/app/models/<model>.npz: the archived turn_*_hr.dat converted losslessly to numpy: R (n_R, 3) lattice vectors in units of the cell of the .wout, degen (n_R,) Wigner-Seitz degeneracies, H (n_R, 33, 33) complex hoppings H_mn(R) = <m,0|H|n,R> in eV on the absolute energy scale of the archive (H(k) = sum_R exp(2 pi i k.R) H(R) / degen(R)).',
    },
    'models': {'TiSb': 'archive folder TiSb (CsTiSb.zip)', 'VSb': 'archive folder VSb (CsVSb.zip)',
               'CrSb': 'archive folder CrSb (CsCrSb.zip)', 'TiBi': 'archive folder TiBi (CsTiBi.zip)',
               'CrBi': 'archive folder CrBi (CsCrBi.zip)', 'CrBi_9.2': 'archive folder CrBi_9.2 (CsCrBi_9.2.zip)'},
    'requests': {'ssmf': {'format': {'kind': 'ssmf', 'model': 'TiSb | VSb | CrSb | TiBi | CrBi | CrBi_9.2',
                                     'electrons': 'electrons per cell in the 33-orbital model, both spins (26-36)',
                                     'U_eV': '0-8', 'J_over_U': '0-0.3'},
                          'cost': 1,
                          'returns': 'metallic_solution; if true, Z and occupation per spin for Wannier functions 1-15 (the d orbitals) and the chemical potential (eV, absolute scale of the model)'}},
    'budget': 20,
    'slave_spin_method': {
        'correlated_orbitals': 'Wannier functions 1-15 (5 d orbitals on each of the 3 transition-metal sites); all other orbitals are uncorrelated',
        'interaction': 'density-density Kanamori: U n_up n_dn per orbital, U - 2J between opposite spins in different orbitals, U - 3J between equal spins in different orbitals',
        'double_counting': 'fully localized limit',
        'solution': 'Z2 slave-spin mean field (paramagnetic, T = 0, 10 meV Fermi smearing, 9x9x6 k mesh); the metallic solution is followed continuously from weak coupling (U ramped at fixed J/U); the three transition-metal sites are symmetry-equivalent',
        'Z': 'quasiparticle weight of each d Wannier function (1 = uncorrelated, 0 = Mott localized)',
    },
}


def sigma_mesh(model, mesh=(6, 6, 4)):
    """Full-range PDOS second moment from eigenvectors on a coarse uniform mesh (an alternative estimator)."""
    ks, w = K.kmesh_tr(*mesh)
    e, v = np.linalg.eigh(K.hk(model, ks))
    wgt = (np.abs(v[:, :5, :]) ** 2) * w[:, None, None]
    m1 = (wgt * e[:, None, :]).sum((0, 2)) / wgt.sum((0, 2))
    return np.sqrt((wgt * (e[:, None, :] - m1[None, :, None]) ** 2).sum((0, 2)) / wgt.sum((0, 2)))


def main():
    stage()
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(CATALOG, indent=1) + '\n')
    R = refs()
    get = lambda job: R[job]
    models = {f: K.load_model(f'/data/models/{k}.npz') for f, (k, s) in FOLDERS.items()}
    counts = {}
    for f, m in models.items():
        ks, w = K.kmesh_tr(24, 24, 16)
        e = np.linalg.eigvalsh(K.hk(m, ks))
        counts[f] = K.electron_count(e, m['efermi'], 0.002, w)
    pristine = {f: int(round(c)) for f, c in counts.items()}
    nominal = {'TiSb': 27, 'VSb': 30, 'CrSb': 33, 'TiBi': 27, 'CrBi': 33, 'CrBi_9.2': 33}
    J = '0.18'
    cr = get(f'CrSb:34:6:{J}')
    imin = int(np.argmin(cr['Z']))
    sig = {f: K.sigma_exact(m)[:5] for f, m in models.items()}

    def share(nA, nB, A='TiSb', B='CrSb', suffix=''):
        L = lambda mdl, n: np.log(zbar(get(f'{mdl}:{n}:6:{J}{suffix}')))
        total = L(B, nB) - L(A, nA)
        fill = 0.5 * ((L(A, nB) - L(A, nA)) + (L(B, nB) - L(B, nA)))
        one_A = (L(A, nB) - L(A, nA)) / total
        one_B = (L(B, nB) - L(B, nA)) / total
        return fill / total, one_A, one_B

    ref = {
        'electron_counts': pristine,
        'CrSb_Z_mean': zbar(cr), 'CrSb_Z_min': float(min(cr['Z'])),
        'most_renormalized_orbital': ORB[imin],
        'filling_share': share(28, 34)[0],
        'hund_ratio': zbar(get('CrSb:34:6:0')) / zbar(cr),
        'CrBi_U4p5_Z_mean': zbar(get(f'CrBi:34:4.5:{J}')),
        'CrBi92_U4p5_Z_mean': zbar(get(f'CrBi_9.2:34:4.5:{J}')),
        'sigma_CrBi': float(sig['CrBi'][imin]), 'sigma_TiSb': float(sig['TiSb'][imin]),
    }
    tol = {'Z_abs': 0.025, 'share_abs': 0.08, 'hund_rel': 0.10, 'sigma_rel': 0.03}
    # ---------------------------------------------------------------- wrong paths (ablations)
    crn = get(f'CrSb:33:6:{J}')
    sh_nom = share(27, 33)
    sh_one = share(28, 34)
    # cross-compound regression ln Zbar = a + b nbar + c sigma_bar over the pristine Sb/Bi compounds
    rows = [(f, n) for f, n in [('TiSb', 28), ('VSb', 31), ('CrSb', 34), ('TiBi', 28)]]
    Xr = np.array([[1.0, np.mean(get(f'{f}:{n}:6:{J}')['n']), np.mean(sig[f])] for f, n in rows])
    yr = np.array([np.log(zbar(get(f'{f}:{n}:6:{J}'))) for f, n in rows])
    coef = np.linalg.lstsq(Xr, yr, rcond=None)[0]
    reg_share = coef[1] * (Xr[2, 1] - Xr[0, 1]) / (yr[2] - yr[0])
    win_CrBi, donly_CrBi = sigma_variants(models['CrBi'])
    win_TiSb, donly_TiSb = sigma_variants(models['TiSb'])
    ablations = {
        'nominal_counts': nominal,
        'CrSb_nominal': {'Z_mean': zbar(crn), 'Z_min': float(min(crn['Z']))},
        'share_nominal': sh_nom[0],
        'hund_nominal': zbar(get('CrSb:33:6:0')) / zbar(crn),
        'CrBi_nominal': {'Sb_lattice_Z_mean': zbar(get(f'CrBi:33:4.5:{J}')),
                         'TiBi_lattice_Z_mean': zbar(get(f'CrBi_9.2:33:4.5:{J}'))},
        'CrBi_swapped': {'Sb_lattice_Z_mean': ref['CrBi92_U4p5_Z_mean'], 'TiBi_lattice_Z_mean': ref['CrBi_U4p5_Z_mean']},
        'sigma_swapped_CrBi': {'CsCr3Bi5_Sb_lattice': float(sig['CrBi_9.2'][imin]), 'CsTi3Sb5': ref['sigma_TiSb']},
        'orbital_w90_default_order': W90_DEFAULT[imin],
        'orbital_paper_frame': 'dyz',
        'share_one_sided_Ti': sh_one[1], 'share_one_sided_Cr': sh_one[2], 'share_regression': float(reg_share),
        'sigma_windowed': {'CsCr3Bi5_Sb_lattice': float(win_CrBi[imin]), 'CsTi3Sb5': float(win_TiSb[imin])},
        'sigma_d_only': {'CsCr3Bi5_Sb_lattice': float(donly_CrBi[imin]), 'CsTi3Sb5': float(donly_TiSb[imin])},
        'sigma_mean_orbital': {'CsCr3Bi5_Sb_lattice': float(np.mean(sig['CrBi'])), 'CsTi3Sb5': float(np.mean(sig['TiSb']))},
    }
    rob = {}
    m12 = lambda job: R.get(job + ':12,12,8')
    if m12(f'CrSb:34:6:{J}'):
        c12 = m12(f'CrSb:34:6:{J}')
        rob['CrSb_mesh12'] = {'Z_mean': zbar(c12), 'Z_min': float(min(c12['Z']))}
        rob['hund_mesh12'] = zbar(m12('CrSb:34:6:0')) / zbar(c12)
        rob['CrBi_mesh12'] = {'Sb_lattice_Z_mean': zbar(m12(f'CrBi:34:4.5:{J}')),
                              'TiBi_lattice_Z_mean': zbar(m12(f'CrBi_9.2:34:4.5:{J}'))}
        rob['share_mesh12'] = share(28, 34, suffix=':12,12,8')[0]
    rob['sigma_mesh_histogram'] = {'CsCr3Bi5_Sb_lattice': float(sigma_mesh(models['CrBi'])[imin]),
                                   'CsTi3Sb5': float(sigma_mesh(models['TiSb'])[imin])}
    (TASK/'environment/lab/private').mkdir(parents=True, exist_ok=True)
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps({'reference': ref, 'tolerances': tol}, indent=1) + '\n')
    hidden = [ref['CrSb_Z_mean'], ref['CrSb_Z_min'], ref['filling_share'], ref['hund_ratio'], ref['CrBi_U4p5_Z_mean'],
              ref['CrBi92_U4p5_Z_mean'], ref['sigma_CrBi'], ref['sigma_TiSb']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(hidden) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': ablations, 'robustness': rob}, indent=1) + '\n')
    calib = {'electron_count_at_archived_EF_24x24x16': counts, 'reference': ref, 'tolerances': tol,
             'ablations': ablations, 'robustness': rob, 'regression_coefficients': coef.tolist(),
             'share_parts': {'one_sided_Ti_path': sh_one[1], 'one_sided_Cr_path': sh_one[2]},
             'raw_refs': R}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    print(json.dumps({'reference': ref, 'ablations': ablations, 'robustness': rob}, indent=1))


if __name__ == '__main__':
    main()
