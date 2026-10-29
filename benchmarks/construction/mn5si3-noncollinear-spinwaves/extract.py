"""Extract everything the task uses from the two unpacked Materials Cloud archives into CACHE/pool.json.

    python3 -B extract.py CACHE      # after fetch.py; CACHE = ~/.cache/ai-theorist/mn5si3-noncollinear-spinwaves

Sources (CC BY 4.0): paxx3-a1k09 (PRB 105, 104404 (2022)) -> CACHE/x22, 4tgt6-1p378 (APL Mater. 11, 081103 (2023)) -> CACHE/x24.
"""
import glob
import json
import os
import sys
from pathlib import Path
import numpy as np

CACHE = Path(sys.argv[1])
X22 = CACHE/'x22/Spin-wave calculations'
X24 = CACHE/'x24'
F2 = X24/'Fig 2'
A = 13.0564948 * 0.529177210903          # DFT lattice constant (bohr -> A), inputcard ALATBASIS
CA = 0.69667149059                        # c/a of the spin-wave calculations (inputcard_MnSiDisper.inp a3)
L_A = np.array([1.0, np.sqrt(3), CA])     # orthorhombic cell in units of a


def rows(path, ncols=None, skip=1):
    out = []
    for line in Path(path).read_text().replace('\r', '\n').splitlines()[skip:]:
        p = line.split()
        if p and (ncols is None or len(p) == ncols):
            out.append([float(v) for v in p])
    return np.array(out)


def tas(path):
    """Raw triple-axis file (CR line ends, tab header)."""
    lines = [l for l in Path(path).read_text().replace('\r\n', '\n').replace('\r', '\n').split('\n') if l.strip()]
    head = lines[0].split('\t')
    data = np.array([[float(v) for v in l.split()] for l in lines[1:]])
    return {h: data[:, k].tolist() for k, h in enumerate(head) if k < data.shape[1]}


def main():
    pos = np.loadtxt(F2/'inputfiles/lattice_MnSi2_noncol.dat')          # 16 magnetic-model sites, fractional
    cart = pos * L_A
    # ---------------------------------------------------------------- full structure from the DFT inputcard
    card = (CACHE/'x22/DFT/inputcard').read_text().splitlines()
    k0 = [k for k, l in enumerate(card) if l.strip() == 'RBASIS'][0]
    sites = []
    for l in card[k0 + 1:k0 + 33]:
        p = l.split()
        f = np.array([float(v) for v in p[:3]]) / np.array([1, np.sqrt(3), 0.69679557694152])
        f = f - np.floor(f + 1e-9)
        sites.append((p[3], f))
    pm = pos - np.floor(pos + 1e-9)
    struct, n2b, nsi = [], 0, 0
    model_label = {}
    for elem, f in sites:
        d = np.abs(((pm - f) + 0.5) % 1 - 0.5).max(1)
        j = int(np.argmin(d))
        if d[j] < 1e-3:
            lab = ('Mn2a_%d' % (j + 1)) if j < 8 else ('Mn1_%d' % (j - 7))
            model_label[j] = lab
        elif elem == 'Mn2':
            n2b += 1; lab = 'Mn2b_%d' % n2b
        else:
            nsi += 1; lab = 'Si_%d' % nsi
        struct.append({'label': lab, 'element': 'Mn' if elem.startswith('Mn') else 'Si', 'wyckoff_type': elem,
                       'frac': (pos[j] - np.floor(pos[j] + 1e-9)).round(5).tolist() if d[j] < 1e-3 else f.round(5).tolist()})
    assert len(model_label) == 16
    struct = {'sites': struct, 'model_index_label': [model_label[j] for j in range(16)]}
    # ---------------------------------------------------------------- exchange tables
    def from_cartesian(fname):
        out = []
        for r in rows(F2/'inputfiles'/fname, 13):
            i, j = int(r[0]), int(r[1])
            d = r[9:12]
            diff = (cart[i] + d - cart[j]) / L_A
            assert np.abs(diff - np.round(diff)).max() < 1e-3
            out.append({'i': i, 'j': j, 'vec_A': (d * A).round(4).tolist(), 'J_meV': round(float(r[8]), 5)})
        return out

    def from_offsets(fname, ncol):
        out = []
        for r in rows(F2/'inputfiles'/fname):
            i, j = int(r[0]), int(r[1])
            d = (np.array(r[2:5]) + pos[j] - pos[i]) * L_A
            out.append({'i': i, 'j': j, 'vec_A': (d * A).round(4).tolist(), 'J_meV': round(float(r[8]), 5)})
        return out
    temp = rows(F2/'inputfiles/pair_MnSi_rad4.5noncol_temp.txt')
    mn1 = []
    for r in temp:
        i, j = int(r[0]), int(r[1])
        if i >= 8 or j >= 8:
            d = (np.array(r[2:5]) + pos[j] - pos[i]) * L_A
            mn1.append({'i': i, 'j': j, 'vec_A': (d * A).round(4).tolist(), 'J_meV': round(float(r[8]), 5)})
    exch = {'noncollinear_r4.5': from_cartesian('pair_MnSi_rad4.5noncol.txt'),
            'collinear_AFMz_r4.5': from_offsets('pairMnSi_AFMz.txt', 9),
            'noncollinear_r17': from_cartesian('pair_MnSi_rad17.txt'),
            'mn1_couplings': mn1}
    for k, v in exch.items():
        d = sorted({round(float(np.linalg.norm(b['vec_A'])), 2) for b in v})
        print(k, len(v), d[:8])
    # archived (buggy-offset) table used by the authors' SWIS run, kept for the build notes
    archived_offsets = rows(F2/'inputfiles/pair_MnSi_rad4.5noncol_temp.txt')[:, [0, 1, 2, 3, 4, 8]].tolist()
    # ---------------------------------------------------------------- neutron data
    disp = {}
    for key, name in (('h', 'Q=(Qh,1,0)'), ('k', 'Q=(0,Qk,1)'), ('l', 'Q=(0,3,Ql)')):
        d = np.loadtxt(F2/'INS_data_Mn5Si3'/(name + '_10K_dispersion.txt'), comments='#')
        disp[key] = d[:, :3].tolist() if d.shape[1] == 3 else d.tolist()
    cuts = {}
    for sub, key in (('Fig4a', 'h'), ('Fig4b', 'k'), ('Fig4c', 'l'), ('Fig4d', 'l')):
        for f in sorted(glob.glob(str(X22/'exp data/Fig4'/sub/'E=*meV.txt'))):
            E = float(os.path.basename(f)[2:-7])
            a = np.loadtxt(f, comments='#')
            cuts['%s_%g' % (key, E)] = {'direction': key, 'energy_meV': E, 'x': a[:, 0].round(4).tolist(),
                                        'counts': a[:, 1].tolist(), 'error': a[:, 2].round(3).tolist(), 'figure': sub}
    f3a = np.loadtxt(X22/'exp data/Fig3/data/Fig3a_data.txt', comments='#')
    escans = {}
    for k, dq in enumerate([0, 0.01, 0.02, 0.03, 0.04, 0.05]):
        escans['gap_10K_dh%.2f' % dq] = {'Q_rlu': [1 + dq, 2, 0], 'T_K': 10, 'B_T': 0, 'E_meV': f3a[:, 0].round(4).tolist(),
                                          'counts': f3a[:, 1 + 2 * k].tolist(), 'error': f3a[:, 2 + 2 * k].round(3).tolist()}
    f3c = np.loadtxt(X22/'exp data/Fig3/data/Fig3c_data.txt', comments='#')
    escans['optical_210_10K'] = {'Q_rlu': [2, 1, 0], 'T_K': 10, 'B_T': 0, 'E_meV': f3c[:, 0].round(4).tolist(),
                                 'counts': f3c[:, 1].tolist(), 'error': f3c[:, 2].round(3).tolist()}
    for T in (2, 10, 20, 30, 40, 50, 60):
        r = tas(X24/'Fig 3a'/('E-scan_%dK.txt' % T))
        escans['030_%dK' % T] = {'Q_rlu': [0, 3, 0], 'T_K': T, 'B_T': 0, 'E_meV': r['EN'], 'counts': r['CNTS'],
                                 'error': r['CNTS_er'], 'monitor': r['M1'][0], 'T_measured_K': r['TT']}
    for B in (0, 2, 4, 6, 8, 10):
        r = tas(X24/'Fig 3c'/('E-scan_%dT.txt' % B))
        escans['120_50K_%dT' % B] = {'Q_rlu': [1, 2, 0], 'T_K': 50, 'B_T': B, 'E_meV': r['EN'], 'counts': r['CNTS'],
                                     'error': r['CNTS_er'], 'monitor': r['M1'][0]}
    gap = np.loadtxt(X22/'exp data/Fig3/data/Fig3b_data.txt', comments='#')
    gap_table = {'10K': {'q_h': gap[:, 0].tolist(), 'E_meV': gap[:, 1].tolist(), 'E_err': gap[:, 2].tolist()},
                 '80K': {'q_h': gap[:, 3].tolist(), 'E_alpha_meV': gap[:, 4].tolist(), 'E_alpha_err': gap[:, 5].tolist(),
                         'E_beta_meV': gap[:, 6].tolist(), 'E_beta_err': gap[:, 7].tolist()}}
    qscans = {}
    for f in sorted(glob.glob(str(X24/'Fig 5b/*.txt'))) + sorted(glob.glob(str(X24/'Fig 5c/*.txt'))) + \
            sorted(glob.glob(str(X24/'Fig 4c/*.txt'))) + sorted(glob.glob(str(X24/'Fig 4a/*.txt'))) + \
            [str(X24/'Fig 4b/2meV_4T_50K.txt'), str(X24/'Fig 4b/2meV_10T_50K.txt')]:
        r = tas(f)
        E = round(float(np.mean(r['EN'])), 2)
        T = int(round(float(np.mean(r['TT'])) / 10) * 10) if np.mean(r['TT']) > 5 else 2
        B = int(round(float(np.mean(r['MAG']))))
        key = 'Qh20_%gmeV_%dK_%dT' % (E, T, abs(B))
        if key in qscans:
            continue
        qscans[key] = {'Q': '(Qh,2,0)', 'E_meV': E, 'T_K': T, 'B_T': abs(B), 'Qh': r['QH'], 'counts': r['CNTS'],
                       'error': r['CNTS_er'], 'T_measured_K': r['TT']}
    phase = rows(X24/'Fig 1c/data/H-T phase diagram_new.txt')
    # ---------------------------------------------------------------- archived ground state (graded canting)
    gs = np.loadtxt(F2/'spirit_output/spinconfig_MnSinoncol_final_h=0.00.ovf', comments='#')
    gs = gs / np.linalg.norm(gs, axis=1, keepdims=True)
    n1 = gs[8]
    canting = float(np.degrees(np.arcsin(np.abs(gs[:8] @ n1).mean())))
    print('archived ground-state canting (deg):', canting, 'Mn1 spread', np.abs(gs[8:] @ n1).min())
    pool = {'structure': struct, 'model_positions': pos.tolist(), 'a_A': A, 'c_over_a': CA, 'exchange': exch,
            'archived_offset_table': archived_offsets, 'dispersion': disp, 'cuts': cuts, 'escans': escans,
            'gap_table': gap_table, 'qscans': qscans,
            'phase_diagram': {'columns': ['T1', 'B_AFM1prime_AFM2', 'T2', 'B_AFM1pp_AFM1prime', 'T3', 'B_AFM2_PM',
                                          'T4', 'B_AFM1_AFM1pp'], 'rows': phase.tolist()},
            'archived_ground_state': gs.tolist(), 'archived_canting_deg': canting}
    (CACHE/'pool.json').write_text(json.dumps(pool))
    print('cuts', sorted(cuts), '\nescans', sorted(escans), '\nqscans', len(qscans), sorted(qscans)[:5])


if __name__ == '__main__':
    main()
