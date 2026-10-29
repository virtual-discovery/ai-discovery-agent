"""Extract every number used by cnt-kohn-anomaly-nonadiabatic from the unpacked archive (no modification).

    python3 -B extract.py [CACHE]     (default ~/.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic)

Writes CACHE/pool/:
  coupling.npz   archived g+q / g-q tables (k grid, energy-ordered bands at k and k+-Q, |g|, omega_Q) for nu = 1, 20, 23
  pit.npz        all 37 archived self-energy tables (Fig2/PiT/<T>): w_q, w grid, the three value columns
  pit_text.json  the verbatim header and table lines of nu = 1, 20, 23 for every temperature
  archive.json   bands (band.eig), phonons (matdyn.freq), the DFPT smearing scan of fig2_d_f_h.py,
                 constants and comments quoted from the Fig. 2 scripts, and the text of the small archive files
"""
import ast
import json
import re
import sys
from pathlib import Path
import numpy as np

CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/cnt-kohn-anomaly-nonadiabatic'
RAW, POOL = CACHE/'raw', CACHE/'pool'
MODES = (1, 20, 23)
TEXT_FILES = ['README.md', 'dft_inputs/scf.in', 'dft_inputs/ph.in', 'dft_inputs/epw.in', 'dft_inputs/specfun/epw.in',
              'dft_inputs/specfun/input.py', 'dft_inputs/specfun/pathq.dat', 'Fig1/coupling_calculator.c',
              'Fig1/fig1.py', 'Fig1/backup/output.csv']


def read_bands(path):
    lines = path.read_text().split('\n')
    nb = int(re.search(r'nbnd=\s*(\d+)', lines[0]).group(1))
    nk = int(re.search(r'nks=\s*(\d+)', lines[0]).group(1))
    vals = np.array(' '.join(lines[1:]).split(), float).reshape(nk, 3 + nb)
    return vals[:, :3], vals[:, 3:]


def main():
    POOL.mkdir(parents=True, exist_ok=True)
    # coupling tables: columns n, m, nu, e_k (eV), e_k+-q (eV), omega_q (meV), |g| (meV); k = -0.5 ... 0.5 (2000 points)
    arrays = {'k': np.linspace(-0.5, 0.5, 2000, endpoint=False)}
    for nu in MODES:
        for tag, proc in (('p', '+'), ('m', '-')):
            e, eq, g = np.zeros((2, 2000)), np.zeros((2, 2000)), np.zeros((2, 2, 2000))
            for n in (1, 2):
                for m in (1, 2):
                    d = np.loadtxt(RAW/f'Fig2/g/g{proc}q.{n}{m}{nu}')
                    assert d.shape == (2000, 7) and (d[:, 0] == n).all() and (d[:, 1] == m).all() and (d[:, 2] == nu).all()
                    if n == 1 and m == 1:
                        wq = d[0, 5]
                    assert np.allclose(d[:, 5], wq)
                    e[n - 1], eq[m - 1], g[n - 1, m - 1] = d[:, 3], d[:, 4], d[:, 6]
            arrays[f'e_{nu}_{tag}'], arrays[f'eq_{nu}_{tag}'], arrays[f'g_{nu}_{tag}'] = e, eq, g
            arrays[f'wq_{nu}_{tag}'] = np.array(wq)
    np.savez_compressed(POOL/'coupling.npz', **arrays)
    # self-energy tables: 8 numeric columns (q index, mode, T, w_q eV, w eV, c6, c7, c8), 36 modes per w
    temps = sorted(float(p.name) for p in (RAW/'Fig2/PiT').iterdir())
    name = {t: ('1578.87' if abs(t - 1578.87) < 1e-6 else str(int(t))) for t in temps}
    wq, c6, c7, c8, text = [], [], [], [], {}
    for t in temps:
        lines = (RAW/'Fig2/PiT'/name[t]).read_text().split('\n')
        header, body = lines[:2], [l for l in lines[2:] if l.strip()]
        d = np.array([l.split() for l in body], float)
        assert d.shape == (36 * 801, 8) and np.allclose(d[:, 2], t)
        d = d.reshape(801, 36, 8)
        w = d[:, 0, 4]
        assert np.allclose(d[:, :, 4], w[:, None]) and np.allclose(d[:, :, 3], d[0, :, 3][None, :])
        assert np.ptp(d[:, :, 6], axis=0).max() == 0.0          # column 7 (1-based) does not depend on w
        wq.append(d[0, :, 3]); c6.append(d[:, :, 5].T); c7.append(d[:, :, 6].T); c8.append(d[:, :, 7].T)
        text[name[t]] = {'header': header, **{str(nu): [body[i * 36 + nu - 1] for i in range(801)] for nu in MODES}}
    np.savez_compressed(POOL/'pit.npz', temps=np.array(temps), w_eV=w, wq_eV=np.array(wq), col6=np.array(c6),
                        col7=np.array(c7), col8=np.array(c8))
    (POOL/'pit_text.json').write_text(json.dumps(text))
    # bands, phonons, DFPT smearing scan and script constants
    kb, eb = read_bands(RAW/'Fig2/e_band/band.eig')
    kp, ep = read_bands(RAW/'Fig2/ph_band/matdyn.freq')
    s = (RAW/'Fig2/fig2_d_f_h.py').read_text()
    sig = ast.literal_eval(re.search(r'SIGMA_VALUES = np\.array\((\[[^\]]*\])\)', s).group(1))
    block = re.search(r'OMEGA_EQ1 = np\.array\((.*?)\) \* CMM1_TO_EV', s, re.S).group(1)
    om = ast.literal_eval(re.sub(r'#[^\n]*', '', block))
    a_b = (RAW/'Fig2/fig2_a_b.py').read_text()
    c_e_g = (RAW/'Fig2/fig2_c_e_g.py').read_text()
    archive = {
        'bands': {'k_alat_2pi': kb[:, 2].round(6).tolist(), 'energies_eV': eb.T.round(6).tolist()},
        'phonons': {'q_alat_2pi': kp[:, 2].round(6).tolist(), 'frequencies_cm1': ep.T.round(4).tolist()},
        'dfpt_scan': {'sigma_mRy': sig, 'modes': [1, 20, 23], 'frequency_cm1': om},
        'fermi_shift_band_plot_eV': float(re.search(r'FERMI_SHIFT = (-?[\d.]+)', a_b).group(1)),
        'fermi_shift_comment': re.search(r'FERMI_SHIFT = [^\n]*', a_b).group(0),
        'k_plus_minus_q_comment': re.search(r'In a 1D system[^\n]*\n[^\n]*', c_e_g).group(0).replace('\n', ' '),
        'files': {f: ((RAW/f) if (RAW/f).exists() else (CACHE/f)).read_text() for f in TEXT_FILES},
    }
    assert len(om) == 3 and all(len(r) == len(sig) for r in om)
    (POOL/'archive.json').write_text(json.dumps(archive))
    print(json.dumps({'temperatures': temps, 'n_files': len(TEXT_FILES), 'sigma_mRy': sig, 'scan_cm1': om,
                      'fermi_shift': archive['fermi_shift_band_plot_eV'], 'comment': archive['k_plus_minus_q_comment']}, indent=1))


if __name__ == '__main__':
    main()
