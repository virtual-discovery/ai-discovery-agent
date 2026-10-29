"""Reference and ablation calculations for 2d-qsh-shc-plateau (runs in the runtime image, one material per call).

    python3 -B reference.py MODELS_DIR SRC_ROOT OUT_DIR material_folder

Uses the lab's own evaluator (qsh_tb.py) on the lab's npz models. Writes OUT_DIR/<name>.json with:
  band edges (300^2 mesh), Fermi energies printed in the archive, Z2 from Wilson loops at several resolutions,
  SHC (x, y, z) at the Wannier mid-gap and at the archive's Fermi energies, for the reference settings
  (wannier_centres, eta = 1e-4 eV, 180^2 + 8x8 refinement) and for the ablation settings.
"""
import json
import re
import sys
import time
from pathlib import Path
import numpy as np
import build_models
import qsh_tb


def z2_count(w):
    """Parity of the crossings of the largest-gap midpoint by the hybrid Wannier centres (Soluyanov-Vanderbilt)."""
    def ref(a):
        a = np.sort(a)
        g = np.diff(np.r_[a, a[0] + 1])
        j = int(np.argmax(g))
        return (a[j] + g[j] / 2) % 1, g[j]
    total, mingap = 0, 1.0
    for i in range(len(w) - 1):
        z0, g0 = ref(w[i])
        z1, _ = ref(w[i + 1])
        mingap = min(mingap, g0)
        b = np.sort(w[i + 1])
        lo, hi = min(z0, z1), max(z0, z1)
        total += int(np.sum((b > lo) & (b < hi))) if hi - lo < 0.5 else int(np.sum((b > hi) | (b < lo)))
    return int(total % 2), float(mingap)


def archive_facts(f):
    out = (f/'aiida.out').read_text()
    win = (f/'aiida.win').read_text()
    ne = int(round(float(re.search(r'number of electrons\s*=\s*(\S+)', out).group(1))))
    ex = re.search(r'exclude_bands\s*=\s*([^\n]+)', win).group(1)
    nex = sum(int(b) - int(a) + 1 for a, b in re.findall(r'(\d+)\s*-\s*(\d+)', ex))
    get = lambda key: float(re.search(r'# %s = (\S+)' % key, win).group(1))
    return {'nelec': ne, 'n_excluded': nex, 'nocc': ne - nex,
            'scf_fermi_eV': float(re.search(r'the Fermi energy is\s+(\S+)', out).group(1)),
            'win_fermi_energy_eV': float(re.search(r'\nfermi_energy\s*=\s*(\S+)', win).group(1)),
            'win_VBM_eV': get('VBM'), 'win_CBM_eV': get('CBM'), 'win_gap_eV': get('Bandgap')}


def main():
    models, src, out, folder = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
    f = src/folder
    name = folder.split('/')[-1].replace('|', '_')
    m = qsh_tb.Model(models/f'{name}.npz', build_models.cell_of(f))
    res = {'name': name, 'folder': folder, 'num_wann': m.nw, 'c_A': m.c, 'area_A2': m.area, **archive_facts(f)}
    t0 = time.time()
    be = m.band_edges(res['nocc'], (300, 300))
    res['band_edges'] = be
    mid = 0.5 * (be['vbm'][0] + be['cbm'][0])
    res['mid_gap_eV'] = mid
    res['win_mid_gap_eV'] = 0.5 * (res['win_VBM_eV'] + res['win_CBM_eV'])
    res['t_edges'] = round(time.time() - t0)
    # Z2 at several Wilson-loop resolutions, both loop axes
    z2 = {}
    for nk, nl in ((21, 40), (51, 80), (101, 100), (101, 200), (201, 200), (201, 300)):
        for axis in (0, 1):
            w = m.wcc(res['nocc'], np.linspace(0, 0.5, nk), nl, axis)
            z2[f'{nk}x{nl}_axis{axis}'] = z2_count(w)
    res['z2'] = z2
    mus = [mid, res['scf_fermi_eV'], res['win_mid_gap_eV'], res['win_fermi_energy_eV']]
    res['mus'] = mus
    runs = {'reference': ((180, 180), 1e-4, 'wannier_centres', (8, 10.0, 6000)),
            'eta5meV': ((180, 180), 5e-3, 'wannier_centres', (8, 10.0, 6000)),
            'lattice': ((180, 180), 1e-4, 'lattice', (8, 10.0, 6000)),
            'mesh120_eta1meV': ((120, 120), 1e-3, 'wannier_centres', None),
            'mesh120_eta5meV': ((120, 120), 5e-3, 'wannier_centres', None),
            'mesh60_eta1meV': ((60, 60), 1e-3, 'wannier_centres', None),
            'mesh240_eta1meV': ((240, 240), 1e-3, 'wannier_centres', None)}
    res['shc'] = {}
    for label, (mesh, eta, conv, ref) in runs.items():
        t0 = time.time()
        sheet, w90, nref = m.shc(mesh, [0, 1, 2], eta, mus, conv, ref)
        res['shc'][label] = {'sheet_xyz': sheet.tolist(), 'w90_xyz': w90.tolist(), 'nref': nref, 'seconds': round(time.time() - t0)}
        print(name, label, np.round(sheet[:, 0], 4).tolist(), nref, round(time.time() - t0), flush=True)
    (out/f'{name}.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps({k: res[k] for k in ('name', 'nocc', 'mid_gap_eV', 'z2')}), flush=True)


if __name__ == '__main__':
    main()
