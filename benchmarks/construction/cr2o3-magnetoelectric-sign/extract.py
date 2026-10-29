"""Parse the unpacked archive (fetch.sh) into pool.json: every archived quantity the lab can serve.

    PYTHONPATH=<dir with h5py> python3 -B extract.py [RAW_DIR]

RAW_DIR defaults to ~/.cache/ai-theorist/cr2o3-magnetoelectric-sign/raw/archive. h5py is needed only to read
the two phonopy .hdf5 files. Nothing is recomputed; the script only parses text/hdf5 outputs and asserts the
consistency checks recorded in README.md.
"""
import glob
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = Path(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/.cache/ai-theorist/cr2o3-magnetoelectric-sign/raw/archive'))
VE = RAW/'VASP_and_Elk'
V = VE/'VASP_calculations'
NAMES = ['Cr1', 'Cr2', 'Cr3', 'Cr4', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6']
NUM = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'


def text(p):
    return Path(p).read_text(errors='replace')


def poscar(p):
    L = text(p).split('\n')
    s = float(L[1].split()[0])
    cell = [[float(x) * s for x in L[i].split()[:3]] for i in (2, 3, 4)]
    species, counts = L[5].split(), [int(x) for x in L[6].split()]
    n = sum(counts)
    frac = [[float(x) for x in L[8 + i].split()[:3]] for i in range(n)]
    return {'comment': L[0].strip(), 'cell_A': cell, 'species': species, 'counts': counts, 'frac': frac}


def amp_value(tag):
    return float(tag.replace('_', '.'))


# ----------------------------------------------------------------------------------------------- VASP Born charges
def born_outcar(p):
    t = text(p)
    def vec(label):
        m = re.findall(label + r'=\(\s*(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s*\)', t)
        return [float(x) for x in m[-1]]
    i = t.index('position of ions in cartesian coordinates')
    pos = [[float(x) for x in l.split()[:3]] for l in t[i:].split('\n')[1:11]]
    lines = [l.strip() for l in t.split('\n') if 'dipole moment' in l or 'p[sp2]' in l][-4:]
    return {'outcar_positions_cartesian_A': pos, 'p_ion': vec(r'p\[ion\]'), 'p_elc': vec(r'p\[elc\]'),
            'p_sp1': vec(r'p\[sp1\]'), 'p_sp2': vec(r'p\[sp2\]'), 'dipole_lines': lines,
            'energy_sigma0_eV': float(re.findall(r'energy\(sigma->0\) =\s*(' + NUM + ')', t)[-1]),
            'finished': 'General timing' in t}


def born_data():
    out = {}
    root = V/'Born_effective_charges'
    for d in sorted(root.glob('Disp_*')):
        tag = d.name[len('Disp_'):]
        if tag == '0_0':
            hdr = [l.strip() for l in text(d/'OUTCAR').split('\n')[:2100] if re.search(r'TITEL|ZVAL   =  |POMASS =  \d', l)]
            out['0_0'] = {'reference': {'poscar': poscar(d/'POSCAR'), **born_outcar(d/'OUTCAR')}, 'outcar_header': hdr}
            continue
        out[tag] = {}
        for n in NAMES:
            for c in 'xyz':
                f = d/f'Pos_{n}_{c}'
                out[tag][f'{n}_{c}'] = {'poscar': poscar(f/'POSCAR'), **born_outcar(f/'OUTCAR')}
    return out


# ----------------------------------------------------------------------------------------------- phonons
def band_gamma(p):
    """Gamma-point block of a phonopy band.yaml: masses, frequencies (THz), real eigenvectors."""
    t = text(p)
    masses = [float(x) for x in re.findall(r'mass:\s*(' + NUM + ')', t)]
    g = t[t.index('phonon:'):]
    g = g[:g.index('- q-position', 20)]
    freqs, vecs = [], []
    for band in re.split(r'\n  - # \d+\n', g)[1:]:
        freqs.append(float(re.search(r'frequency:\s*(' + NUM + ')', band).group(1)))
        comps = re.findall(r'- \[\s*(' + NUM + r'),\s*(' + NUM + r')\s*\]', band)
        v = np.array([[float(a), float(b)] for a, b in comps]).reshape(10, 3, 2)
        assert np.abs(v[..., 1]).max() < 1e-8
        vecs.append(v[..., 0].round(10).tolist())
    symbols = re.findall(r'- symbol:\s*(\w+)', t)
    return {'symbols': symbols, 'masses_amu': masses[:10], 'frequencies_THz': freqs, 'eigenvectors': vecs}


def fc_eigen():
    import h5py
    d = V/'Phonon_calculation/Getting_eigenvectors_force_constant_matrix'
    vals = {k: float(v[()]) for k, v in h5py.File(d/'Cr_ph_eigenval_pol.hdf5', 'r').items()}
    vecs = {k: np.asarray(v[()]).tolist() for k, v in h5py.File(d/'Cr_ph_eigenvec_pol.hdf5', 'r').items()}
    return {'values': vals, 'vectors': vecs}


# ----------------------------------------------------------------------------------------------- Elk moments
def elk_run(d):
    t = text(d/'INFO.OUT')
    blocks = t.split('Moments :')
    last = 'Moments :' + blocks[-1].split('\n \n')[0]
    tot = [float(x) for x in re.findall(r'total moment\s+:\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + ')', t)[-1]]
    mt = [float(x) for x in re.findall(r'total in muffin-tins\s+:\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + ')', t)[-1]]
    inter = [float(x) for x in re.findall(r'interstitial\s+:\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + ')', t)[-1]]
    atoms = re.findall(r'atom\s+(\d+)\s+:\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + ')', last)
    rms = re.findall(r'RMS change in Kohn-Sham potential \(target\) :\s+(' + NUM + r')\s+\(\s*(' + NUM + r')\)', t)
    etot = re.findall(r'total energy\s+:\s+(' + NUM + ')', t)
    inp = text(d/'elk.in')
    a0 = inp.index('atoms')
    atoms_block = inp[a0:inp.index('reducebf')].rstrip()
    sidecars = {}
    for f in sorted(d.glob('*.out')):
        lines = [l for l in text(f).split('\n') if 'I_MPI_PMI_LIBRARY' not in l and l.strip()]
        sidecars[f.name] = lines[-6:]
    return {'elk_in': inp, 'elk_in_atoms_block': atoms_block,
            'moments_block': last.rstrip(),
            'total_moment_muB': tot, 'total_in_muffin_tins_muB': mt, 'interstitial_muB': inter,
            'muffin_tin_atom_moments_muB': [[float(x) for x in a[1:]] for a in atoms],
            'scf_loops': len(re.findall(r'Loop number', t)),
            'messages': sorted(set(m for m in ('Reached self-consistent loops maximum', 'Convergence targets achieved',
                                               'Self-consistent loop stopped') if m in t)),
            'last_rms_potential_change': float(rms[-1][0]) if rms else None,
            'rms_target': float(rms[-1][1]) if rms else None,
            'final_total_energy_Ha': float(etot[-1]),
            'run_completed': 'Elk version 8.4.21 stopped' in t,
            'other_files_in_folder': sidecars}


def elk_data():
    out = {}
    for md in sorted((VE/'Elk_calculations/Induced_magnetic_moments').glob('Eu_*_x')):
        out[md.name] = {d.name[len('Amplitude_'):]: elk_run(d/'Calc_force') for d in sorted(md.glob('Amplitude_*'))}
    return out


# ----------------------------------------------------------------------------------------------- VASP moments
def vasp_tables(t, label):
    """Last per-ion table printed under `label` (e.g. 'magnetization (x)')."""
    i = t.rindex(' ' + label + '\n')
    rows = []
    for l in t[i:].split('\n')[4:14]:
        rows.append([float(x) for x in l.split()[1:5]])
    return rows


def vasp_run(d):
    t = text(d/'OUTCAR')
    tot = re.findall(r'number of electron\s+(' + NUM + r')\s+magnetization\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + ')', t)
    out = {'poscar': poscar(d/'POSCAR'),
           'total_magnetization_muB': [float(x) for x in tot[-1][1:]],
           'total_magnetization_line': [l.strip() for l in t.split('\n') if 'number of electron' in l][-1],
           'electronic_steps': len(re.findall(r'-+ Iteration', t)),
           'finished_normally': 'General timing' in t,
           'last_energy_sigma0_eV': float(re.findall(r'energy\(sigma->0\) =\s*(' + NUM + ')', t)[-1])}
    for c in 'xyz':
        out[f'ion_magnetization_{c}'] = vasp_tables(t, f'magnetization ({c})')
        if f' orbital moment ({c})' in t:
            out[f'ion_orbital_moment_{c}'] = vasp_tables(t, f'orbital moment ({c})')
    return out


def vasp_moment_data():
    root = V/'Induced_magnetic_moments'
    out = {'INCAR': text(root/'INCAR'), 'KPOINTS': text(root/'KPOINTS')}
    for md in sorted(root.glob('Eu_*_x')):
        out[md.name] = {d.name[len('Amplitude_'):]: vasp_run(d) for d in sorted(md.glob('Amplitude_*'))}
    return out


# ----------------------------------------------------------------------------------------------- QE inputs
def qe_data():
    root = RAW/'Quantum_espresso'
    out = {'phonons_at_gamma': {f.name: text(f) for f in sorted((root/'phonons_at_gamma').iterdir())
                                if not f.name.startswith('.') and not f.name.endswith('~')}}
    for md in sorted(root.glob('mode_E_u_*')):
        out[md.name] = {d.name[len('displacement_amplitude_'):]: {f.name: text(f) for f in sorted(d.iterdir())
                                                                  if not f.name.startswith('.') and not f.name.endswith('~')}
                        for d in sorted(md.glob('displacement_amplitude_*'))}
    return out


def main():
    structures = {
        'relaxation_POSCAR': poscar(V/'Structural_relaxation/POSCAR'),
        'relaxation_CONTCAR': poscar(V/'Structural_relaxation/CONTCAR'),
        'relaxation_INCAR': text(V/'Structural_relaxation/INCAR'),
        'phonon_POSCAR': poscar(V/'Phonon_calculation/POSCAR'),
        'phonon_INCAR': text(V/'Phonon_calculation/INCAR'),
        'born_INCAR': text(V/'Born_effective_charges/INCAR'),
        'FORCE_SETS': text(V/'Phonon_calculation/FORCE_SETS'),
        'bands_conf': text(V/'Phonon_calculation/Getting_eigenvectors_force_constant_matrix/bands.conf'),
        'elk_species_note': 'Cr.in and O.in species files (Elk_calculations/) are standard species files.',
    }
    ph = V/'Phonon_calculation/Getting_eigenvectors_force_constant_matrix'
    pool = {'structures': structures, 'born': born_data(), 'fc_eigen': fc_eigen(),
            'phonon_gamma': {k: band_gamma(ph/f'band_{k}.yaml') for k in ('CrO', 'Cr', 'H')},
            'elk': elk_data(), 'vasp_moments': vasp_moment_data(), 'qe': qe_data(),
            'readme': text(RAW.parent.parent/'README.txt') if (RAW.parent.parent/'README.txt').exists() else None}
    # --- consistency checks recorded in README.md
    b = pool['born']
    for key, rec in b['-0_005'].items():          # mislabeled folder: POSCARs say -0.005, OUTCARs are other runs
        k = NAMES.index(key.split('_')[0]); c = 'xyz'.index(key.split('_')[1])
        zval = 14.0 if k < 4 else 6.0
        actual = 0.005 if key.startswith(('O5', 'O6')) else 0.05     # O5/O6: +0.005 runs; all others: +0.05 runs
        assert abs(rec['p_ion'][c] + zval * actual) < 1e-6, key
        assert abs(rec['outcar_positions_cartesian_A'][k][c] - b['0_0']['reference']['outcar_positions_cartesian_A'][k][c]
                   - actual) < 2e-5, key
    for tag in ('0_005', '0_02', '0_05', '-0_02', '-0_05'):
        u = amp_value(tag)
        for key, rec in b[tag].items():
            k = NAMES.index(key.split('_')[0]); c = 'xyz'.index(key.split('_')[1])
            zval = 14.0 if k < 4 else 6.0
            assert abs(rec['p_ion'][c] + zval * u) < 1e-6, (tag, key)
    e = pool['elk']['Eu_3_x']['0_03']
    assert 'Convergence targets achieved' in e['messages'] and e['scf_loops'] < 200
    (HERE/'pool.json').write_text(json.dumps(pool))
    print('pool.json', round((HERE/'pool.json').stat().st_size / 1e6, 2), 'MB')


if __name__ == '__main__':
    main()
