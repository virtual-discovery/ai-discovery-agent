"""Reduce the fetched archive members to compact arrays (pool/) for make_task.py.

Needs h5py (not in the runtime image); run once in a throwaway virtualenv:
    python3 -m venv /tmp/h5env && /tmp/h5env/bin/pip install h5py numpy
    /tmp/h5env/bin/python -B extract.py [RAW] [POOL]
Defaults: RAW = ~/.cache/ai-theorist/skutterudite-wigner-coherence/raw, POOL = .../pool.

Per phono3py run it writes pool/<run>.npz with the archived per-mode arrays on the irreducible
q-mesh (frequency, gamma, gamma_isotope, heat_capacity, qpoint, weight, temperature) and the velocity
data: for phono3py 2.1.0 runs G[q,s,s'] = (1/3) sum_a |V^a_ss'(q)|^2 at the irreducible q (the trace of
the archived star-summed gv_by_gv_operator divided by the weight); for the 1.19 runs v2[q,s] = (1/3)
|v_s(q)|^2 from the archived group velocities. runs.json holds metadata, launch commands and the
thermal-conductivity tables parsed from each run's log; harmonic.json holds element-resolved PDOS,
mean-square displacements and Gamma-point frequencies.
"""
import json
import re
import sys
from pathlib import Path
import h5py
import numpy as np
from fetch import RUNS, RUN_FILES, HARMONIC

CACHE = Path.home()/'.cache/ai-theorist/skutterudite-wigner-coherence'
RAW = Path(sys.argv[1]) if len(sys.argv) > 1 else CACHE/'raw'
POOL = Path(sys.argv[2]) if len(sys.argv) > 2 else CACHE/'pool'
COMPOUND = {r: r.split('/')[0] for r in RUNS}


def log_tables(text):
    """Last K_P/K_C/K_T (or kappa) tables printed in a phono3py log: {label: {T: xx}}."""
    tables = {}
    for label, T, xx in re.findall(r'^(K_P|K_C|K_T)[=\s]\s*([\d.]+)\s+([-\d.]+)', text, re.M):
        tables.setdefault(label, {})[str(int(float(T)))] = float(xx)
    return tables


def header(text):
    keep = []
    for line in text.splitlines()[:60]:
        s = line.strip()
        if re.match(r'^\d+\.\d+\.\d+$', s):
            keep.append('phono3py version ' + s)
        elif any(k in s for k in ('Supercell (dim)', 'Phonon supercell', 'Mesh sampling', 'BZ integration', 'Reading fc',
                                  'Cutoff frequency', 'Run mode')):
            keep.append(s)
    seen, out = set(), []
    for k in keep:
        if k not in seen:
            seen.add(k); out.append(k)
    mass = re.search(r'Mass variance parameters:([^\n]*)', text)
    if mass:
        out.append('Mass variance parameters:' + mass.group(1).rstrip())
    return out


def yaml_cell(text):
    lat = re.search(r'primitive_cell:\s*\n\s*lattice:\s*\n((?:\s*- \[.*\n){3})', text).group(1)
    lattice = [[float(x) for x in re.findall(r'-?\d+\.\d+', row.split('#')[0])] for row in lat.strip().splitlines()]
    block = text.split('primitive_cell:')[1].split('supercell')[0]
    symbols = re.findall(r'symbol: (\w+)', block)
    masses = [float(m) for m in re.findall(r'mass: ([\d.]+)', block)]
    version = re.search(r'version: ([\d.]+)', text).group(1)
    factor = float(re.search(r'frequency_unit_conversion_factor: ([\d.]+)', text).group(1))
    conf = dict(re.findall(r'^\s{4}(\w+): "([^"]*)"', text, re.M))
    return {'lattice_A': lattice, 'symbols': symbols, 'masses_amu': masses, 'version': version, 'configuration': conf,
            'frequency_unit_conversion_factor': factor}


def main():
    POOL.mkdir(parents=True, exist_ok=True)
    meta = {}
    for run, folder in RUNS.items():
        base = RAW/folder
        h5 = next(f for f in RUN_FILES[run] if f.startswith('kappa'))
        f = h5py.File(base/h5, 'r')
        d = {k: f[k][()] for k in ('frequency', 'gamma', 'heat_capacity', 'qpoint', 'weight', 'temperature', 'mesh')}
        d['weight'] = d['weight'].astype(np.int64)
        d['mesh'] = d['mesh'].astype(np.int64)
        d['kappa_unit_conversion'] = np.array(f['kappa_unit_conversion'][()])
        if 'gamma_isotope' in f:
            d['gamma_isotope'] = f['gamma_isotope'][()]
        if 'gv_by_gv_operator' in f:
            g = f['gv_by_gv_operator'][()]
            tr = (g[..., 0] + g[..., 1] + g[..., 2]) / 3.0
            assert np.abs(tr.imag).max() < 1e-8 * np.abs(tr.real).max()
            d['G'] = tr.real / d['weight'][:, None, None]
            d['kappa_P_RTA'] = f['kappa_P_RTA'][()][:, 0]
            d['kappa_C'] = f['kappa_C'][()][:, 0]
        else:
            gv = f['gv_by_gv'][()]
            d['v2'] = (gv[..., 0] + gv[..., 1] + gv[..., 2]) / 3.0 / d['weight'][:, None]
            d['kappa_hdf5'] = f['kappa'][()][:, 0]
        np.savez_compressed(POOL/(run.replace('/', '__') + '.npz'), **d)
        files = {}
        for name in RUN_FILES[run]:
            p = base/name
            if name.endswith('.out'):
                text = p.read_text(errors='replace')
                files[name] = {'kind': 'phono3py log', 'header': header(text), 'tables': log_tables(text)}
            elif name.endswith('.sh'):
                files[name] = {'kind': 'launch script',
                               'commands': [l.strip() for l in p.read_text().splitlines() if l.strip().startswith('phono3py')
                                            or l.strip().startswith('cp ')]}
            elif name == 'tc_data.dat':
                files[name] = {'kind': 'post-processed table', 'text': p.read_text()}
            elif name.endswith('.py'):
                text = p.read_text()
                start = text.index('def calc_tc_raw_matrix')
                end = text.index("root='data_save'", start)
                files[name] = {'kind': 'post-processing script (excerpt: the function that wrote tc_data.dat)',
                               'text': text[start:end]}
        cell = yaml_cell((base/'phono3py.yaml').read_text())
        meta[run] = {'compound': COMPOUND[run], 'archive_folder': 'materials_cloud_archive_submission/' + folder,
                     'kappa_file': h5, 'kappa_md5_group': None, 'cell': cell,
                     'n_irreducible_q': int(len(d['weight'])), 'n_bands': int(d['frequency'].shape[1]),
                     'mesh': d['mesh'].tolist(), 'temperatures_K': d['temperature'].astype(int).tolist(),
                     'has_gamma_isotope': 'gamma_isotope' in d, 'has_velocity_operator': 'G' in d, 'files': files}
    # byte-identical kappa files (archive duplicates)
    import hashlib
    digests = {r: hashlib.md5((RAW/RUNS[r]/meta[r]['kappa_file']).read_bytes()).hexdigest() for r in RUNS}
    for r in RUNS:
        meta[r]['kappa_md5_group'] = digests[r]
    (POOL/'runs.json').write_text(json.dumps(meta, indent=1) + '\n')
    # harmonic data
    harm = {}
    for comp, (pdos_dir, msd_dir) in HARMONIC.items():
        poscar = (RAW/msd_dir/'POSCAR').read_text().splitlines()
        if poscar[5].split()[0].isdigit():      # VASP 4 style: species on the comment line
            species, counts = poscar[0].split(), [int(x) for x in poscar[5].split()]
        else:
            species, counts = poscar[5].split(), [int(x) for x in poscar[6].split()]
        order = [s for s, n in zip(species, counts) for _ in range(n)]
        pd = np.loadtxt(RAW/pdos_dir/'projected_dos.dat')
        pdos = {s: pd[:, 1:][:, [i for i, o in enumerate(order) if o == s]].sum(axis=1) for s in species}
        y = (RAW/msd_dir/'thermal_displacements_calcs/mesh_50/thermal_displacements.yaml').read_text()
        msd = {}
        for block in y.split('- temperature:')[1:]:
            T = str(int(float(block.split()[0])))
            rows = [[float(x) for x in r] for r in re.findall(r'\[\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\s*\]', block)]
            msd[T] = {s: float(np.mean([sum(rows[i]) for i, o in enumerate(order) if o == s])) for s in species}
        pyaml = (RAW/pdos_dir/'phonopy.yaml').read_text()
        settings = {'calculator': re.search(r'calculator: "?(\w+)', pyaml).group(1),
                    'frequency_unit_conversion_factor': float(re.search(r'frequency_unit_conversion_factor: ([\d.]+)', pyaml).group(1)),
                    'length_unit': re.search(r'length: "?(\w+)', pyaml).group(1)}
        cmd = RAW/pdos_dir/'cmd.sh'
        msd_cmd = RAW/msd_dir/'thermal_displacements_calcs/mesh_50/1_cmd_MSD.sh'
        harm[comp] = {'pdos_frequency_axis': pd[:, 0].tolist(), 'pdos_by_element': {s: v.tolist() for s, v in pdos.items()},
                      'pdos_run': {'command': cmd.read_text().strip() if cmd.exists() else None, 'phonopy_yaml': settings,
                                   'file': 'projected_dos.dat (first column as written; per-atom columns summed by element)'},
                      'msd_A2': msd, 'msd_run': {'command': msd_cmd.read_text().strip()},
                      'atoms': dict(zip(species, counts))}
    for tag, path in (('YbFe4Sb12-AFM', 'Fe_based/YbFe4Sb12/AFM/2_phonons/post_processing_matdyn_crystal/matdyn__crystal.freq'),
                      ('YbFe4Sb12-NM', 'Fe_based/YbFe4Sb12/NM/2_phonons/post_processing_matdyn_crystal/matdyn__crystal.freq')):
        lines = (RAW/path).read_text().splitlines()
        nb = int(re.search(r'nbnd=\s*(\d+)', lines[0]).group(1))
        vals = []
        for line in lines[2:]:
            vals += [float(x) for x in line.split()]
            if len(vals) >= nb:
                break
        harm.setdefault(tag, {})['gamma_frequencies_dfpt_cm1'] = vals[:nb]
    (POOL/'harmonic.json').write_text(json.dumps(harm) + '\n')
    print('wrote', sorted(p.name for p in POOL.iterdir()))


if __name__ == '__main__':
    main()
