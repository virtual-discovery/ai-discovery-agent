"""Fresh Quantum ESPRESSO 6.7 calculations for aluminium slabs and bulk aluminium.

Writes the pw.x input for a validated request, runs `mpirun -np <ranks> pw.x` in the job folder, extracts the
planar-averaged potentials with pp.x, and returns the parsed physical results. Every number is computed live
from the inputs written here; nothing is looked up.
"""
import math
import os
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RY_EV = 13.605693122994
BOHR_A = 0.529177210903
UPF = 'Al.pbe-n-kjpaw_psl.1.0.0.UPF'
MASS = 26.9815385
A0 = 4.04                      # protocol fcc lattice constant (Angstrom)
DUAL = 8                       # ecutrho = 8 ecutwfc (SSSP recommendation for this PAW file)
NP = int(os.environ.get('LAB_MPI_RANKS', '2'))
SQ3 = math.sqrt(3.0)
SQ2 = math.sqrt(2.0)

# Surface cells (in-plane vectors in units of A0), in-plane stacking shifts (crystal) and interlayer spacing.
FACES = {
    '111': {'cell': [[1 / SQ2, 0.0], [1 / (2 * SQ2), SQ3 / (2 * SQ2)]],
            'stack': [(0.0, 0.0), (1 / 3, 1 / 3), (2 / 3, 2 / 3)], 'd': 1 / SQ3},
    '100': {'cell': [[1 / SQ2, 0.0], [0.0, 1 / SQ2]], 'stack': [(0.0, 0.0), (0.5, 0.5)], 'd': 0.5},
    '110': {'cell': [[1.0, 0.0], [0.0, 1 / SQ2]], 'stack': [(0.0, 0.0), (0.5, 0.5)], 'd': 1 / (2 * SQ2)},
}


def surface_cell(face, a=A0):
    f = FACES[face]
    vec = np.array(f['cell']) * a
    area = abs(vec[0, 0] * vec[1, 1] - vec[0, 1] * vec[1, 0])
    return vec, area, f['d'] * a


def slab_geometry(face, layers, vacuum, a=A0):
    """Symmetric slab of `layers` bulk-terminated planes centred in the cell; vacuum = gap between outer planes."""
    vec, area, d = surface_cell(face, a)
    stack = FACES[face]['stack']
    c = (layers - 1) * d + vacuum
    cell = np.zeros((3, 3))
    cell[0, :2], cell[1, :2], cell[2, 2] = vec[0], vec[1], c
    z = c / 2 + (np.arange(layers) - (layers - 1) / 2) * d
    frac = [[stack[n % len(stack)][0], stack[n % len(stack)][1], z[n] / c] for n in range(layers)]
    return cell, frac, area, d, c


def oriented_bulk_geometry(face, a=A0):
    vec, area, d = surface_cell(face, a)
    stack = FACES[face]['stack']
    n = len(stack)
    cell = np.zeros((3, 3))
    cell[0, :2], cell[1, :2], cell[2, 2] = vec[0], vec[1], n * d
    frac = [[stack[i][0], stack[i][1], i / n] for i in range(n)]
    return cell, frac, area, d


def primitive_geometry(a=A0):
    h = a / 2
    cell = np.array([[-h, 0.0, h], [0.0, h, h], [-h, h, 0.0]])
    return cell, [[0.0, 0.0, 0.0]]


def pw_input(calc, cell, frac, kmesh, smearing, ecut, fix=None, conv_thr=1e-9, nbnd=None):
    nat = len(frac)
    lines = ["&CONTROL", f" calculation='{calc}', prefix='al', outdir='./work', pseudo_dir='{HERE}'",
             " disk_io='low', tprnfor=.true."]
    if calc == 'relax':
        lines.append(" nstep=40, forc_conv_thr=2.0d-4, etot_conv_thr=1.0d-5")
    lines += ["/", "&SYSTEM", f" ibrav=0, nat={nat}, ntyp=1",
              f" ecutwfc={ecut}, ecutrho={DUAL * ecut}",
              f" occupations='smearing', smearing='mv', degauss={smearing}"]
    if nbnd:
        lines.append(f" nbnd={nbnd}")
    lines += ["/", "&ELECTRONS", f" conv_thr={conv_thr:.1e}, electron_maxstep=250",
              " mixing_mode='local-TF', mixing_beta=0.5, mixing_ndim=12", "/"]
    if calc == 'relax':
        lines += ["&IONS", " ion_dynamics='bfgs'", "/"]
    lines += ["ATOMIC_SPECIES", f"Al {MASS} {UPF}", "CELL_PARAMETERS angstrom"]
    lines += ['  '.join(f'{v:.10f}' for v in row) for row in cell]
    lines.append("ATOMIC_POSITIONS crystal")
    for i, p in enumerate(frac):
        flags = '' if fix is None else ('  0 0 1' if fix[i] else '  0 0 0')
        lines.append(f"Al {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}{flags}")
    lines += ["K_POINTS automatic", f"{kmesh[0]} {kmesh[1]} {kmesh[2]} 0 0 0", ""]
    return '\n'.join(lines)


def run_program(args, folder, infile, outfile, timeout):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    with open(folder/outfile, 'w') as out, open(folder/(outfile + '.err'), 'w') as err:
        proc = subprocess.run(args + ['-in', infile], cwd=folder, stdout=out, stderr=err, env=env, timeout=timeout)
    text = (folder/outfile).read_text(errors='ignore')
    if proc.returncode or 'JOB DONE.' not in text:
        raise RuntimeError(f'{args[-1] if args else "program"} did not finish (return code {proc.returncode})')
    return text


def mpirun(program, ranks):
    return ['mpirun', '--allow-run-as-root', '--oversubscribe', '--bind-to', 'none', '-np', str(ranks), program, '-npool', str(ranks)]


def parse_pw(text):
    converged = re.findall(r'convergence has been achieved in\s+(\d+) iterations', text)
    if not converged or 'convergence NOT achieved' in text:
        raise ValueError('pw.x SCF did not converge within the iteration limit; no result is returned')
    free = float(re.findall(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text)[-1])
    ts = re.findall(r'smearing contrib\. \(-TS\)\s+=\s+([-\d.]+)\s+Ry', text)
    fermi = float(re.findall(r'the Fermi energy is\s+([-\d.]+)\s+ev', text)[-1])
    nks = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    return {'free_energy_Ry': free, 'free_energy_eV': free * RY_EV,
            'smearing_contribution_minus_TS_eV': float(ts[-1]) * RY_EV if ts else None,
            'fermi_energy_eV': fermi, 'irreducible_k_points': nks,
            'scf_iterations': [int(x) for x in converged]}


def band_summary(folder):
    """Number of bands and the lowest energy reached by the highest computed band over the k-points (eV)."""
    root = ET.parse(folder/'work/al.save/data-file-schema.xml').getroot()
    bs = root.find('./output/band_structure')
    top = [float(k.findtext('eigenvalues').split()[-1]) for k in bs.findall('ks_energies')]
    return {'bands_computed': int(bs.findtext('nbnd')), 'highest_band_minimum_eV': round(min(top) * 2 * RY_EV, 4)}


def parse_forces(text):
    block = text[text.rfind('Forces acting on atoms'):]
    rows = re.findall(r'atom\s+\d+\s+type\s+\d+\s+force =\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', block)
    return [[float(v) * RY_EV / BOHR_A for v in row] for row in rows]


def parse_final_positions(text, nat):
    blocks = re.findall(r'ATOMIC_POSITIONS \(crystal\)\n((?:Al.*\n){%d})' % nat, text)
    if not blocks:
        return None
    return [[float(v) for v in line.split()[1:4]] for line in blocks[-1].strip().split('\n')]


def read_cube_planar(path):
    """Planar average (over the surface cell) of a pp.x cube file; returns z (Angstrom) and V (eV)."""
    lines = Path(path).read_text().split('\n')
    nat = abs(int(lines[2].split()[0]))
    n = [int(lines[3 + i].split()[0]) for i in range(3)]
    dz = float(lines[5].split()[3]) * BOHR_A
    data = np.array(' '.join(lines[6 + nat:]).split(), dtype=float)
    grid = data[:n[0] * n[1] * n[2]].reshape(n[0], n[1], n[2])
    v = grid.mean(axis=(0, 1)) * RY_EV
    return np.arange(n[2]) * dz, v


def potentials(folder, timeout):
    """Electrostatic (ionic + Hartree) and total local (also exchange-correlation) potentials."""
    out = {}
    for name, num in {'es': 11, 'tot': 1}.items():
        (folder/f'pp_{name}.in').write_text(
            f"&INPUTPP\n prefix='al', outdir='./work', filplot='{name}.dat', plot_num={num}\n/\n"
            f"&PLOT\n nfile=1, filepp(1)='{name}.dat', weight(1)=1.0, iflag=3, output_format=6, "
            f"fileout='{name}.cube'\n/\n")
        run_program(['pp.x'], folder, f'pp_{name}.in', f'pp_{name}.out', timeout)
        out[name] = read_cube_planar(folder/f'{name}.cube')
        for suffix in ('.dat', '.cube'):
            (folder/(name + suffix)).unlink(missing_ok=True)
    return out['es'][0], out['es'][1], out['tot'][1]


def run_slab(req, folder, timeout):
    face, layers, vacuum = req['face'], int(req['layers']), float(req['vacuum_A'])
    k, smearing, ecut = int(req['kpoints']), float(req['smearing_Ry']), int(req['ecutwfc_Ry'])
    relax = req['relax']
    a = float(req.get('lattice_A', A0))      # construction-only override; validate() never accepts it
    cell, frac, area, d, c = slab_geometry(face, layers, vacuum, a)
    fix = None
    if relax == 'outer_two_both_sides':      # z of the two outermost planes on each side; symmetry is kept
        fix = [i < 2 or i >= layers - 2 for i in range(layers)]
    calc = 'scf' if fix is None else 'relax'
    (folder/'pw.in').write_text(pw_input(calc, cell, frac, (k, k, 1), smearing, ecut, fix,
                                         conv_thr=2e-10 * layers, nbnd=req.get('nbnd')))
    t0 = time.monotonic()
    text = run_program(mpirun('pw.x', NP), folder, 'pw.in', 'pw.out', timeout)
    res = parse_pw(text)
    pos = frac
    if calc == 'relax':
        if 'bfgs converged' not in text:
            raise ValueError('BFGS relaxation did not converge within 40 steps; no result is returned')
        pos = parse_final_positions(text, layers) or frac
        res['relaxation_steps'] = len(res['scf_iterations'])
    zpos = np.array([p[2] * c for p in pos])
    res.update(band_summary(folder))
    z, v_es, v_tot = potentials(folder, timeout)
    shutil.rmtree(folder/'work', ignore_errors=True)
    spacing = np.diff(zpos)
    res.update({
        'kind': 'slab', 'face': face, 'layers': layers, 'vacuum_A': vacuum, 'kpoints': f'{k}x{k}x1 Gamma-centred',
        'smearing': f'Marzari-Vanderbilt {smearing} Ry', 'ecutwfc_Ry': ecut, 'ecutrho_Ry': DUAL * ecut,
        'relax': relax, 'lattice_constant_A': a,
        'surface_cell_area_A2': round(area, 6), 'cell_height_A': round(c, 6),
        'bulk_interlayer_spacing_A': round(d, 6),
        'atom_z_A': [round(x, 5) for x in zpos],
        'interlayer_spacings_A': [round(x, 5) for x in spacing],
        'forces_eV_per_A': [[round(f, 5) for f in row] for row in parse_forces(text)],
        'potential_z_A': [round(x, 4) for x in z],
        'electrostatic_potential_eV': [round(x, 5) for x in v_es],
        'total_local_potential_eV': [round(x, 5) for x in v_tot],
        'pw_input': (folder/'pw.in').read_text().replace(str(HERE), '/lab'),
        'wall_seconds': round(time.monotonic() - t0, 1)})
    return res


def run_bulk(req, folder, timeout):
    k, smearing, ecut = req['kpoints'], float(req['smearing_Ry']), int(req['ecutwfc_Ry'])
    a = float(req.get('lattice_A', A0))      # construction-only override; validate() never accepts it
    if req['cell'] == 'primitive':
        cell, frac = primitive_geometry(a)
        mesh = (k, k, k)
        extra = {'cell': 'primitive fcc (1 atom)'}
    else:
        cell, frac, area, d = oriented_bulk_geometry(req['face'], a)
        mesh = (k, k, int(req['kz']))
        extra = {'cell': f'oriented bulk: ({req["face"]}) surface cell, {len(frac)} layers', 'face': req['face'],
                 'surface_cell_area_A2': round(area, 6), 'bulk_interlayer_spacing_A': round(d, 6)}
    (folder/'pw.in').write_text(pw_input('scf', cell, frac, mesh, smearing, ecut, conv_thr=1e-11))
    t0 = time.monotonic()
    text = run_program(mpirun('pw.x', NP), folder, 'pw.in', 'pw.out', timeout)
    res = parse_pw(text)
    res.update(band_summary(folder))
    shutil.rmtree(folder/'work', ignore_errors=True)
    res.update({'kind': 'bulk', **extra, 'atoms': len(frac), 'kpoints': 'x'.join(map(str, mesh)) + ' Gamma-centred',
                'smearing': f'Marzari-Vanderbilt {smearing} Ry', 'ecutwfc_Ry': ecut, 'ecutrho_Ry': DUAL * ecut,
                'lattice_constant_A': a, 'free_energy_eV_per_atom': res['free_energy_eV'] / len(frac),
                'pw_input': (folder/'pw.in').read_text().replace(str(HERE), '/lab'),
                'wall_seconds': round(time.monotonic() - t0, 1)})
    return res


def calculate(req, workdir, timeout=1500):
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    try:
        res = run_slab(req, folder, timeout) if req['kind'] == 'slab' else run_bulk(req, folder, timeout)
    finally:      # keep inputs and outputs, never the wavefunctions or charge density
        shutil.rmtree(folder/'work', ignore_errors=True)
        for scratch in folder.glob('*.cube'):
            scratch.unlink()
    for key in ('free_energy_eV', 'fermi_energy_eV'):
        if not math.isfinite(res[key]):
            raise RuntimeError('non-finite result')
    return res
