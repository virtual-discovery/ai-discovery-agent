"""Fresh Quantum ESPRESSO 6.7 calculations for zincblende InP (pw.x SCF and ph.x DFPT at Gamma).

Every result is computed live from the inputs written here; nothing is looked up.
"""
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RY_EV = 13.605693122994
BOHR_A = 0.529177210903
SPECIES = [('In', 114.818, 'In.pbe-dn-rrkjus_psl.0.2.2.UPF'), ('P', 30.973762, 'P.pbe-n-rrkjus_psl.1.0.0.UPF')]
DUAL = 8            # ecutrho = 8 * ecutwfc (SSSP precision recommendation for these ultrasoft files)
NBND = 13           # 9 occupied bands (In 4d10 5s2 5p1, P 3s2 3p3) plus 4 empty bands for the band edges
NP = int(os.environ.get('LAB_MPI_RANKS', '2'))


def cell_vectors(lattice_A, strain_zz):
    a = lattice_A / 2.0
    cell = np.array([[0.0, a, a], [a, 0.0, a], [a, a, 0.0]])
    cell[:, 2] *= 1.0 + strain_zz
    return cell


def pw_input(lattice_A, strain_zz, kmesh, kshift, ecutwfc):
    cell = cell_vectors(lattice_A, strain_zz)
    lines = ["&CONTROL", " calculation='scf', prefix='inp', outdir='./work'", f" pseudo_dir='{HERE}'",
             " tprnfor=.true., tstress=.true., disk_io='low'", "/",
             "&SYSTEM", " ibrav=0, nat=2, ntyp=2", f" ecutwfc={ecutwfc}, ecutrho={DUAL * ecutwfc}",
             f" nbnd={NBND}, occupations='fixed'", "/",
             "&ELECTRONS", " conv_thr=1.0d-12, mixing_beta=0.5", "/", "ATOMIC_SPECIES"]
    lines += [f"{el} {mass} {upf}" for el, mass, upf in SPECIES]
    lines += ["CELL_PARAMETERS angstrom"] + ['  '.join(f'{v:.10f}' for v in row) for row in cell]
    lines += ["ATOMIC_POSITIONS crystal", "In 0.00 0.00 0.00", "P  0.25 0.25 0.25",
              "K_POINTS automatic", f"{kmesh} {kmesh} {kmesh} {kshift} {kshift} {kshift}", ""]
    return '\n'.join(lines), cell


PH_INPUT = """phonons and dielectric response at Gamma
&INPUTPH
 prefix='inp', outdir='./work', fildyn='inp.dyn', tr2_ph=1.0d-15, epsil=.true., trans=.{trans}.
/
0.0 0.0 0.0
"""


def run_mpi(program, folder, infile, outfile, timeout):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    cmd = ['mpirun', '--allow-run-as-root', '-np', str(NP), program, '-npool', str(NP), '-in', infile]
    with open(folder/outfile, 'w') as out, open(folder/(outfile + '.err'), 'w') as err:
        proc = subprocess.run(cmd, cwd=folder, stdout=out, stderr=err, env=env, timeout=timeout)
    text = (folder/outfile).read_text(errors='ignore')
    if proc.returncode or 'JOB DONE.' not in text:
        raise RuntimeError(f'{program} did not finish (return code {proc.returncode})')
    return text


def parse_pw(text):
    if 'convergence has been achieved' not in text:
        raise RuntimeError('pw.x SCF did not converge')
    energy = float(re.findall(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text)[-1])
    block = text[text.rfind('total   stress'):].split('\n')
    pressure = float(re.findall(r'P=\s*([-\d.]+)', block[0])[0])
    stress = [[float(x) for x in block[i].split()[3:6]] for i in (1, 2, 3)]
    homo, lumo = map(float, re.findall(r'highest occupied, lowest unoccupied level \(ev\):\s+([-\d.]+)\s+([-\d.]+)', text)[-1])
    nks = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    nsym = re.findall(r'(\d+) Sym\. Ops\.', text)
    return {'total_energy_Ry': energy, 'total_energy_eV': energy * RY_EV,
            'pressure_GPa': pressure / 10.0, 'stress_GPa': [[v / 10.0 for v in row] for row in stress],
            'highest_occupied_eV': homo, 'lowest_unoccupied_eV': lumo, 'band_gap_on_mesh_eV': round(lumo - homo, 4),
            'irreducible_k_points': nks, 'symmetry_operations': int(nsym[-1]) if nsym else None}


def parse_dyn(path):
    """Read the Gamma dynamical matrix, dielectric tensor and Born charges from the ph.x fildyn file."""
    lines = Path(path).read_text().split('\n')
    i = next(k for k, line in enumerate(lines) if 'Dynamical  Matrix in cartesian axes' in line)
    phi = np.zeros((6, 6))
    j = i + 4
    for _ in range(4):
        s1, s2 = (int(x) - 1 for x in lines[j].split())
        for r in range(3):
            vals = [float(x) for x in lines[j + 1 + r].split()]
            phi[3 * s1 + r, 3 * s2:3 * s2 + 3] = vals[0::2]
        j += 4
    k = next(k for k, line in enumerate(lines) if 'Dielectric Tensor:' in line)
    eps = [[float(x) for x in lines[k + 2 + r].split()] for r in range(3)]
    k = next(k for k, line in enumerate(lines) if 'Effective Charges E-U' in line)
    zstar = []
    for atom in range(2):
        start = k + 2 + atom * 4
        zstar.append([[float(x) for x in lines[start + 1 + r].split()] for r in range(3)])
    freqs = [float(x) for x in re.findall(r'freq \(\s*\d+\) =\s+[-\d.]+ \[THz\] =\s+([-\d.]+) \[cm-1\]', '\n'.join(lines))]
    return phi, eps, zstar, freqs


def parse_ph_out(text):
    """Dielectric tensor and Born charges (dF/dE) printed by ph.x (used when no phonons are computed)."""
    lines = text.split('\n')
    k = max(i for i, line in enumerate(lines) if 'Dielectric constant in cartesian axis' in line)
    eps = [[float(x) for x in lines[k + 2 + r].replace('(', ' ').replace(')', ' ').split()] for r in range(3)]
    k = max(i for i, line in enumerate(lines) if 'Effective charges (d Force / dE) in cartesian axis' in line)
    zstar = []
    for atom in range(2):
        start = k + 2 + atom * 4
        zstar.append([[float(x) for x in lines[start + 1 + r].replace('(', ' ').replace(')', ' ').split()[1:]]
                      for r in range(3)])
    return eps, zstar


def calculate(request, workdir, timeout=3000):
    """Run pw.x (and ph.x for kind 'dfpt') for one validated request; return the physical results."""
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    lattice, strain = float(request['lattice_A']), float(request.get('strain_zz', 0.0))
    kmesh, kshift = int(request['kmesh']), int(request['kshift'])
    ecut = int(request.get('ecutwfc_Ry', 50))
    text, cell = pw_input(lattice, strain, kmesh, kshift, ecut)
    (folder/'scf.in').write_text(text)
    t0 = time.monotonic()
    result = {'kind': request['kind'], 'lattice_A': lattice, 'strain_zz': strain,
              'cell_vectors_A': cell.round(10).tolist(), 'cell_volume_A3': float(abs(np.linalg.det(cell))),
              'atoms_crystal': {'In': [0.0, 0.0, 0.0], 'P': [0.25, 0.25, 0.25]},
              'kmesh': f'{kmesh}x{kmesh}x{kmesh} shift {kshift} {kshift} {kshift}',
              'ecutwfc_Ry': ecut, 'ecutrho_Ry': DUAL * ecut}
    result.update(parse_pw(run_mpi('pw.x', folder, 'scf.in', 'scf.out', timeout)))
    result['pw_seconds'] = round(time.monotonic() - t0, 1)
    if request['kind'] in ('dfpt', 'dielectric'):
        (folder/'ph.in').write_text(PH_INPUT.format(trans='true' if request['kind'] == 'dfpt' else 'false'))
        t1 = time.monotonic()
        ph_text = run_mpi('ph.x', folder, 'ph.in', 'ph.out', timeout)
        if request['kind'] == 'dfpt':
            phi, eps, zstar, freqs = parse_dyn(folder/'inp.dyn')
            result.update({'force_constants_gamma_Ry_per_bohr2': phi.round(10).tolist(),
                           'atomic_masses_amu': {el: mass for el, mass, _ in SPECIES},
                           'gamma_frequencies_analytic_raw_cm-1': freqs})
        else:
            eps, zstar = parse_ph_out(ph_text)
        result.update({'dielectric_tensor_eps_inf': eps, 'born_charges_raw': {'In': zstar[0], 'P': zstar[1]},
                       'ph_seconds': round(time.monotonic() - t1, 1)})
    shutil.rmtree(folder/'work', ignore_errors=True)     # keep inputs and outputs, drop wavefunctions and scratch
    result['wall_seconds'] = round(time.monotonic() - t0, 1)
    if not all(math.isfinite(v) for v in [result['total_energy_Ry'], result['pressure_GPa']]):
        raise RuntimeError('non-finite result')
    return result


if __name__ == '__main__':      # construction helper: python3 qe_driver.py '<request json>' workdir
    import sys
    print(json.dumps(calculate(json.loads(sys.argv[1]), sys.argv[2]), indent=1))
