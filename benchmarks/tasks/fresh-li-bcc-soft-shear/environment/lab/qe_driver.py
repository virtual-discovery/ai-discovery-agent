"""Fresh Quantum ESPRESSO 6.7 pw.x calculations for (strained) bcc lithium, one atom per primitive cell.

Every result is computed live from the inputs written here; nothing is looked up.
"""
import itertools
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
RYB3_GPA = 14710.507848260711        # 1 Ry/bohr^3 in GPa
UPF = 'li_pbe_v1.4.uspp.F.UPF'        # SSSP 1.3 precision (GBRV 1.4), 1s2 2s1 in valence
MASS = 6.94
DUAL = 8                              # ecutrho = 8 * ecutwfc (SSSP recommendation 40/320 Ry)
NP = int(os.environ.get('LAB_MPI_RANKS', '2'))
BASE = 0.5 * np.array([[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]])   # bcc primitive vectors / a (rows)


def strain_tensor(voigt):
    """Symmetric strain tensor from a Voigt 6-vector with engineering shears (e4 = 2 eps_yz, ...)."""
    e1, e2, e3, e4, e5, e6 = (float(x) for x in voigt)
    return np.array([[e1, e6 / 2, e5 / 2], [e6 / 2, e2, e4 / 2], [e5 / 2, e4 / 2, e3]])


def cell_vectors(lattice_A, voigt):
    """Rows are the deformed primitive vectors r' = (1 + eps) r in Angstrom."""
    return lattice_A * BASE @ (np.eye(3) + strain_tensor(voigt)).T


_MATS = [np.array(m).reshape(3, 3) for m in itertools.product((-1, 0, 1), repeat=9)]
_MATS = [m for m in _MATS if round(abs(np.linalg.det(m))) == 1]


def symmetry_operations(voigt, tol=1e-7):
    """Number of point operations of the deformed bcc lattice (= space group of the one-atom crystal)."""
    cell = cell_vectors(1.0, voigt)
    g = cell @ cell.T
    return sum(1 for m in _MATS if np.abs(m @ g @ m.T - g).max() < tol)


def pw_input(cell, kmesh, degauss, ecut):
    lines = ["&CONTROL", " calculation='scf', prefix='li', outdir='./work'", f" pseudo_dir='{HERE}'",
             " tprnfor=.true., tstress=.true., disk_io='low'", "/",
             "&SYSTEM", " ibrav=0, nat=1, ntyp=1", f" ecutwfc={ecut}, ecutrho={DUAL * ecut}",
             f" occupations='smearing', smearing='mv', degauss={degauss}", "/",
             "&ELECTRONS", " conv_thr=1.0d-12, mixing_beta=0.7", "/",
             "ATOMIC_SPECIES", f"Li {MASS} {UPF}", "CELL_PARAMETERS angstrom"]
    lines += ['  '.join(f'{x:.12f}' for x in row) for row in cell]
    lines += ["ATOMIC_POSITIONS crystal", "Li 0.0 0.0 0.0", "K_POINTS automatic", f"{kmesh} {kmesh} {kmesh} 0 0 0", ""]
    return '\n'.join(lines)


def run_mpi(folder, timeout):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    cmd = ['mpirun', '--allow-run-as-root', '-np', str(NP), 'pw.x', '-npool', str(NP), '-in', 'scf.in']
    with open(folder/'scf.out', 'w') as out, open(folder/'scf.err', 'w') as err:
        proc = subprocess.run(cmd, cwd=folder, stdout=out, stderr=err, env=env, timeout=timeout)
    text = (folder/'scf.out').read_text(errors='ignore')
    if proc.returncode or 'JOB DONE.' not in text:
        raise RuntimeError(f'pw.x did not finish (return code {proc.returncode})')
    return text


def parse_pw(text):
    if 'convergence has been achieved' not in text:
        raise RuntimeError('pw.x SCF did not converge')
    free = float(re.findall(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text)[-1])
    internal = float(re.findall(r'internal energy E=F\+TS\s+=\s+([-\d.]+)\s+Ry', text)[-1])
    block = text[text.rfind('total   stress'):].split('\n')
    stress = np.array([[float(x) for x in block[r].split()[0:3]] for r in (1, 2, 3)]) * RYB3_GPA
    return {'free_energy_Ry': free, 'free_energy_eV': free * RY_EV,
            'internal_energy_Ry': internal, 'internal_energy_eV': internal * RY_EV,
            'stress_GPa': stress.round(6).tolist(), 'pressure_GPa': round(float(np.trace(stress) / 3), 6),
            'fermi_energy_eV': float(re.findall(r'the Fermi energy is\s+([-\d.]+)\s+ev', text)[-1]),
            'irreducible_k_points': int(re.findall(r'number of k points=\s+(\d+)', text)[-1]),
            'symmetry_operations': int(re.findall(r'(\d+) Sym\. Ops\.', text)[-1]),
            'scf_iterations': int(re.findall(r'convergence has been achieved in\s+(\d+) iterations', text)[-1])}


def calculate(request, workdir, timeout=3000):
    """Run one validated SCF request; return the physical results."""
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    lattice = float(request['lattice_A'])
    voigt = [float(x) for x in request.get('strain_voigt', [0.0] * 6)]
    kmesh, degauss = int(request['kmesh']), float(request['smearing_Ry'])
    ecut = int(request.get('ecutwfc_Ry', 40))
    cell = cell_vectors(lattice, voigt)
    (folder/'scf.in').write_text(pw_input(cell, kmesh, degauss, ecut))
    t0 = time.monotonic()
    result = {'lattice_A': lattice, 'strain_voigt': voigt, 'cell_vectors_A': cell.round(12).tolist(),
              'cell_volume_A3': float(abs(np.linalg.det(cell))),
              'kmesh': f'{kmesh}x{kmesh}x{kmesh} Gamma-centred', 'smearing': f'Marzari-Vanderbilt {degauss} Ry',
              'ecutwfc_Ry': ecut, 'ecutrho_Ry': DUAL * ecut}
    result.update(parse_pw(run_mpi(folder, timeout)))
    shutil.rmtree(folder/'work', ignore_errors=True)     # keep input and output, drop wavefunctions
    result['wall_seconds'] = round(time.monotonic() - t0, 1)
    if not all(math.isfinite(v) for v in [result['free_energy_Ry'], result['pressure_GPa']]):
        raise RuntimeError('non-finite result')
    return result


if __name__ == '__main__':      # construction helper: python3 qe_driver.py '<request json>' workdir
    import sys
    print(json.dumps(calculate(json.loads(sys.argv[1]), sys.argv[2]), indent=1))
