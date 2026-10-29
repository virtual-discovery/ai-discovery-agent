"""Exploration driver for diamond Ge (pw.x scf/relax with strain and sublattice displacement)."""
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HERE = Path(os.environ.get('GE_PSEUDO_DIR', '/pseudo'))
BOHR = 0.529177210903
RY_EV = 13.605693122994
NP = int(os.environ.get('LAB_MPI_PROCS', '2'))
PSEUDOS = {'gbrv': ('ge_pbe_v1.4.uspp.F.UPF', 14), 'psl4': ('Ge.pbe-n-rrkjus_psl.1.0.0.UPF', 4),
           'psl14': ('Ge.pbe-dn-rrkjus_psl.1.0.0.UPF', 14)}
MASS = 72.630


def geometry(a, strain, disp):
    e1, e2, e3, e4, e5, e6 = strain
    eps = np.array([[e1, e6 / 2, e5 / 2], [e6 / 2, e2, e4 / 2], [e5 / 2, e4 / 2, e3]])
    F = np.eye(3) + eps
    A0 = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
    A = A0 @ F.T
    u = np.array(disp, float)
    r1 = -u / 2
    r2 = F @ (0.25 * a * np.ones(3)) + u / 2
    frac = np.linalg.solve(A.T, np.array([r1, r2]).T).T
    return A, frac


def pw_input(req):
    upf, zval = PSEUDOS[req['pseudo']]
    A, frac = geometry(req['lattice_A'], req['strain_voigt'], req['displacement_A'])
    nocc = zval  # 2 atoms x zval / 2
    occ = req['occupations']
    system = ['ibrav=0, nat=2, ntyp=1', f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={8 * req['ecutwfc_Ry']}",
              f'nbnd={nocc + 4}']
    if occ == 'fixed':
        system.append("occupations='fixed'")
    else:
        system.append(f"occupations='smearing', smearing='mv', degauss={req['degauss_Ry']}")
    if req.get('nosym'):
        system.append('nosym=.true.')
    control = [f"calculation='{req['kind']}', prefix='ge', outdir='./work'", f"pseudo_dir='{HERE}'",
               "tprnfor=.true., tstress=.true., disk_io='low'"]
    if req['kind'] == 'relax':
        control.append(f"forc_conv_thr={req['forc_conv_thr']}, etot_conv_thr={req.get('etot_conv_thr', 1e-5)}, nstep=40")
    lines = ['&CONTROL', *control, '/', '&SYSTEM', *system, '/', '&ELECTRONS',
             f"conv_thr={req.get('conv_thr', 1e-10)}, mixing_beta=0.5, electron_maxstep=200", '/']
    if req['kind'] == 'relax':
        lines += ['&IONS', "ion_dynamics='bfgs'", '/']
    lines += ['ATOMIC_SPECIES', f'Ge {MASS} {upf}', 'CELL_PARAMETERS angstrom']
    lines += ['  '.join(f'{v:.12f}' for v in row) for row in A]
    lines += ['ATOMIC_POSITIONS crystal'] + ['Ge ' + '  '.join(f'{v:.12f}' for v in row) for row in frac]
    n, s = req['kmesh'], req['kshift']
    lines += ['K_POINTS automatic', f'{n} {n} {n} {s} {s} {s}', '']
    return '\n'.join(lines)


def parse(folder, req):
    text = (folder/'pw.out').read_text(errors='ignore')
    if 'JOB DONE.' not in text:
        raise RuntimeError('pw.x did not finish')
    if 'convergence NOT achieved' in text:
        raise RuntimeError('SCF not converged')
    e = float(re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)[-1])
    block = text[text.rfind('total   stress'):].split('\n')
    stress = [[float(x) for x in block[i].split()[3:6]] for i in (1, 2, 3)]
    forces = [[float(x) for x in m.groups()] for m in
              re.finditer(r'atom\s+\d+ type\s+\d+\s+force =\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', text)][-2:]
    nk = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    nsym = re.findall(r'(\d+) Sym\. Ops\.', text)
    iters = [int(x) for x in re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)]
    root = ET.parse(folder/'work/ge.save/data-file-schema.xml').getroot()
    out = root.find('output')
    cell = np.array([[float(x) for x in out.find('atomic_structure/cell').find(t).text.split()] for t in ('a1', 'a2', 'a3')]) * BOHR
    pos = np.array([[float(x) for x in at.text.split()] for at in out.find('atomic_structure/atomic_positions')]) * BOHR
    bands = out.find('band_structure')
    res = {'energy_Ry': e, 'stress_kbar': stress, 'forces_Ry_bohr': forces, 'nk_irr': nk,
           'nsym': int(nsym[-1]) if nsym else None, 'scf_iterations': iters,
           'positions_A': pos.tolist(), 'cell_A': cell.tolist()}
    ef = bands.find('fermi_energy')
    if ef is not None:
        res['fermi_eV'] = float(ef.text) * 2 * RY_EV
    ho = bands.find('highestOccupiedLevel')
    if ho is not None:
        res['homo_eV'] = float(ho.text) * 2 * RY_EV
    lu = bands.find('lowestUnoccupiedLevel')
    if lu is not None:
        res['lumo_eV'] = float(lu.text) * 2 * RY_EV
    best = None
    for ks in bands.findall('ks_energies'):
        k = np.array([float(x) for x in ks.find('k_point').text.split()])
        if best is None or np.linalg.norm(k) < best[0]:
            best = (np.linalg.norm(k), k, [float(x) * 2 * RY_EV for x in ks.find('eigenvalues').text.split()])
    res['k_nearest_gamma'] = best[1].tolist()
    res['bands_nearest_gamma_eV'] = best[2]
    return res


def run(req, workdir):
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    (folder/'pw.in').write_text(pw_input(req))
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMPI_MCA_btl_vader_single_copy_mechanism='none')
    t0 = time.monotonic()
    with open(folder/'pw.out', 'w') as out, open(folder/'pw.err', 'w') as err:
        subprocess.run(['mpirun', '--allow-run-as-root', '-np', str(NP), 'pw.x', '-npool', str(NP), '-in', 'pw.in'],
                       cwd=folder, stdout=out, stderr=err, env=env, timeout=3600)
    res = parse(folder, req)
    res['wall_s'] = round(time.monotonic() - t0, 1)
    shutil.rmtree(folder/'work', ignore_errors=True)
    return res


if __name__ == '__main__':
    batch = json.load(open(sys.argv[1]))
    outfile = sys.argv[2]
    done = set()
    if os.path.exists(outfile):
        done = {json.loads(l)['tag'] for l in open(outfile) if l.strip()}
    for item in batch:
        if item['tag'] in done:
            continue
        row = {'tag': item['tag'], 'request': item['request']}
        try:
            row['result'] = run(item['request'], os.path.join(sys.argv[3], item['tag']))
        except Exception as exc:
            row['error'] = f'{type(exc).__name__}: {exc}'
        with open(outfile, 'a') as f:
            f.write(json.dumps(row) + '\n')
        r = row.get('result', {})
        print(item['tag'], r.get('wall_s'), r.get('nk_irr'), r.get('energy_Ry'), row.get('error', ''), flush=True)
