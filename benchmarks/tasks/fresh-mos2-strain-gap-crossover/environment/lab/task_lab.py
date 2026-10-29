"""Lab logic for fresh-mos2-strain-gap-crossover: new Quantum ESPRESSO 6.7 calculations for monolayer 1H-MoS2.

validate(request) checks a structured request against public.json and returns its charge, computed from the
canonical size of the request: the number of symmetry-distinct k points of the self-consistent mesh (counted here
exactly as pw.x counts them for the strained, displaced cell), times 3 for a relaxation, plus the number of
band-structure k points, times a factor for spin-orbit coupling and for the cutoff.
run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x` (threads pinned to 1) for the scf or
relaxation, pp.x for the planar-averaged electrostatic potential, and pw.x again for the requested band k points;
it returns the raw physical results and keeps pw.in/pw.out (and the pp.x and bands files) as evidence.
Nothing hidden is ever served: grade() compares the terminal submission with reference values computed at
construction time through this same run() and stored in private/ (benchmarks/construction/<task id>/).
"""
import itertools
import json
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
PUBLIC = json.loads((HERE/'public.json').read_text())
SPEC = PUBLIC['request']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
BOHR_A = 0.529177210903
RY_EV = 13.605693122994
SQ3 = math.sqrt(3.0)
FIELDS = ('kind', 'relativistic', 'lattice_A', 'strain', 'sulfur_height_A', 'sulfur_shift_A', 'cell_height_A',
          'ecutwfc_Ry', 'kmesh', 'band_kpoints')
SPECIES = ('Mo', 'S')
N_ELECTRONS = 26                       # Mo 4s2 4p6 4d5 5s1 (14) + 2 x S 3s2 3p4 (6)


# ---------------------------------------------------------------- geometry and symmetry
def geometry(a, eps2, h, shift, c):
    """Cell (rows, Angstrom), Cartesian and crystal positions (Mo, S top, S bottom) of the strained monolayer.

    The unstrained cell is a1 = a(1, 0, 0), a2 = a(-1/2, sqrt3/2, 0), a3 = (0, 0, c): x is a zigzag direction,
    y the armchair direction. The in-plane strain maps r -> (1 + eps) r (cell vectors and ideal sites); Mo sits at
    crystal (1/3, 2/3, 1/2), the S atoms at the strained ideal site (2/3, 1/3) displaced in-plane by `shift`
    (Cartesian, both S atoms together) and at z_Mo +- h.
    """
    F = np.eye(3)
    F[:2, :2] += np.array(eps2, float)
    cell = np.array([[a, 0.0, 0.0], [-a / 2, a * SQ3 / 2, 0.0], [0.0, 0.0, c]]) @ F.T
    mo = np.array([1 / 3, 2 / 3, 0.5]) @ cell
    s0 = np.array([2 / 3, 1 / 3, 0.5]) @ cell + np.array([shift[0], shift[1], 0.0])
    cart = np.array([mo, s0 + [0, 0, h], s0 - [0, 0, h]])
    frac = np.linalg.solve(cell.T, cart.T).T
    return cell, cart, frac


def _lattice_ops(cell, tol=1e-6):
    """Integer matrices M (crystal basis, acting on column coordinates) that preserve the metric."""
    G = cell @ cell.T
    ops = []
    for entries in itertools.product((-1, 0, 1), repeat=9):
        M = np.array(entries, float).reshape(3, 3)
        if abs(abs(np.linalg.det(M)) - 1) > 1e-9:
            continue
        if np.abs(M.T @ G @ M - G).max() < tol * np.abs(G).max():
            ops.append(M)
    return ops


def crystal_ops(cell, frac, species, tol=1e-5):
    """Point operations (crystal basis) of the lattice and of the crystal (with the fractional translation that
    maps the Mo atom onto itself, the only Mo of the cell)."""
    lattice = _lattice_ops(cell)
    crystal = []
    for M in lattice:
        t = frac[0] - M @ frac[0]
        ok = True
        for i, x in enumerate(frac):
            y = M @ x + t
            if not any(species[j] == species[i] and np.abs((y - frac[j]) - np.round(y - frac[j])).max() < tol
                       for j in range(len(frac))):
                ok = False
                break
        if ok:
            crystal.append(M)
    return lattice, crystal


def irreducible_k(cell, frac, n):
    """Number of k points pw.x keeps for an n x n x 1 Gamma-centred mesh: orbits of the mesh under the crystal
    point group and time reversal (pw.x reduces by the lattice group and re-expands stars for the crystal)."""
    _, crystal = crystal_ops(cell, frac, ('Mo', 'S', 'S'))
    kops = [np.round(np.linalg.inv(M).T) for M in crystal]
    grid = [np.array([i / n, j / n, 0.0]) for i in range(n) for j in range(n)]
    key = lambda k: tuple(np.round((k - np.floor(k + 1e-9)) * n).astype(int) % n)       # noqa: E731
    seen, orbits = set(), 0
    for k in grid:
        if key(k) in seen:
            continue
        orbits += 1
        for M in kops:
            for kk in (M @ k, -(M @ k)):
                seen.add(key(kk))
    return orbits, len(crystal)


# ---------------------------------------------------------------- request checking and charge
def _number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def normalize(request):
    """Return the canonical request, or raise ValueError with a usable message (invalid requests are free)."""
    if not isinstance(request, dict):
        raise ValueError('JSON object required; see GET /task for the request schema')
    unknown, missing = set(request) - set(FIELDS), [f for f in FIELDS if f not in request]
    if unknown:
        raise ValueError(f'unknown fields {sorted(unknown)}; allowed: {list(FIELDS)}')
    if missing:
        raise ValueError(f'missing fields {missing}: every field is required (no defaults)')
    allowed = SPEC['allowed']
    out = {}
    for key in ('kind', 'relativistic', 'cell_height_A', 'ecutwfc_Ry', 'kmesh'):
        value = request[key]
        if isinstance(value, bool) or value not in allowed[key]:
            raise ValueError(f'{key} must be one of {allowed[key]}')
        out[key] = value
    if out['kind'] == 'relax' and out['relativistic'] != 'scalar':
        raise ValueError('relaxations are scalar-relativistic only (kind "relax" needs relativistic "scalar")')
    lo, hi = SPEC['lattice_A_range']
    if not _number(request['lattice_A']) or not lo <= request['lattice_A'] <= hi:
        raise ValueError(f'lattice_A must be a number in [{lo}, {hi}]')
    out['lattice_A'] = round(float(request['lattice_A']), 6)
    eps = request['strain']
    if (not isinstance(eps, list) or len(eps) != 2 or not all(isinstance(r, list) and len(r) == 2 for r in eps)
            or not all(_number(v) for r in eps for v in r)):
        raise ValueError('strain must be a 2x2 list of numbers (symmetric in-plane strain tensor, x zigzag, y armchair)')
    eps = np.round(np.array(eps, float), 6)
    if abs(eps[0, 1] - eps[1, 0]) > 1e-9:
        raise ValueError('strain must be symmetric (eps_xy = eps_yx)')
    if np.abs(eps).max() > SPEC['max_abs_strain_component']:
        raise ValueError(f"each strain component must lie within +-{SPEC['max_abs_strain_component']}")
    if np.any((np.abs(eps) > 0) & (np.abs(eps) < SPEC['min_nonzero_strain_component'])):
        raise ValueError(f"nonzero strain components must have magnitude >= {SPEC['min_nonzero_strain_component']}")
    out['strain'] = (eps + 0.0).tolist()
    lo, hi = SPEC['sulfur_height_A_range']
    if not _number(request['sulfur_height_A']) or not lo <= request['sulfur_height_A'] <= hi:
        raise ValueError(f'sulfur_height_A must be a number in [{lo}, {hi}]')
    out['sulfur_height_A'] = round(float(request['sulfur_height_A']), 6)
    sh = request['sulfur_shift_A']
    if not isinstance(sh, list) or len(sh) != 2 or not all(_number(v) for v in sh):
        raise ValueError('sulfur_shift_A must be a list of two numbers [dx, dy] (Angstrom)')
    sh = np.round(np.array(sh, float), 6)
    if np.abs(sh).max() > SPEC['max_abs_sulfur_shift_A']:
        raise ValueError(f"each sulfur_shift_A component must lie within +-{SPEC['max_abs_sulfur_shift_A']} A")
    out['sulfur_shift_A'] = (sh + 0.0).tolist()
    kp = request['band_kpoints']
    if not isinstance(kp, list) or len(kp) > SPEC['max_band_kpoints']:
        raise ValueError(f"band_kpoints must be a list of at most {SPEC['max_band_kpoints']} points [k1, k2]")
    pts = []
    for k in kp:
        if not isinstance(k, list) or len(k) != 2 or not all(_number(v) and abs(v) <= 1.0 for v in k):
            raise ValueError('each band k point must be [k1, k2]: crystal coordinates (units of the strained '
                             'reciprocal vectors b1, b2), each within [-1, 1]')
        pts.append([round(float(v), 6) for v in k])
    out['band_kpoints'] = pts
    return out


def size_and_charge(req):
    cell, _, frac = geometry(req['lattice_A'], req['strain'], req['sulfur_height_A'], req['sulfur_shift_A'],
                             req['cell_height_A'])
    nk, nsym = irreducible_k(cell, frac, req['kmesh'])
    cost = SPEC['charge']
    size = nk * (cost['relax_multiplier'] if req['kind'] == 'relax' else 1)
    size += cost['band_kpoint_weight'] * len(req['band_kpoints'])
    units = size * cost['relativistic_factor'][req['relativistic']] * cost['cutoff_factor'][str(req['ecutwfc_Ry'])]
    return nk, nsym, max(1, math.ceil(units / cost['k_points_per_unit'] - 1e-9))


def validate(request):
    req = normalize(request)
    nk, nsym, charge = size_and_charge(req)
    if charge > SPEC['charge']['max_units_per_request']:
        raise ValueError(f"request too large: {charge} units (maximum {SPEC['charge']['max_units_per_request']} "
                         f"per request); use fewer band k points, a smaller mesh or cutoff")
    return charge


# ---------------------------------------------------------------- pw.x / pp.x input, execution, parsing
def pw_input(req, calculation, cell, frac, kpoints):
    fx = SPEC['fixed']
    rel = req['relativistic']
    nbnd = (N_ELECTRONS // 2 + fx['empty_bands']) * (2 if rel == 'full' else 1)
    control = [f"calculation='{calculation}', prefix='mos2', outdir='./work', pseudo_dir='{HERE}'",
               "tprnfor=.true., tstress=.true., disk_io='low', verbosity='high'"]
    if calculation == 'relax':
        control.append(f"forc_conv_thr={fx['relax_forc_conv_thr_Ry_per_bohr']}, "
                       f"etot_conv_thr={fx['relax_etot_conv_thr_Ry']}, nstep={fx['relax_max_steps']}")
    system = ['ibrav=0, nat=3, ntyp=2', f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={fx['dual'] * req['ecutwfc_Ry']}",
              f"nbnd={nbnd}, occupations='fixed'"]
    if rel == 'full':
        system.append('noncolin=.true., lspinorb=.true.')
    lines = ['&CONTROL', *control, '/', '&SYSTEM', *system, '/', '&ELECTRONS',
             f"conv_thr={fx['conv_thr_Ry']}, mixing_mode='local-TF', mixing_beta={fx['mixing_beta']}, "
             f"electron_maxstep={fx['electron_maxstep']}", '/']
    if calculation == 'relax':
        lines += ['&IONS', "ion_dynamics='bfgs'", '/']
    lines += ['ATOMIC_SPECIES']
    for el in SPECIES:
        lines.append(f"{el} {fx['mass_amu'][el]} {SPEC['pseudopotentials'][rel][el]}")
    lines += ['CELL_PARAMETERS angstrom'] + ['  '.join(f'{v:.12f}' for v in row) for row in cell]
    lines += ['ATOMIC_POSITIONS crystal']
    lines += [f'{el} ' + '  '.join(f'{v:.12f}' for v in row) for el, row in zip(('Mo', 'S', 'S'), frac)]
    if kpoints is None:
        n = req['kmesh']
        lines += ['K_POINTS automatic', f'{n} {n} 1 0 0 0', '']
    else:
        lines += ['K_POINTS crystal', str(len(kpoints))] + [f'{k[0]:.10f} {k[1]:.10f} 0.0 1.0' for k in kpoints] + ['']
    return '\n'.join(lines)


def _run(args, folder, infile, outfile, timeout):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    with (folder/outfile).open('w') as out, (folder/(outfile + '.err')).open('w') as err:
        try:
            proc = subprocess.run(args + ['-in', infile], cwd=folder, stdout=out, stderr=err, env=env, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise ValueError(f'{outfile}: exceeded the time limit; the job is charged') from None
    text = (folder/outfile).read_text(errors='ignore')
    if proc.returncode or 'JOB DONE.' not in text:
        raise ValueError(f'{outfile}: the program did not finish; the job is charged')
    return text


def _mpi(program):
    return ['mpirun', '--allow-run-as-root', '--oversubscribe', '-np', str(CPUS), program, '-npool', str(CPUS)]


def parse_scf(text, req, cell):
    if 'convergence NOT achieved' in text:
        raise ValueError('SCF did not converge; the job is charged (see pw.out)')
    energies = [float(x) for x in re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)]
    block = text[text.rfind('total   stress'):].split('\n')
    stress_kbar = [[float(x) for x in block[i].split()[3:6]] for i in (1, 2, 3)]
    forces = [[[float(x) for x in m.groups()] for m in
               re.finditer(r'atom\s+\d+ type\s+\d+\s+force =\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', chunk)][:3]
              for chunk in text.split('Forces acting on atoms')[1:]]
    nk = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    nsym = re.findall(r'(\d+) Sym\. Ops\.', text)
    iterations = [int(x) for x in re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)]
    homo_lumo = re.findall(r'highest occupied, lowest unoccupied level \(ev\):\s+([-\d.]+)\s+([-\d.]+)', text)
    res = {'total_energy_Ry': energies[-1],
           'stress_GPa': [[round(v / 10.0, 6) for v in row] for row in stress_kbar],
           'forces_Ry_per_bohr': forces[-1],
           'scf_mesh': {'irreducible_k_points': nk, 'symmetry_operations': int(nsym[-1]) if nsym else None},
           'scf_iterations': iterations if req['kind'] == 'relax' else iterations[-1]}
    if homo_lumo:
        res['scf_mesh']['highest_occupied_level_eV'] = float(homo_lumo[-1][0])
        res['scf_mesh']['lowest_unoccupied_level_eV'] = float(homo_lumo[-1][1])
    if req['kind'] == 'relax':
        res['ionic_steps'] = [{'total_energy_Ry': e, 'max_force_component_Ry_per_bohr': max(abs(v) for at in f for v in at)}
                              for e, f in zip(energies, forces)]
        res['relax_converged'] = 'bfgs converged' in text
        if not res['relax_converged']:
            raise ValueError('BFGS relaxation did not converge; the job is charged (see pw.out)')
    return res


def final_positions(folder, cell):
    root = ET.parse(folder/'work/mos2.save/data-file-schema.xml').getroot()
    out = root.find('output')
    pos_bohr = np.array([[float(x) for x in at.text.split()] for at in out.find('atomic_structure/atomic_positions')])
    cart = pos_bohr * BOHR_A
    return cart, np.linalg.solve(cell.T, cart.T).T


def planar_average(path):
    """Planar average over the cell of a pp.x cube file: z (Angstrom) and V (eV)."""
    lines = Path(path).read_text().split('\n')
    nat = abs(int(lines[2].split()[0]))
    n = [int(lines[3 + i].split()[0]) for i in range(3)]
    dz = float(lines[5].split()[3]) * BOHR_A
    data = np.array(' '.join(lines[6 + nat:]).split(), dtype=float)
    grid = data[:n[0] * n[1] * n[2]].reshape(n[0], n[1], n[2])
    return np.arange(n[2]) * dz, grid.mean(axis=(0, 1)) * RY_EV


def bands_from_xml(folder, cell):
    root = ET.parse(folder/'work/mos2.save/data-file-schema.xml').getroot()
    out = root.find('output')
    alat_A = float(out.find('atomic_structure').get('alat')) * BOHR_A
    B = np.linalg.inv(cell).T            # rows b_i / (2 pi)
    pts = []
    for ks in out.find('band_structure').findall('ks_energies'):
        kc = np.array([float(x) for x in ks.findtext('k_point').split()]) / alat_A     # Cartesian, 1/A (/2pi)
        crys = np.linalg.solve(B.T, kc)
        pts.append({'k_crystal': [round(v, 6) for v in crys[:2]],
                    'k_cartesian_inverse_A': [round(v * 2 * math.pi, 6) for v in kc[:2]],
                    'energies_eV': [round(float(x) * 2 * RY_EV, 5) for x in ks.findtext('eigenvalues').split()]})
    return pts


def run(request, workdir):
    req = normalize(request)
    nk_expected, nsym_expected, _ = size_and_charge(req)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    timeout = SPEC['fixed']['timeout_s']
    cell, cart, frac = geometry(req['lattice_A'], req['strain'], req['sulfur_height_A'], req['sulfur_shift_A'],
                                req['cell_height_A'])
    start = time.monotonic()
    try:
        (folder/'pw.in').write_text(pw_input(req, req['kind'], cell, frac, None))
        text = _run(_mpi('pw.x'), folder, 'pw.in', 'pw.out', timeout)
        res = parse_scf(text, req, cell)
        cart, frac = final_positions(folder, cell)
        (folder/'pp.in').write_text("&INPUTPP\n prefix='mos2', outdir='./work', filplot='ves.dat', plot_num=11\n/\n"
                                    "&PLOT\n nfile=1, filepp(1)='ves.dat', weight(1)=1.0, iflag=3, output_format=6, "
                                    "fileout='ves.cube'\n/\n")
        _run(['pp.x'], folder, 'pp.in', 'pp.out', timeout)
        z, v = planar_average(folder/'ves.cube')
        bands = []
        if req['band_kpoints']:
            (folder/'bands.in').write_text(pw_input(req, 'bands', cell, frac, req['band_kpoints']))
            _run(_mpi('pw.x'), folder, 'bands.in', 'bands.out', timeout)
            bands = bands_from_xml(folder, cell)
    finally:      # keep inputs and outputs as evidence, never the wavefunctions, density or cube files
        shutil.rmtree(folder/'work', ignore_errors=True)
        for scratch in ('ves.dat', 'ves.cube'):
            (folder/scratch).unlink(missing_ok=True)
    mo, st, sb = cart
    res.update({
        'cell_vectors_A': [[round(v, 10) for v in row] for row in cell],
        'final_positions_cartesian_A': {'Mo': [round(v, 6) for v in mo], 'S_top': [round(v, 6) for v in st],
                                        'S_bottom': [round(v, 6) for v in sb]},
        'final_positions_crystal': [[round(v, 10) for v in row] for row in frac],
        'final_sulfur_height_A': round((st[2] - sb[2]) / 2, 6),
        'final_sulfur_shift_A': [round(v, 6) for v in ((st[:2] + sb[:2]) / 2 - (np.array([2 / 3, 1 / 3, 0.5]) @ cell)[:2])],
        'occupied_bands': N_ELECTRONS if req['relativistic'] == 'full' else N_ELECTRONS // 2,
        'electrostatic_potential_planar_average': {'z_A': [round(x, 4) for x in z], 'V_eV': [round(x, 5) for x in v]},
        'band_kpoints_results': bands,
        'settings': req, 'charged_k_points': nk_expected, 'charged_symmetry_operations': nsym_expected,
        'wall_seconds': round(time.monotonic() - start, 1)})
    return res


# ---------------------------------------------------------------- grading
def _point(x):
    return isinstance(x, list) and len(x) == 2 and all(_number(v) for v in x)


def k_matches(k, ref, frame, tol):
    """True when crystal point k lies within tol (1/A) of a symmetry image of ref, modulo reciprocal vectors.

    frame: {'b_inverse_A': rows b1, b2 (Cartesian, 1/A, including 2 pi), 'k_ops': 2x2 integer matrices acting on
    crystal k (the crystal point group; time reversal is added here)}."""
    b = np.array(frame['b_inverse_A'], float)
    k, ref = np.array(k, float), np.array(ref, float)
    for M in frame['k_ops']:
        for sign in (1, -1):
            d = k - sign * (np.array(M, float) @ ref)
            d -= np.round(d)
            for n1 in (-1, 0, 1):
                for n2 in (-1, 0, 1):
                    if np.linalg.norm((d + [n1, n2]) @ b) <= tol:
                        return True
    return False


def grade(payload, trace, validation):
    private = json.loads((HERE/'private/hidden.json').read_text())
    ref, tol, frames = private['reference'], private['tolerances'], private['frames']
    fields = ['a0_A', 'direct_gap_strain_range_percent', 'edge_leaving_K_at_lower_end', 'edge_leaving_K_at_upper_end',
              'deformation_potentials_eV', 'uniaxial_zigzag_compression', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    rng, dp, uni = payload['direct_gap_strain_range_percent'], payload['deformation_potentials_eV'], payload['uniaxial_zigzag_compression']
    edges = {f: payload[f] for f in ('edge_leaving_K_at_lower_end', 'edge_leaving_K_at_upper_end')}
    form = (_number(payload['a0_A']) and _point(rng)
            and all(isinstance(e, dict) and e.get('band') in ('valence', 'conduction') and _point(e.get('k_crystal'))
                    for e in edges.values())
            and isinstance(dp, dict) and all(_number(dp.get(k)) for k in ('valence_K', 'conduction_K', 'gap_K'))
            and isinstance(uni, dict) and _number(uni.get('transverse_strain_percent')) and _number(uni.get('gap_eV'))
            and _point(uni.get('vbm_k_crystal')) and _point(uni.get('cbm_k_crystal')))
    if not form:
        return False, {'error': 'submission does not follow the deliverable format (numbers, [lower, upper] strains, '
                                'band names "valence"/"conduction", k points as [k1, k2])'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    errors['a0_A'] = payload['a0_A'] - ref['a0_A']
    checks['a0_A'] = abs(errors['a0_A']) <= tol['a0_A']
    # The ends of the direct-gap range are compared as lattice parameters a0 (1 + strain), with the submitted a0: an
    # admissible a0 estimate then does not move the band crossings, which are graded on their own.
    for i, end in enumerate(('lower', 'upper')):
        a_sub = payload['a0_A'] * (1 + rng[i] / 100)
        a_ref = ref['a0_A'] * (1 + ref['direct_gap_strain_range_percent'][i] / 100)
        errors[f'direct_gap_{end}_end_lattice_A'] = a_sub - a_ref
        checks[f'direct_gap_{end}_end'] = abs(a_sub - a_ref) <= tol['crossover_lattice_A']
    for f, e in edges.items():
        r = ref[f]
        checks[f] = e['band'] == r['band'] and k_matches(e['k_crystal'], r['k_crystal'], frames['biaxial'], tol['k_inverse_A'])
    for k in ('valence_K', 'conduction_K', 'gap_K'):
        errors['dp_' + k] = dp[k] / ref['deformation_potentials_eV'][k] - 1
        checks['deformation_potential_' + k] = abs(errors['dp_' + k]) <= tol['deformation_potential_relative'][k]
    ru = ref['uniaxial_zigzag_compression']
    errors['uniaxial_transverse_strain_percent'] = uni['transverse_strain_percent'] - ru['transverse_strain_percent']
    checks['uniaxial_transverse_strain'] = abs(errors['uniaxial_transverse_strain_percent']) <= tol['transverse_strain_percent']
    errors['uniaxial_gap_eV'] = uni['gap_eV'] - ru['gap_eV']
    checks['uniaxial_gap'] = abs(errors['uniaxial_gap_eV']) <= tol['gap_eV']
    for key in ('vbm_k_crystal', 'cbm_k_crystal'):
        checks['uniaxial_' + key.split('_')[0]] = k_matches(uni[key], ru[key], frames['uniaxial'], tol['k_inverse_A'])
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'note': 'relative errors for deformation potentials, absolute otherwise'}
