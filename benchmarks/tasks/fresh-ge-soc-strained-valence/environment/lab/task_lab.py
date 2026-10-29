"""Lab logic for fresh-ge-soc-strained-valence: new Quantum ESPRESSO 6.7 calculations for (strained) diamond Ge.

validate(request) checks a structured request against public.json (every field required, no defaults) and returns
its charge, computed from the canonical size of the request: the number of symmetry-distinct k points pw.x will use
for the SCF (counted here exactly as pw.x counts them for the cubic or tetragonal cell), times a factor for
spin-orbit coupling and cutoff, plus the band k points of the optional non-self-consistent run.
run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x -npool CPUS` (one thread per rank), then
optionally a 'bands' run at the requested k points and bands.x with lsym for symmetry labels, and returns raw
physical results; the inputs and outputs are kept as evidence, wavefunctions are deleted. Nothing hidden is ever
served: grade() compares the terminal submission with reference values computed at construction time through
this same run() and stored in private/ (see benchmarks/construction/fresh-ge-soc-strained-valence/).
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
FIXED = SPEC['fixed']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
BOHR_A = 0.529177210903
RY_EV = 13.605693122994
FIELDS = ('kind', 'a_par_A', 'c_A', 'spin_orbit', 'ecutwfc_Ry', 'kmesh', 'kshift', 'occupations', 'degauss_Ry',
          'band_kpoints_inverse_A', 'symmetry_labels')


# ---------------------------------------------------------------- symmetry and canonical request size
def _cubic_ops():
    ops = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, (p, s) in enumerate(zip(perm, signs)):
                R[i, p] = s
            ops.append(R)
    return ops


OPS = _cubic_ops()


def cell_of(a, c):
    """Primitive vectors (rows, Angstrom) of the fcc lattice of the conventional cell a x a x c."""
    return np.array([[0, a / 2, c / 2], [a / 2, 0, c / 2], [a / 2, a / 2, 0]], float)


def irreducible_k(a, c, n, shift):
    """Number of k points pw.x keeps for an n x n x n Monkhorst-Pack mesh (kpoint_grid + irreducible_BZ).

    The diamond space group keeps all point operations of the (cubic or tetragonal) lattice; with inversion in the
    group, time reversal adds nothing."""
    cell = cell_of(a, c)
    D = np.diag([1.0, 1.0, c / a])
    lattice = [R for R in OPS if np.abs(R @ D @ R.T - D).max() < 1e-12]
    B = np.linalg.inv(cell).T
    frac = [np.round(np.linalg.inv(B.T) @ R @ B.T) for R in lattice]
    grid = (np.array(list(itertools.product(range(n), repeat=3)), float) + shift / 2) / n
    key = lambda k: tuple(np.round((k - np.floor(k + 1e-9)) * 2 * n).astype(int) % (2 * n))   # noqa: E731
    index = {key(k): i for i, k in enumerate(grid)}
    parent = list(range(len(grid)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, k in enumerate(grid):
        for M in frac:
            for kk in (M @ k, -(M @ k)):
                j = index.get(key(kk))
                if j is not None:
                    x, y = find(i), find(j)
                    if x != y:
                        parent[max(x, y)] = min(x, y)
    return len({find(i) for i in range(len(grid))}), len(lattice)


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
        raise ValueError(f'missing fields {missing}: every field is required (no defaults; use null or [] where allowed)')
    allowed = SPEC['allowed']
    out = {}
    for key in ('kind', 'spin_orbit', 'ecutwfc_Ry', 'kmesh', 'kshift', 'occupations', 'symmetry_labels'):
        value = request[key]
        ok = isinstance(value, bool) if key in ('spin_orbit', 'symmetry_labels') else not isinstance(value, bool)
        if not ok or value not in allowed[key]:
            raise ValueError(f'{key} must be one of {allowed[key]}')
        out[key] = value
    for key in ('a_par_A', 'c_A'):
        lo, hi = SPEC[key + '_range']
        if not _number(request[key]) or not lo <= request[key] <= hi:
            raise ValueError(f'{key} must be a number in [{lo}, {hi}] (Angstrom)')
        out[key] = round(float(request[key]), 5)
    if abs(out['c_A'] / out['a_par_A'] - 1) > SPEC['max_abs_c_over_a_minus_1'] + 1e-12:
        raise ValueError(f"|c_A / a_par_A - 1| must not exceed {SPEC['max_abs_c_over_a_minus_1']}")
    if out['occupations'] == 'fixed':
        if request['degauss_Ry'] is not None:
            raise ValueError('degauss_Ry must be null for fixed occupations')
        out['degauss_Ry'] = None
    else:
        if isinstance(request['degauss_Ry'], bool) or request['degauss_Ry'] not in allowed['degauss_Ry']:
            raise ValueError(f"degauss_Ry must be one of {allowed['degauss_Ry']} for smearing")
        out['degauss_Ry'] = request['degauss_Ry']
    ks = request['band_kpoints_inverse_A']
    kmax = SPEC['max_abs_k_component_inverse_A']
    if (not isinstance(ks, list) or len(ks) > SPEC['max_band_kpoints']
            or not all(isinstance(k, list) and len(k) == 3 and all(_number(v) and abs(v) <= kmax for v in k) for k in ks)):
        raise ValueError(f"band_kpoints_inverse_A must be a list of at most {SPEC['max_band_kpoints']} points "
                         f"[kx, ky, kz] (Cartesian, inverse Angstrom, components within +-{kmax}); [] for none")
    out['band_kpoints_inverse_A'] = [[round(float(v), 5) + 0.0 for v in k] for k in ks]
    if out['symmetry_labels'] and not ks:
        raise ValueError('symmetry_labels requires a non-empty band_kpoints_inverse_A')
    return out


def size_and_charge(req):
    nk, nsym = irreducible_k(req['a_par_A'], req['c_A'], req['kmesh'], req['kshift'])
    cost = SPEC['charge']
    mode = 'spin_orbit' if req['spin_orbit'] else 'scalar'
    work = nk * cost['scf_factor'][mode][str(req['ecutwfc_Ry'])]
    nb = len(req['band_kpoints_inverse_A'])
    if nb:
        work += cost['band_run_overhead'] + nb * cost['band_factor'][mode][str(req['ecutwfc_Ry'])]
    units = work / cost['seconds_per_unit']
    return nk, nsym, max(1, math.ceil(units - 1e-9))


def validate(request):
    req = normalize(request)
    nk, nsym, charge = size_and_charge(req)
    if charge > SPEC['charge']['max_units_per_request']:
        raise ValueError(f"request too large: {nk} symmetry-distinct k points and {len(req['band_kpoints_inverse_A'])} "
                         f"band k points -> {charge} units; the maximum per request is "
                         f"{SPEC['charge']['max_units_per_request']} (use a smaller mesh, cutoff or band list)")
    return charge


# ---------------------------------------------------------------- pw.x input, execution, parsing
def write_input(req, folder, calculation, kpoints=None):
    a, c = req['a_par_A'], req['c_A']
    soc = req['spin_orbit']
    system = ['ibrav=0, nat=2, ntyp=1', f'celldm(1)={a / BOHR_A:.10f}',
              f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={FIXED['dual'] * req['ecutwfc_Ry']}",
              f"nbnd={FIXED['bands_spin_orbit'] if soc else FIXED['bands_scalar']}"]
    if soc:
        system.append('noncolin=.true., lspinorb=.true.')
    if req['occupations'] == 'fixed':
        system.append("occupations='fixed'")
    else:
        system.append(f"occupations='smearing', smearing='mv', degauss={req['degauss_Ry']}")
    control = [f"calculation='{calculation}', prefix='ge', outdir='./work'", f"pseudo_dir='{HERE}'",
               "disk_io='low'"]
    if calculation == 'scf':
        control.append('tprnfor=.true., tstress=.true.')
    electrons = [f"conv_thr={FIXED['conv_thr_Ry']}, mixing_beta={FIXED['mixing_beta']}, "
                 f"electron_maxstep={FIXED['electron_maxstep']}"]
    if calculation == 'bands':
        electrons.append('diago_full_acc=.true.')
    r = c / a
    lines = ['&CONTROL', *control, '/', '&SYSTEM', *system, '/', '&ELECTRONS', *electrons, '/',
             'ATOMIC_SPECIES', f"Ge {FIXED['mass_amu']} {SPEC['pseudopotentials']['spin_orbit' if soc else 'scalar']['file']}", 'CELL_PARAMETERS alat',
             f'0.0 0.5 {r / 2:.12f}', f'0.5 0.0 {r / 2:.12f}', '0.5 0.5 0.0',
             'ATOMIC_POSITIONS crystal', 'Ge 0.0 0.0 0.0', 'Ge 0.25 0.25 0.25']
    if kpoints is None:
        n, s = req['kmesh'], req['kshift']
        lines += ['K_POINTS automatic', f'{n} {n} {n} {s} {s} {s}', '']
    else:   # tpiba units of 2 pi / a_par
        lines += ['K_POINTS tpiba', str(len(kpoints))]
        lines += ['%.10f %.10f %.10f 1' % tuple(v * a / (2 * math.pi) for v in k) for k in kpoints] + ['']
    name = 'scf.in' if calculation == 'scf' else 'bands.in'
    (folder/name).write_text('\n'.join(lines))
    return name


def _mpirun(folder, program, infile, outfile, npool=True):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    cmd = ['mpirun', '--allow-run-as-root', '--oversubscribe', '-np', str(CPUS), program]
    if npool:
        cmd += ['-npool', str(CPUS)]
    cmd += ['-in', infile]
    with (folder/outfile).open('w') as out, (folder/(outfile + '.err')).open('w') as err:
        try:
            proc = subprocess.run(cmd, cwd=folder, stdout=out, stderr=err, env=env, timeout=FIXED['timeout_s'])
        except subprocess.TimeoutExpired:
            raise ValueError(f'{program} exceeded the time limit; the job is charged') from None
    text = (folder/outfile).read_text(errors='ignore')
    if proc.returncode or 'JOB DONE.' not in text:
        raise ValueError(f'{program} did not finish; the job is charged (see {outfile})')
    return text


def _xml_bands(folder):
    root = ET.parse(folder/'work/ge.save/data-file-schema.xml').getroot()
    out = root.find('output')
    bands = out.find('band_structure')
    alat_A = float(out.find('atomic_structure').get('alat')) * BOHR_A
    ks = []
    for item in bands.findall('ks_energies'):
        k = np.array([float(x) for x in item.findtext('k_point').split()]) * 2 * math.pi / alat_A
        ks.append((k, [round(float(x) * 2 * RY_EV, 5) for x in item.findtext('eigenvalues').split()]))
    return bands, ks


def parse_scf(folder, req):
    text = (folder/'scf.out').read_text(errors='ignore')
    if 'convergence NOT achieved' in text or 'convergence has been achieved' not in text:
        raise ValueError(f"SCF did not converge within {FIXED['electron_maxstep']} iterations; the job is charged")
    energy = float(re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)[-1])
    block = text[text.rfind('total   stress'):].split('\n')
    stress_kbar = [[float(x) for x in block[i].split()[3:6]] for i in (1, 2, 3)]
    bands, ks = _xml_bands(folder)
    res = {'total_energy_Ry': energy,
           'stress_GPa': [[round(v / 10.0, 5) for v in row] for row in stress_kbar],
           'pressure_GPa': round(sum(stress_kbar[i][i] for i in range(3)) / 30.0, 5),
           'irreducible_k_points': int(re.findall(r'number of k points=\s+(\d+)', text)[-1]),
           'symmetry_operations': int(re.findall(r'(\d+) Sym\. Ops\.', text)[-1]),
           'scf_iterations': int(re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)[-1])}
    if req['occupations'] == 'smearing':
        res['fermi_energy_eV'] = round(float(bands.findtext('fermi_energy')) * 2 * RY_EV, 5)
    else:
        res['highest_occupied_level_eV'] = round(float(bands.findtext('highestOccupiedLevel')) * 2 * RY_EV, 5)
        lu = bands.findtext('lowestUnoccupiedLevel')
        if lu is not None:
            res['lowest_unoccupied_level_eV'] = round(float(lu) * 2 * RY_EV, 5)
    k, e = min(ks, key=lambda item: float(np.linalg.norm(item[0])))
    res['mesh_point_nearest_gamma'] = {'k_cartesian_inverse_A': [round(float(v), 6) + 0.0 for v in k],
                                       'band_energies_eV': e}
    return res


def parse_labels(text):
    """Little group and irreducible-representation labels per band k point from bands.x (lsym) output."""
    out = []
    for chunk in text.split('xk=(')[1:]:
        k = [float(x) for x in chunk.split(')')[0].replace(',', ' ').split()]
        group = re.search(r'Band symmetry,\s*(.*?)\s*(double )?point group', chunk)
        rows = []
        for m in re.finditer(r'e\(\s*(\d+)\s*-\s*(\d+)\)\s*=\s*([-\d.]+)\s+eV\s+(\d+)\s+-->\s+(.*)', chunk):
            rows.append({'bands': [int(m.group(1)), int(m.group(2))], 'energy_eV': float(m.group(3)),
                         'label': m.group(5).strip()})
        if 'zone border point and non-symmorphic' in chunk and not rows:
            rows = 'not classified by bands.x (zone-border point of a non-symmorphic group)'
        out.append({'k_tpiba': k, 'point_group': group.group(1).strip() if group else None, 'bands': rows})
    return out


def run(request, workdir):
    req = normalize(request)
    nk_expected, nsym_expected, charge = size_and_charge(req)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    write_input(req, folder, 'scf')
    _mpirun(folder, 'pw.x', 'scf.in', 'scf.out')
    result = parse_scf(folder, req)
    if result['irreducible_k_points'] != nk_expected:
        result['note_k_count'] = f'charged for {nk_expected} k points, pw.x used {result["irreducible_k_points"]}'
    result['scf_wall_seconds'] = round(time.monotonic() - start, 1)
    if req['band_kpoints_inverse_A']:
        write_input(req, folder, 'bands', req['band_kpoints_inverse_A'])
        _mpirun(folder, 'pw.x', 'bands.in', 'bands.out')
        _, ks = _xml_bands(folder)
        result['band_kpoints'] = [{'k_cartesian_inverse_A': kin, 'band_energies_eV': e}
                                  for kin, (_, e) in zip(req['band_kpoints_inverse_A'], ks)]
        if req['symmetry_labels']:
            (folder/'bands_pp.in').write_text("&BANDS\n prefix='ge', outdir='./work', filband='bands.dat', "
                                              "lsym=.true.\n/\n")
            text = _mpirun(folder, 'bands.x', 'bands_pp.in', 'bands_pp.out', npool=False)
            labels = parse_labels(text)
            for item, lab in zip(result['band_kpoints'], labels):
                item['symmetry'] = {'point_group': lab['point_group'], 'levels': lab['bands']}
    result['wall_seconds'] = round(time.monotonic() - start, 1)
    result['cell_vectors_A'] = [[round(v, 10) for v in row] for row in cell_of(req['a_par_A'], req['c_A'])]
    result['settings'] = req
    result['charged_k_points'] = nk_expected
    shutil.rmtree(folder/'work', ignore_errors=True)        # keep inputs/outputs as evidence, drop wavefunctions
    for f in folder.glob('bands.dat*'):
        f.unlink()
    return result


# ---------------------------------------------------------------- grading
def grade(payload, trace, validation):
    private = json.loads((HERE/'private/hidden.json').read_text())
    ref, tol = private['reference'], private['tolerances']
    fields = list(ref)
    missing = [f for f in fields + ['evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bad = [f for f in fields if not _number(payload[f])]
    if bad:
        return False, {'error': f'fields must be finite numbers: {bad}'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for field in fields:
        err = payload[field] - ref[field]
        checks[field] = abs(err) <= tol[field]
        errors[field] = err
    return all(checks.values()), {'checks': checks, 'errors': errors, 'note': 'absolute errors in the field units'}
