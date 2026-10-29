"""Parse the 490 archived Quantum ESPRESSO runs of Materials Cloud record gpxwh-dx645 into one pool.

    python3 -B extract.py [CACHE]      (run in ai-theorist-benchmark-runtime:20261007; CACHE defaults to
                                        ~/.cache/ai-theorist/hydrogen-plasma-facing-materials)

Writes CACHE/pool.json (every run: input settings, cell, initial and final positions, energies, convergence
messages, run header) and ./derived.json (energies and settings only, small, kept in the repository).
Nothing is recomputed; values are read from pwscf.in / pwscf.out as released.
"""
import json
import re
import sys
from collections import Counter
from pathlib import Path

CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/hydrogen-plasma-facing-materials'
RAW = CACHE/'raw'
HERE = Path(__file__).resolve().parent
FLOAT = r'[-+]?\d*\.?\d+(?:[dDeE][-+]?\d+)?'


def num(text):
    return float(text.replace('d', 'e').replace('D', 'e'))


def ident(rel):
    """Catalog id for an archive folder."""
    p = rel.split('/')
    if p[0] == 'H-IFE':
        if p[1] == 'hydrogen':
            return 'H-IFE/' + p[2]
        if p[1] == 'supercell':
            return f'H-IFE/{p[2]}/supercell'
        return f'H-IFE/{p[2]}/site-{p[3]}'
    kind = {'single_atom': 'atom', 'unitcell': 'unitcell', 'slab': 'slab', 'slab_supercell': 'slab_supercell',
            'slab_supercell_surface_vacancy': 'vacancy'}[p[1]]
    if kind == 'unitcell':
        return f'SBE/{p[2]}/unitcell'
    if kind == 'vacancy':
        return f'SBE/{p[2]}/vacancy/{p[3]}/{p[4]}'
    return f'SBE/{p[2]}/{kind}/{p[3]}'


def block3(text, start):
    lines = text[start:].split('\n')[1:4]
    return [[float(x) for x in l.split()[:3]] for l in lines]


def parse(folder):
    inp = (folder/'pwscf.in').read_text()
    out = (folder/'pwscf.out').read_text(errors='ignore')
    get = lambda key: (re.search(r'\b%s\s*=\s*([^\s,]+)' % key, inp) or [None, None])[1]
    nat = int(get('nat'))
    assert re.search(r'CELL_PARAMETERS\s+\{?angstrom', inp), folder
    assert re.search(r'ATOMIC_POSITIONS\s+\{?crystal', inp), folder
    cell = block3(inp, inp.index('CELL_PARAMETERS'))
    i = inp.index('ATOMIC_POSITIONS')
    pos0 = [[l.split()[0]] + [round(float(x), 7) for x in l.split()[1:4]] for l in inp[i:].split('\n')[1:1 + nat]]
    species = dict(re.findall(r'^\s*([A-Z][a-z]?)\s+[\d.]+\s+(\S+\.(?:UPF|upf))\s*$', inp, re.M))
    kp = inp[inp.index('K_POINTS'):].split('\n')
    kpts = 'gamma' if 'gamma' in kp[0].lower() else kp[1].strip()
    mags = {int(k): num(v) for k, v in re.findall(r'starting_magnetization\((\d+)\)\s*=\s*(%s)' % FLOAT, inp)}
    s = lambda key: (get(key) or '').strip("'\"") or None
    record = {'input': {
        'calculation': s('calculation'), 'ecutwfc_Ry': num(get('ecutwfc')), 'ecutrho_Ry': num(get('ecutrho')),
        'k_points': kpts, 'occupations': s('occupations'), 'smearing': s('smearing'),
        'degauss_Ry': num(get('degauss')) if get('degauss') else None, 'nspin': int(get('nspin') or 1),
        'starting_magnetization': [mags[k] for k in sorted(mags)] or None,
        'conv_thr_Ry': num(get('conv_thr')), 'etot_conv_thr_Ry': num(get('etot_conv_thr')) if get('etot_conv_thr') else None,
        'forc_conv_thr_Ry_per_bohr': num(get('forc_conv_thr')) if get('forc_conv_thr') else None,
        'electron_maxstep': int(get('electron_maxstep')) if get('electron_maxstep') else None,
        'mixing_beta': num(get('mixing_beta')) if get('mixing_beta') else None,
        'pseudopotentials': species, 'nat': nat, 'species_counts': dict(Counter(p[0] for p in pos0))},
        'cell_A': cell, 'positions_crystal_initial': pos0}
    head = re.search(r'Program PWSCF (v\.\S+) starts on\s+(.*?)\s*$', out, re.M)
    procs = re.search(r'running on\s+(\d+) processors', out)
    wall = re.findall(r'PWSCF\s+:\s+(.*?)CPU\s+(.*?)WALL', out)
    record['run'] = {'program': 'PWSCF ' + head.group(1) if head else None,
                     'started': re.sub(r'\s+', ' ', head.group(2)) if head else None,
                     'mpi_processes': int(procs.group(1)) if procs else None,
                     'wall_time': wall[-1][1].strip() if wall else None}
    energies = [float(x) for x in re.findall(r'^!\s+total energy\s+=\s+(%s) Ry' % FLOAT, out, re.M)]
    final = re.findall(r'Final (?:energy|enthalpy)\s+=\s+(%s) Ry' % FLOAT, out)
    relax = None
    m = re.search(r'bfgs converged in\s+(\d+) scf cycles and\s+(\d+) bfgs steps', out)
    if m:
        relax = f'bfgs converged in {m.group(1)} scf cycles and {m.group(2)} bfgs steps'
    elif 'The maximum number of steps has been reached' in out:
        relax = 'The maximum number of steps has been reached.'
    force = re.findall(r'Total force =\s+(%s)' % FLOAT, out)
    mag = re.findall(r'total magnetization\s+=\s+(%s) Bohr' % FLOAT, out)
    fermi = re.findall(r'the Fermi energy is\s+(%s) ev' % FLOAT, out)
    press = re.findall(r'P=\s*(%s)' % FLOAT, out)
    o = {'last_total_energy_Ry': energies[-1] if energies else None, 'n_scf_cycles': len(energies),
         'final_energy_Ry': float(final[-1]) if final else None, 'relaxation_message': relax,
         'scf_not_converged_warnings': out.count('convergence NOT achieved'),
         'final_total_force_Ry_per_bohr': float(force[-1]) if force else None,
         'total_magnetization_muB_per_cell': float(mag[-1]) if mag else None,
         'fermi_energy_eV': float(fermi[-1]) if fermi else None,
         'pressure_kbar': float(press[-1]) if press else None, 'job_done': 'JOB DONE' in out}
    fin = out.rfind('Begin final coordinates')
    blocks = list(re.finditer(r'ATOMIC_POSITIONS \(crystal\)\n', out))
    if record['input']['calculation'] in ('relax', 'vc-relax') and blocks:
        b = blocks[-1].end()
        o['positions_crystal_final'] = [[l.split()[0]] + [round(float(x), 7) for x in l.split()[1:4]]
                                        for l in out[b:].split('\n')[:nat]]
        o['positions_final_source'] = 'final coordinates' if fin >= 0 else 'last ionic step (relaxation not converged)'
    if record['input']['calculation'] == 'vc-relax':
        c = list(re.finditer(r'CELL_PARAMETERS \(angstrom\)\n', out))
        o['cell_A_final'] = block3(out, c[-1].start()) if c else None
    record['output'] = o
    return record


def main():
    pool = {}
    for f in sorted(RAW.rglob('pwscf.in')):
        rel = str(f.parent.relative_to(RAW))
        r = parse(f.parent)
        r['archive_path'] = rel
        pool[ident(rel)] = r
    assert len(pool) == 490, len(pool)
    (CACHE/'pool.json').write_text(json.dumps(pool, separators=(',', ':')))
    keep = ('calculation', 'ecutwfc_Ry', 'ecutrho_Ry', 'k_points', 'nspin', 'nat', 'species_counts', 'pseudopotentials')
    derived = {k: {'archive_path': v['archive_path'], **{x: v['input'][x] for x in keep},
                   **{x: v['output'][x] for x in ('last_total_energy_Ry', 'final_energy_Ry', 'relaxation_message',
                                                   'n_scf_cycles', 'total_magnetization_muB_per_cell')},
                   'run': v['run']} for k, v in pool.items()}
    (HERE/'derived.json').write_text(json.dumps(derived, indent=0, sort_keys=True) + '\n')
    print(len(pool), 'runs ->', CACHE/'pool.json', f"({(CACHE/'pool.json').stat().st_size/1e6:.1f} MB)")


if __name__ == '__main__':
    main()
