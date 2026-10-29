"""Parse the released Calderon et al. (2023) archive into pool.json (everything the lab can serve).

    python3 -B extract.py [RAW_ROOT]

RAW_ROOT defaults to ~/.cache/ai-theorist/wurtzite-ferroelectric-switching/raw/AlBN-system-Science
(created by fetch.sh). Nothing is recomputed: values are read from the authors' QE outputs.
"""
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path.home()/'.cache/ai-theorist/wurtzite-ferroelectric-switching/raw/AlBN-system-Science'
SYSTEMS = ['AlN', 'AlBN']
RUNS = ['NEB', 'refined_NEB', 'refined_NEB2']
RY = 13.605693122994   # eV (QE 6.7 value)


def read_axsf(path):
    """Final image structures of a QE NEB run (prefix.axsf, written at the last iteration)."""
    lines = path.read_text().split('\n')
    n_img = int(lines[0].split()[1])
    i = [k for k, l in enumerate(lines) if l.strip() == 'PRIMVEC'][0]
    cell = [[float(x) for x in lines[i + j].split()] for j in (1, 2, 3)]
    images = []
    k = 0
    while len(images) < n_img:
        while not lines[k].strip().startswith('PRIMCOORD'):
            k += 1
        nat = int(lines[k + 1].split()[0])
        species, pos, forces = [], [], []
        for line in lines[k + 2:k + 2 + nat]:
            p = line.split()
            species.append(p[0])
            pos.append([float(x) for x in p[1:4]])
            forces.append([float(x) for x in p[4:7]])
        images.append({'species': species, 'positions_A': pos, 'forces_xsf': forces})
        k += 2 + nat
    return cell, images


def parse_neb_out(path):
    text = path.read_text(errors='ignore')
    settings = {}
    for key in ['restart_mode', 'opt_scheme', 'num_of_images', 'nstep_path', 'CI_scheme', 'first_last_opt',
                'use_freezing', 'ds', 'k_max', 'k_min', 'path_thr', 'initial path length']:
        m = re.findall(r'^\s*' + re.escape(key) + r'\s*=\s*(.+?)\s*$', text, re.M)
        if m:
            settings[key] = m[-1]
    blocks = list(re.finditer(r'image\s+energy \(eV\)\s+error \(eV/A\)\s+frozen\s*\n\s*\n((?:\s+\d+\s+-?\d+\.\d+\s+\d+\.\d+\s+[TF]\s*\n)+)', text))
    last = blocks[-1].group(1)
    table = []
    for line in last.strip().split('\n'):
        p = line.split()
        table.append({'image': int(p[0]), 'energy_eV': float(p[1]), 'error_eV_per_A': float(p[2]), 'frozen': p[3]})
    act = re.findall(r'activation energy \(->\) =\s+(-?\d+\.\d+) eV', text)
    iters = re.findall(r'-+ iteration\s+(\d+) -+', text)
    conv = re.findall(r'neb: convergence achieved in\s+\d+ iterations', text)
    tail = [l for l in text.rstrip().split('\n')[-40:]]
    return {'settings': settings, 'final_table': table, 'n_iteration_blocks': len(blocks),
            'last_iteration_label': int(iters[-1]) if iters else None,
            'activation_energy_forward_history_eV': [float(a) for a in act[-6:]],
            'convergence_message': conv[-1] if conv else None, 'tail': tail}


def engine_settings(neb_in):
    text = neb_in.read_text()
    out = {}
    for key in ['ecutwfc', 'ecutrho', 'occupations', 'nat', 'ntyp', 'nosym', 'conv_thr', 'mixing_beta']:
        m = re.search(r'^\s*' + key + r'\s*=\s*(\S+)', text, re.M)
        if m:
            out[key] = m.group(1).strip(',')
    m = re.search(r'K_POINTS automatic\s*\n\s*(.+)', text)
    out['K_POINTS automatic'] = m.group(1).strip()
    out['ATOMIC_SPECIES'] = re.findall(r'^\s*([A-Z][a-z]?)\s+[\d.]+\s+(\S+\.upf)', text, re.M)
    return out


def parse_pw_in(path):
    text = path.read_text()
    m = re.search(r'CELL_PARAMETERS angstrom\s*\n(.+)\n(.+)\n(.+)', text)
    cell = [[float(x) for x in m.group(i).split()] for i in (1, 2, 3)]
    i = text.index('ATOMIC_POSITIONS angstrom')
    species, pos = [], []
    for line in text[i:].split('\n')[1:]:
        p = line.split()
        if len(p) < 4:
            break
        species.append(p[0])
        pos.append([float(x) for x in p[1:4]])
    settings = {}
    for key in ['calculation', 'lberry', 'gdir', 'nppstr', 'ecutwfc', 'ecutrho', 'occupations', 'conv_thr']:
        mm = re.search(r'^\s*' + key + r'\s*=\s*(\S+)', text, re.M)
        if mm:
            settings[key] = mm.group(1)
    settings['K_POINTS automatic'] = re.search(r'K_POINTS automatic\s*\n\s*(.+)', text).group(1).strip()
    return cell, species, pos, settings


def parse_berry(out_path):
    text = out_path.read_text(errors='ignore')
    i = text.index('POLARIZATION CALCULATION')
    j = text.index('The polarization direction is', i)
    block = text[i - 70:j + 80]
    block = block[block.index('\n') + 1:] if '\n' in block[:70] else block
    m = re.search(r'P =\s+(-?\d+\.\d+)\s+\(mod\s+(\d+\.\d+)\)\s+C/m\^2', block)
    ion = re.search(r'Ionic Phase:\s+(-?\d+\.\d+)', block)
    ele = re.search(r'Electronic Phase:\s+(-?\d+\.\d+)', block)
    tot = re.search(r'TOTAL PHASE:\s+(-?\d+\.\d+)\s+MOD_TOT:\s+(\d+)', block)
    return {'text': block.rstrip(), 'P_C_m2': float(m.group(1)), 'mod_C_m2': float(m.group(2)),
            'ionic_phase': float(ion.group(1)), 'electronic_phase': float(ele.group(1)),
            'total_phase': float(tot.group(1)), 'mod_tot': int(tot.group(2))}


def scf_energy(out_path):
    text = out_path.read_text(errors='ignore')
    e = re.findall(r'^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry', text, re.M)
    return float(e[-1]), 'convergence has been achieved' in text


def parse_vc(path):
    text = path.read_text(errors='ignore')
    i = text.rindex('Begin final coordinates')
    j = text.index('End final coordinates', i)
    final = text[i:j + len('End final coordinates')]
    e = re.findall(r'^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry', text, re.M)
    p = re.findall(r'P=\s*(-?\d+\.\d+)', text)
    nat = int(re.search(r'number of atoms/cell\s+=\s+(\d+)', text).group(1))
    conv = re.search(r'bfgs converged in.*', text)
    return {'final_coordinates_block': final, 'nat': nat, 'final_total_energy_Ry': float(e[-1]),
            'final_pressure_kbar': float(p[-1]), 'bfgs': conv.group(0).strip() if conv else None}


def main():
    pool = {'neb_runs': {}, 'berry': {}, 'vc_relax': {}, 'born_charges': {}}
    for s in SYSTEMS:
        pool['neb_runs'][s] = {}
        for r in RUNS:
            d = RAW/s/'NEB'/r
            info = parse_neb_out(d/'AlBN.neb.out')
            if s == 'AlBN' and r == 'refined_NEB2':
                # The restart log (AlBN.neb.out) continues the earlier 75-image run (AlBN.neb1.out).
                info['earlier_log_activation_history_eV'] = parse_neb_out(d/'AlBN.neb1.out')['activation_energy_forward_history_eV']
            info['engine'] = engine_settings(d/'neb.in')
            info['scheduler_out'] = (d/'scheduler.out').read_text(errors='ignore').strip().split('\n')[-6:]
            cell, images = read_axsf(d/'AlBN.axsf')
            info['cell_A'] = cell
            info['images'] = images
            assert len(images) == len(info['final_table']), (s, r)
            pool['neb_runs'][s][r] = info
        pool['berry'][s] = {}
        for d in sorted((RAW/s/'polarization/nscf').glob('*_image')):
            idx = int(d.name.split('_')[0])
            cell, species, pos, settings = parse_pw_in(d/'nscf.in')
            b = parse_berry(next(d.glob('*.nscf.out')))
            e, conv = scf_energy(next((RAW/s/'polarization/scf'/d.name).glob('*.scf.out')))
            _, _, _, scf_settings = parse_pw_in(RAW/s/'polarization/scf'/d.name/'scf.in')
            pool['berry'][s][idx] = {'cell_A': cell, 'species': species, 'positions_A': pos, 'nscf_settings': settings,
                                     'scf_settings': scf_settings, 'scf_total_energy_Ry': e, 'scf_converged': conv, **b}
    pool['vc_relax'] = {'AlN_primitive': parse_vc(RAW/'AlN/vc-relax/AlBN.vc-relax.out'),
                        'AlN_32atom_initial_state': parse_vc(RAW/'AlN/NEB/initial_state/AlBN.vc-relax.out'),
                        'AlBN_32atom_initial_state': parse_vc(RAW/'AlBN/NEB/initial_state/AlBN.vc-relax.out'),
                        'AlBN_192atom': parse_vc(RAW/'AlBN/443_vc-relax/AlBN.vc-relax.out')}
    with open(RAW/'Polarization and Coercive Field/Born effective charge (AlBN).csv') as f:
        pool['born_charges']['AlBN'] = [{'atom': row[0], 'Z33_e': float(row[1])} for row in list(csv.reader(f))[1:]]
    # Consistency checks between the released files
    for s in SYSTEMS:
        neb = pool['neb_runs'][s]['refined_NEB2']
        assert len(neb['images']) == len(pool['berry'][s])
        for k, img in enumerate(neb['images']):
            b = pool['berry'][s][k]
            a, c = np.array(img['positions_A']), np.array(b['positions_A'])
            for sp in set(img['species']):   # same atoms, possibly listed in another order
                ia = [i for i, x in enumerate(img['species']) if x == sp]
                ic = [i for i, x in enumerate(b['species']) if x == sp]
                dist = np.linalg.norm(a[ia][:, None] - c[ic][None], axis=2)
                assert len(ia) == len(ic) and dist.min(axis=1).max() < 2e-5, (s, k, sp)
            d_e = b["scf_total_energy_Ry"] * RY - neb["final_table"][k]["energy_eV"]
            if abs(d_e) > 2e-3:
                print("scf vs NEB energy mismatch", s, k, round(d_e, 4))
        auth = np.loadtxt(RAW/f'Polarization and Coercive Field/{s}_polarization.txt', skiprows=1)
        raw = np.array([pool['berry'][s][k]['P_C_m2'] for k in range(len(neb['images']))])
        assert np.allclose(auth[:, 1], raw, atol=1e-6), s
    out = HERE/'pool.json'
    out.write_text(json.dumps(pool) + '\n')
    print('wrote', out, round(out.stat().st_size / 2**20, 2), 'MB')


if __name__ == '__main__':
    main()
