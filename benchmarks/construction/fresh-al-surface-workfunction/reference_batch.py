"""Write reference_batch.json: every construction request used by make_task.py (name -> lab request).

    python3 -B reference_batch.py

Slab requests are exactly what an agent could send to the lab, except the reference-only 13-, 14- and 18-layer
slabs on the 18x18 mesh (beyond the lab's 12-layer limit) and the lattice-constant scan ('lattice_A' is a
construction-only override that validate() rejects). Names starting with dflt_ omit nbnd (pw.x default bands).
"""
import json
import math
from pathlib import Path


def nb(n):
    return int(math.ceil(2.5 * n)) + 4


def slab(face, n, k=18, vac=15.0, sm=0.02, relax='none', bands=True):
    req = {'kind': 'slab', 'face': face, 'layers': n, 'vacuum_A': vac, 'kpoints': k, 'smearing_Ry': sm,
           'ecutwfc_Ry': 30, 'relax': relax}
    if bands:
        req['nbnd'] = nb(n)
    return req


def bulk(cell, k, face=None, kz=None, sm=0.02, a=None):
    req = {'kind': 'bulk', 'cell': cell, 'kpoints': k, 'smearing_Ry': sm, 'ecutwfc_Ry': 30}
    if face:
        req.update(face=face, kz=kz)
    if a:
        req['lattice_A'] = a
    return req


def batch():
    b = {}
    for a in (4.00, 4.02, 4.03, 4.04, 4.05, 4.06, 4.08):
        b[f'eos_a{a:.3f}'] = bulk('primitive', 24, a=a)
    for k in (12, 16, 18, 20, 24, 28, 32, 36, 40):
        b[f'bulk_prim_k{k}'] = bulk('primitive', k)
    for kz in (6, 12, 18, 24):
        b[f'bulk_or111_k18_kz{kz}'] = bulk('oriented', 18, '111', kz)
    for k in (12, 24, 30, 36):
        b[f'bulk_or111_k{k}_kz12'] = bulk('oriented', k, '111', 12)
    for k in (18, 24):
        for sm in (0.01, 0.04):
            b[f'bulk_or111_k{k}_kz12_sm{int(sm * 100):02d}'] = bulk('oriented', k, '111', 12, sm)
    for kz in (9, 18, 27):
        b[f'bulk_or100_k18_kz{kz}'] = bulk('oriented', 18, '100', kz)
    for k in (12, 24, 30):
        b[f'bulk_or100_k{k}_kz18'] = bulk('oriented', k, '100', 18)
    for face in ('111', '100'):
        for n in range(3, 15):
            b[f's{face}_N{n}_k18'] = slab(face, n)
        b[f'h{face}_N18_k18'] = slab(face, 18)
        for k in (12, 24, 30) + ((36,) if face == '111' else ()):
            b[f's{face}_N7_k{k}'] = slab(face, 7, k)
        for n in range(8, 13):
            b[f's{face}_N{n}_k12'] = slab(face, n, 12)
        for n in range(8, 13 if face == '111' else 12):
            b[f's{face}_N{n}_k24'] = slab(face, n, 24)
    for sm in (0.01, 0.04):
        b[f's111_N7_k18_sm{int(sm * 100):02d}'] = slab('111', 7, sm=sm)
        b[f's111_N10_k18_sm{int(sm * 100):02d}'] = slab('111', 10, sm=sm)
    for v in (8.0, 11.0, 20.0):
        b[f's111_N7_k18_v{int(v)}'] = slab('111', 7, vac=v)
    for n in (7, 9, 10, 11, 12):
        b[f'r111_N{n}_k18_both'] = slab('111', n, relax='outer_two_both_sides')
    for n in (7, 9, 10, 11):
        b[f'r100_N{n}_k18_both'] = slab('100', n, relax='outer_two_both_sides')
    for n in range(3, 13):
        b[f'dflt_s111_N{n}_k18'] = slab('111', n, bands=False)
    for n in range(7, 13):
        b[f'dflt2_s100_N{n}_k18'] = slab('100', n, bands=False)
    for n, bands in ((7, 18), (7, 22), (7, 26), (10, 24), (10, 28), (10, 32)):
        b[f'nb_111_N{n}_k18_nb{bands}'] = dict(slab('111', n), nbnd=bands)
    for bands in (14, 18):
        b[f'nb_100_N5_k18_nb{bands}'] = dict(slab('100', 5), nbnd=bands)
    b['nb_100_N5_k18_dflt'] = slab('100', 5, bands=False)
    return [{'name': k, 'request': v} for k, v in b.items()]


if __name__ == '__main__':
    out = Path(__file__).with_name('reference_batch.json')
    out.write_text(json.dumps(batch(), indent=1) + '\n')
    print(len(batch()), 'requests ->', out)
