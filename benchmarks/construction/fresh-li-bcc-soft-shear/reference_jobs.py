"""Write the reference-calculation lists (requests in the lab's own schema) run by run_batch.sh.

    python3 -B reference_jobs.py      -> jobs_ref.json (all reference runs, in the order they were made)
"""
import json
from pathlib import Path

A0, AC, AX = 3.4348, 3.30, 3.477      # PBE zero-pressure a0 (k48, 0.01 Ry EOS), compressed cell, experiment (78 K)


def vc(d):
    """Volume-conserving tetragonal partner strain: (1 + d)^2 (1 + e3) = 1."""
    return round(1 / (1 + d) ** 2 - 1, 8)


def req(a, k, s, voigt=None, ecut=None):
    r = {'lattice_A': a, 'kmesh': k, 'smearing_Ry': s}
    if voigt is not None:
        r['strain_voigt'] = voigt
    if ecut:
        r['ecutwfc_Ry'] = ecut
    return r


def cub(a, k, s, ecut=None):
    return [f'c_a{a}_k{k}_s{s}' + (f'_e{ecut}' if ecut else ''), req(a, k, s, ecut=ecut)]


def tet(a, k, s, d, ecut=None):
    return [f't_a{a}_k{k}_s{s}_d{d:+g}' + (f'_e{ecut}' if ecut else ''), req(a, k, s, [d, d, vc(d), 0, 0, 0], ecut)]


def tetlin(a, k, s, d):
    """Linear (not volume-conserving) tetragonal strain (d, d, -2d)."""
    return [f'tl_a{a}_k{k}_s{s}_d{d:+g}', req(a, k, s, [d, d, -2 * d, 0, 0, 0])]


def tri(a, k, s, e):
    return [f'r_a{a}_k{k}_s{s}_e{e:+g}', req(a, k, s, [0, 0, 0, e, e, e])]


def uni(a, k, s, e):
    return [f'u_a{a}_k{k}_s{s}_e{e:+g}', req(a, k, s, [0, 0, e, 0, 0, 0])]


def pair(f, a, k, s, d):
    return [f(a, k, s, d), f(a, k, s, -d)]


def batches():
    b = {}
    jobs = []
    for k, s in [(16, 0.02), (24, 0.02), (32, 0.02), (32, 0.01), (40, 0.01), (24, 0.005), (48, 0.02), (48, 0.01),
                 (48, 0.005), (56, 0.01), (56, 0.005)]:
        jobs.append(cub(A0, k, s))
    for k, s in [(16, 0.02), (24, 0.02), (24, 0.005), (32, 0.01), (32, 0.02), (40, 0.01), (48, 0.02), (48, 0.01),
                 (48, 0.005)]:
        jobs += pair(tet, A0, k, s, 0.01)
    for k, s in [(16, 0.02), (24, 0.02), (32, 0.02), (32, 0.01), (40, 0.01), (48, 0.01), (48, 0.005)]:
        jobs += pair(tri, A0, k, s, 0.01)
    for k, s in [(32, 0.02), (40, 0.01), (48, 0.01), (48, 0.005)]:
        jobs.append(cub(AC, k, s))
        jobs += pair(tet, AC, k, s, 0.01)
    for k, s in [(32, 0.02), (48, 0.01)]:
        jobs += pair(tri, AC, k, s, 0.01)
    jobs.append(cub(AX, 48, 0.01))
    jobs += pair(tet, AX, 48, 0.01, 0.01) + pair(tri, AX, 48, 0.01, 0.01)
    for k, s in [(56, 0.01), (56, 0.005)]:
        jobs += pair(tet, A0, k, s, 0.01)
    b['batch1'] = jobs
    jobs = [cub(a, 48, 0.01) for a in (3.39, 3.41, 3.43, 3.45, 3.47, 3.49)]     # reference EOS
    for d in (0.005, 0.02, 0.03):
        jobs += pair(tet, A0, 48, 0.01, d)
    jobs += pair(tri, A0, 48, 0.01, 0.02) + pair(uni, A0, 48, 0.01, 0.01)
    for a in (3.39, 3.41, 3.43, 3.45, 3.47):
        jobs += [cub(a, 32, 0.01), cub(a, 16, 0.02)]
    jobs += [cub(AC, 16, 0.02)] + pair(tet, AC, 16, 0.02, 0.01) + pair(tri, AC, 16, 0.02, 0.01)
    jobs += [cub(AC, 24, 0.005)] + pair(tet, AC, 24, 0.005, 0.01) + pair(tri, AC, 24, 0.005, 0.01)
    jobs += [cub(A0, 16, 0.005)] + pair(tet, A0, 16, 0.005, 0.01) + pair(tri, A0, 16, 0.005, 0.01)
    jobs += pair(tetlin, AC, 48, 0.01, 0.01)
    jobs += [cub(AC, 56, 0.01)] + pair(tet, AC, 56, 0.01, 0.01) + pair(tri, AC, 48, 0.005, 0.01)
    jobs += pair(tri, A0, 56, 0.01, 0.01)
    b['batch2'] = jobs
    jobs = [cub(AC, 64, 0.01)] + pair(tet, AC, 64, 0.01, 0.01) + [cub(AC, 56, 0.005)] + pair(tet, AC, 56, 0.005, 0.01)
    jobs += pair(tri, A0, 56, 0.005, 0.01)
    b['batch3'] = jobs
    jobs = pair(tri, AC, 56, 0.01, 0.01) + pair(tri, AC, 56, 0.005, 0.01)
    jobs += pair(tri, A0, 64, 0.005, 0.01) + pair(tri, A0, 64, 0.01, 0.01)
    b['batch4'] = jobs
    b['batch5'] = pair(tri, AC, 64, 0.01, 0.01) + pair(tri, AC, 64, 0.005, 0.01)
    b['batch6'] = [cub(A0, 40, 0.01, 60), tet(A0, 40, 0.01, 0.01, 60), tet(A0, 40, 0.01, -0.01, 60)]   # cutoff check
    return b


if __name__ == '__main__':
    all_jobs = [j for jobs in batches().values() for j in jobs]
    names = [n for n, _ in all_jobs]
    assert len(names) == len(set(names))
    Path(__file__).with_name('jobs_ref.json').write_text(json.dumps(all_jobs, indent=0) + '\n')
    print(len(all_jobs), 'jobs')
