"""The reference calculations for fresh-h-chain-peierls-vs-afm (requests only; see make_reference.py).

Every request is an ordinary lab request and runs through task_lab.run, exactly as an agent's would.
Tags: R<spacing>_d<delta>_<nm|afm|fm>[_<variant>].
"""
STANDARD = dict(k_points=48, smearing='mv', degauss_Ry=0.005, box_A=10, ecutwfc_Ry=80)
AFM, FM = [0.5, -0.5], [0.5, 0.5]


def request(R, d, mag, atoms=2, **settings):
    bonds = [round(R * (1 - d), 6), round(R * (1 + d), 6)] if atoms == 2 else [R] * atoms
    return {'kind': 'scf', 'bond_lengths_A': bonds, 'starting_magnetization': mag, **{**STANDARD, **settings}}


PLAN = {}


def add(tag, *args, **kw):
    PLAN[tag] = request(*args, **kw)


# Uniform chains, standard settings
for R in (1.0, 1.1, 1.25, 1.5, 1.75, 2.0):
    add(f'R{R}_d0_nm', R, 0, None)
    add(f'R{R}_d0_afm', R, 0, AFM)
for R in (1.25, 1.75, 2.0):
    add(f'R{R}_d0_fm', R, 0, FM)
for R in (1.25, 1.75):
    add(f'R{R}_d0_fm_1atom', R, 0, [0.5], atoms=1, k_points=96)
# Uniform chains, converged settings (weak-coupling SDW at 1.25; magnetic energies at 1.75 and 2.0)
add('R1.25_d0_afm_mv002_k192', 1.25, 0, AFM, degauss_Ry=0.002, k_points=192)
add('R1.25_d0_nm_mv002_k192', 1.25, 0, None, degauss_Ry=0.002, k_points=192)
add('R1.25_d0_afm_mv002_k96', 1.25, 0, AFM, degauss_Ry=0.002, k_points=96)
add('R1.25_d0_nm_mv002_k96', 1.25, 0, None, degauss_Ry=0.002, k_points=96)
add('R1.25_d0_afm_mv005_k96', 1.25, 0, AFM, k_points=96)
add('R1.25_d0_afm_g005_k96', 1.25, 0, AFM, k_points=96, smearing='gaussian')
add('R1.25_d0_afm_fd005_k96', 1.25, 0, AFM, k_points=96, smearing='fd')
add('R1.25_d0_afm_fd002_k96', 1.25, 0, AFM, k_points=96, smearing='fd', degauss_Ry=0.002)
add('R1.25_d0_afm_mv01', 1.25, 0, AFM, degauss_Ry=0.01)
add('R1.25_d0_afm_mv02', 1.25, 0, AFM, degauss_Ry=0.02)
add('R1.25_d0_nm_mv02', 1.25, 0, None, degauss_Ry=0.02)
add('R1.25_d0_afm_small_start_k96', 1.25, 0, [0.1, -0.1], k_points=96)
for R in (1.75, 2.0):
    add(f'R{R}_d0_nm_mv002_k96', R, 0, None, degauss_Ry=0.002, k_points=96)
    add(f'R{R}_d0_afm_mv002_k96', R, 0, AFM, degauss_Ry=0.002, k_points=96)
# Equilibrium dimerization (nonmagnetic; an AFM start at the minimum checks that it stays nonmagnetic)
for R, deltas in ((1.0, (0.2, 0.25, 0.3)), (1.25, (0.35, 0.39, 0.4, 0.41, 0.45)),
                  (1.75, (0.5, 0.55, 0.56, 0.57, 0.58, 0.6)), (2.0, (0.6, 0.625, 0.65))):
    for d in deltas:
        add(f'R{R}_d{d}_nm', R, d, None)
for R, d in ((1.0, 0.25), (1.25, 0.4), (1.75, 0.57), (2.0, 0.625)):
    add(f'R{R}_d{d}_afm', R, d, AFM)
# Coexistence boundary: AFM starts in dimerized chains
for R, deltas in ((1.5, (0.02, 0.03, 0.04, 0.05)), (1.75, (0.05, 0.1, 0.105, 0.11, 0.115, 0.12, 0.14, 0.16, 0.2, 0.3)),
                  (2.0, (0.05, 0.1, 0.15, 0.19, 0.2, 0.205, 0.21, 0.22, 0.3))):
    for d in deltas:
        add(f'R{R}_d{d}_afm', R, d, AFM)
# Settings robustness of the boundary and of the dimerization energy
add('R1.75_d0.11_afm_mv002_k96', 1.75, 0.11, AFM, degauss_Ry=0.002, k_points=96)
add('R2.0_d0.2_afm_mv002_k96', 2.0, 0.2, AFM, degauss_Ry=0.002, k_points=96)
add('R1.75_d0.11_afm_start1', 1.75, 0.11, [1.0, -1.0])
add('R1.75_d0.57_nm_ec60', 1.75, 0.57, None, ecutwfc_Ry=60)
add('R1.75_d0_afm_ec60', 1.75, 0, AFM, ecutwfc_Ry=60)
add('R1.75_d0.57_nm_ec100', 1.75, 0.57, None, ecutwfc_Ry=100)
add('R1.75_d0_afm_ec100', 1.75, 0, AFM, ecutwfc_Ry=100)
add('R1.75_d0.57_nm_box8', 1.75, 0.57, None, box_A=8)
add('R1.75_d0_afm_box8', 1.75, 0, AFM, box_A=8)
add('R1.75_d0.57_nm_box12', 1.75, 0.57, None, box_A=12)
add('R1.75_d0_afm_box12', 1.75, 0, AFM, box_A=12)
# 4-atom cells at R = 1.25 A (same k density as 2 atoms at 192 k): the SDW again; an up-up-down-down start relaxes to NM
PLAN['R1.25_4atom_afm_mv002_k96'] = {'kind': 'scf', 'bond_lengths_A': [1.25] * 4, 'starting_magnetization': [0.5, -0.5, 0.5, -0.5],
                                     **dict(STANDARD, degauss_Ry=0.002, k_points=96)}
PLAN['R1.25_4atom_uudd_mv002_k96'] = {'kind': 'scf', 'bond_lengths_A': [1.25] * 4, 'starting_magnetization': [0.5, 0.5, -0.5, -0.5],
                                      **dict(STANDARD, degauss_Ry=0.002, k_points=96)}
