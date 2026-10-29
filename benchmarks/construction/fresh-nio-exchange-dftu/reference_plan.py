"""The reference calculations for fresh-nio-exchange-dftu (requests only; run by run_reference.py).

Every request is an ordinary lab request and runs through task_lab.run, exactly as an agent's would.
Tags: U<U>_<cell>_<state>_k<mesh>[_<variant>]; cells pr (primitive_1fu), rh (rhombohedral_2fu),
te (tetragonal_2fu); states fm, af1 (AFM-I, [001] planes), af2 (AFM-II, [111] planes).
"""
CELL = {'pr': 'primitive_1fu', 'rh': 'rhombohedral_2fu', 'te': 'tetragonal_2fu', 'cu': 'cubic_4fu'}
START = {('pr', 'fm'): [0.5], ('rh', 'fm'): [0.5, 0.5], ('rh', 'af2'): [0.5, -0.5],
         ('te', 'fm'): [0.5, 0.5], ('te', 'af1'): [0.5, -0.5]}
STANDARD = dict(kshift=0, ecutwfc_Ry=50, smearing='gaussian', degauss_Ry=0.005)
PLAN = {}


def add(U, cell, state, k, variant='', mag=None, **settings):
    mesh = [k] * 3 if isinstance(k, int) else list(k)
    tag = f'U{U}_{cell}_{state}_k{"".join(map(str, mesh))}' + (f'_{variant}' if variant else '')
    PLAN[tag] = {'kind': 'scf', 'cell': CELL[cell], 'starting_magnetization': mag or START[(cell, state)],
                 'hubbard_U_eV': U, 'kpoints': mesh, **{**STANDARD, **settings}}


# Exchange at U = 4 and 8 eV: k ladders in the cell of each pair, cutoff checks, cross-cell FM references
for U in (4.0, 8.0):
    for k in (3, 4, 6, 7):
        add(U, 'rh', 'fm', k)
        add(U, 'rh', 'af2', k)
    for mesh in ((4, 4, 3), (4, 4, 4), (6, 6, 4), (6, 6, 6), (8, 8, 6), (8, 8, 8), (10, 10, 8)):
        add(U, 'te', 'fm', mesh)
        add(U, 'te', 'af1', mesh)
    for k in (8, 12):
        add(U, 'pr', 'fm', k)
    for ec in (40, 75):
        add(U, 'rh', 'fm', 4, f'e{ec}', ecutwfc_Ry=ec)
        add(U, 'rh', 'af2', 4, f'e{ec}', ecutwfc_Ry=ec)
        add(U, 'te', 'fm', (8, 8, 6), f'e{ec}', ecutwfc_Ry=ec)
        add(U, 'te', 'af1', (8, 8, 6), f'e{ec}', ecutwfc_Ry=ec)
    add(U, 'rh', 'fm', 4, 'm1', mag=[1.0, 1.0])
    add(U, 'rh', 'af2', 4, 'm1', mag=[1.0, -1.0])
    add(U, 'te', 'af1', (8, 8, 6), 'm01', mag=[0.1, -0.1])
# Smearing at U = 4 eV (FM gap 0.77 eV)
for sm, dg in (('mv', 0.02), ('fd', 0.02), ('gaussian', 0.05), ('fd', 0.05)):
    add(4.0, 'rh', 'fm', 6, f'{sm}{dg}', smearing=sm, degauss_Ry=dg)
    add(4.0, 'rh', 'af2', 6, f'{sm}{dg}', smearing=sm, degauss_Ry=dg)
    add(4.0, 'te', 'fm', (8, 8, 6), f'{sm}{dg}', smearing=sm, degauss_Ry=dg)
    add(4.0, 'te', 'af1', (8, 8, 6), f'{sm}{dg}', smearing=sm, degauss_Ry=dg)
# Plain PBE (U = 0): AFM-II insulator, metallic FM (k, smearing and starting-moment sensitivity)
for k in (4, 6, 7):
    add(0.0, 'rh', 'af2', k)
for k in (4, 6):
    add(0.0, 'rh', 'fm', k)
add(0.0, 'rh', 'fm', 6, 'm01', mag=[0.1, 0.1])
add(0.0, 'rh', 'af2', 6, 'm01', mag=[0.1, -0.1])
for k in (8, 12, 16):
    add(0.0, 'pr', 'fm', k)
    add(0.0, 'pr', 'fm', k, 'g002', degauss_Ry=0.002)
for sm, dg in (('mv', 0.02), ('mv', 0.005), ('fd', 0.02), ('gaussian', 0.02)):
    add(0.0, 'pr', 'fm', 12, f'{sm}{dg}', smearing=sm, degauss_Ry=dg)
add(0.0, 'pr', 'fm', 16, 'm1', mag=[1.0])
add(0.0, 'pr', 'fm', 12, 'm1', mag=[1.0])
add(0.0, 'te', 'fm', (8, 8, 6))
add(0.0, 'te', 'af1', (8, 8, 6))
# Where the FM state becomes insulating
for U in (1.0, 1.5, 1.8, 2.0, 2.1, 2.2, 2.3, 2.5, 3.0, 3.5):
    add(U, 'pr', 'fm', 12)
for U in (2.0, 2.1, 2.2):
    add(U, 'pr', 'fm', 16)
for U in (2.0, 3.0, 6.0):
    add(U, 'pr', 'fm', 8)
add(2.0, 'rh', 'fm', 6)
add(2.0, 'rh', 'af2', 6)
# Shifted meshes miss the FM band extrema near the gap closing (gap on the mesh too large)
for U in (2.0, 2.2, 2.5):
    add(U, 'pr', 'fm', 12, 's1', kshift=1)
add(2.2, 'pr', 'fm', 16, 's1', kshift=1)
