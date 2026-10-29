import os
import pickle
import numpy as np
from datetime import datetime

import wannierberri as wberri
from wannierberri.symmetry import Rotation, Mirror

prefix = 'InMnTi2'
fermi = 15.684
deltae = 3
step = 0.001
numstep = int((2*deltae)//step + 1)
# Requires .chk, .eig and .mmn files
system=wberri.System_w90(prefix ,berry=True)

generators = [
    # 'Inversion',
    'TimeReversal',
    'C2z',
    'C2y',
    Rotation(3, [1,1,1]),
    Mirror([1,-1,0])
    ]
system.set_symmetry(generators)

grid=wberri.Grid(system,length=300)
print(grid)

parallel=wberri.Parallel(
        method="ray",
        #cluster=True,
        num_cpus=32
        )
print(parallel)

"""
wberri.integrate(
    system, grid,
    Efermi=np.linspace(fermi-deltae, fermi+deltae, numstep),
    smearEf=10, # 10K
    quantities=["ahc","dos","cumdos"],
    parameters = {
        # 'tetra':True
        'spin':True,
        },
    parallel = parallel,
    adpt_num_iter=10,
    fout_name="InMnTi2"
    )
# """

# """
wberri.integrate(
    system, grid,
    omega = np.linspace(0., 7., 701),
    Efermi=np.linspace(fermi-deltae, fermi+deltae, numstep),
    smearEf=50, # ~5meV
    quantities=['opt_conductivity'],
    parameters = {
        # 'tetra':True
        'smr_fixed_width': 0.005,
        'smr_type':'Gaussian',
        #'spin':True,
        },
    parallel = parallel,
    adpt_num_iter=10,
    fout_name="InMnTi2"
    )
# """

"""
tab_result = wberri.tabulate(
    system, grid,
    quantities=[
        "berry", 
        # "spin"
        ],
    frmsf_name="InMnTi2",
    parallel = parallel,
    ibands=np.arange(15,23),
    Ef0=fermi
    )

with open(f'{prefix}_berry.pickle', 'wb') as f:
      pickle.dump(tab_result, f)

# """

