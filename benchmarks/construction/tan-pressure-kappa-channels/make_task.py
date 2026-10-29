"""Build the lab data, catalog, private references and ablation values for tan-pressure-kappa-channels.

    docker run --rm --cpuset-cpus <cpus> --user $(id -u):$(id -g) -e HOME=/tmp -e OPENBLAS_NUM_THREADS=1 \
        -v $REPO:/repo -v $CACHE/raw/IFCs-data-pressure-TaN:/raw -v $CACHE/fourph:/s \
        -w /repo/benchmarks/construction/tan-pressure-kappa-channels \
        ai-theorist-benchmark-runtime:20261007 python3 -B make_task.py [stage]

stage "data" writes environment/lab/data and the agent catalog; stage "reference" (default: both) computes the
hidden references, the oracle-route values, robustness variants and ablations with the lab's own code path
(tanbte.py) and writes private/hidden.json, hidden_values.json, solution/ablation_values.json and calibration.json.
"""
import json
import shutil
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import tanbte as tb     # noqa: E402

REPO = HERE.parents[2]
TASK = REPO/'benchmarks/tasks/tan-pressure-kappa-channels'
LAB = TASK/'environment/lab'
RAW = Path('/raw')
FOURPH = Path('/s')
PRESSURES = [0, 20, 40, 60, 80, 100]
MESHES = [8, 10, 12, 14, 16, 18, 20]
PHEL_RELEASED = [0, 20, 40, 60, 100]          # BTE.sr_phel_metal withheld at 80 GPa
FOURPH_RELEASED = [(0, 300), (20, 300), (60, 300), (100, 300), (20, 600), (100, 600)]
FOURPH_RUNS = {(0, 300): [1, 2], (20, 300): [1], (40, 300): [1], (60, 300): [1], (80, 300): [1, 2], (100, 300): [1],
               (20, 600): [1], (60, 600): [1], (100, 600): [1]}
BUDGET = 30
COSTS = {'archived_file': 1, 'harmonic': 1, 'rates': 2, 'kappa': 1}


def fourph_rates(p, T):
    runs = [np.load(FOURPH/('r4_%d_16_0.1_8192_%d_%d.npz' % (p, s, T))) for s in FOURPH_RUNS[(p, T)]]
    return {k: np.mean([r[k] for r in runs], axis=0) for k in ('pp', 'pm', 'mm')}


def crystal(p, with_fc3=True):
    pc = tb.read_poscar(RAW/str(p)/'POSCAR')
    fc2 = tb.read_fc2(RAW/str(p)/'FORCE_CONSTANTS_2ND', 2, (5, 5, 5))
    fc3 = tb.read_fcn(RAW/str(p)/'FORCE_CONSTANTS_3RD', 3) if with_fc3 else None
    return pc, fc2, fc3


def stage_data():
    data = LAB/'data'
    if data.exists():
        shutil.rmtree(data)
    (data/'archive').mkdir(parents=True)
    shutil.copy(HERE/'tanbte.py', LAB/'tanbte.py')
    archived = {}
    for p in PRESSURES:
        pc, fc2, fc3 = crystal(p)
        np.savez_compressed(data/('P%d.npz' % p), A=pc['A'], frac=pc['frac'], fc2=fc2,
                            fc3_R=fc3['R'], fc3_atoms=fc3['atoms'], fc3_fc=fc3['fc'])
        (data/'archive'/str(p)).mkdir()
        files = ['POSCAR', 'CONTROL'] + (['BTE.sr_phel_metal'] if p in PHEL_RELEASED else [])
        for f in files:
            shutil.copy(RAW/str(p)/f, data/'archive'/str(p)/f)
        shutil.copy(RAW/'README.md', data/'archive'/str(p)/'README.md')
        archived[str(p)] = files + ['README.md']
    for p, T in FOURPH_RELEASED:
        np.savez_compressed(data/('fourph_P%d_T%d.npz' % (p, T)), **fourph_rates(p, T))
    cr = tb.Crystal(*crystal(0)[:2])
    g = tb.Grid(cr, (16, 16, 16))
    meta = {'pressures_GPa': PRESSURES, 'meshes': MESHES, 'archived_files': archived,
            'fourph_available': [list(x) for x in FOURPH_RELEASED], 'nirr16': int(len(g.irr))}
    (data/'meta.json').write_text(json.dumps(meta, indent=1) + '\n')
    requests = {
        'archived_file': {'kind': 'archived_file', 'pressure_GPa': 0, 'file': 'CONTROL', 'cost': COSTS['archived_file']},
        'harmonic': {'kind': 'harmonic', 'pressure_GPa': 0, 'mesh': 16, 'cost': COSTS['harmonic']},
        'rates': {'kind': 'rates', 'pressure_GPa': 0, 'mesh': 16, 'T_K': 300, 'scalebroad': 0.1, 'cost': COSTS['rates']},
        'kappa': {'kind': 'kappa', 'pressure_GPa': 0, 'mesh': 16, 'T_K': 300, 'scalebroad': 0.1,
                  'channels': ['3ph', 'isotope'], 'solver': 'iterative', 'extra_rates_ps': '<optional 270 x 6 list>',
                  'cost': '%d (mesh <= 16) or %d (mesh 18, 20)' % (COSTS['kappa'], 2 * COSTS['kappa'])}}
    public = {'task_id': 'tan-pressure-kappa-channels', 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
              'requests': {k: {'cost': COSTS[k]} for k in COSTS}, 'catalog': '/app/catalog.json'}
    (LAB/'public.json').write_text(json.dumps(public, indent=1) + '\n')
    tests = json.loads((TASK/'tests/public.json').read_text())
    tests['solver_budget'] = BUDGET
    (TASK/'tests/public.json').write_text(json.dumps(tests, indent=1) + '\n')
    structures = {}
    for p in PRESSURES:
        pc = tb.read_poscar(RAW/str(p)/'POSCAR')
        structures[str(p)] = {'a_A': round(float(np.linalg.norm(pc['A'][0])), 6),
                              'c_A': round(float(np.linalg.norm(pc['A'][2])), 6)}
    catalog = {
        'material': 'theta-TaN (WC-type, P-6m2, 2 atoms/cell), DFT force constants released by Kundu et al. '
                    '(PRL 132, 116301 (2024); Zenodo record 10564835, CC BY 4.0)',
        'pressures_GPa': PRESSURES,
        'lattice_parameters': structures,
        'archive_contents': {
            'per pressure (as archived)': ['POSCAR', 'CONTROL (ShengBTE/FourPhonon input)', 'FORCE_CONSTANTS_2ND',
                                           'FORCE_CONSTANTS_3RD', 'FORCE_CONSTANTS_4TH', 'BTE.sr_phel_metal', 'README.md'],
            'served as files (archived_file)': archived,
            'used by the lab but not served as files': 'FORCE_CONSTANTS_2ND and _3RD at every pressure; '
                                                       'FORCE_CONSTANTS_4TH where four-phonon lookups exist',
            'withheld': {'FORCE_CONSTANTS_4TH': [40, 80], 'BTE.sr_phel_metal': [80],
                         'four-phonon lookups': 'everything except the (GPa, K) pairs listed in fourph_available'}},
        'fourph_available_GPa_K': [list(x) for x in FOURPH_RELEASED],
        'meshes': MESHES,
        'grid_16': {
            'convention': 'Gamma-centred 16x16x16 grid, q = (i, j, k)/16 in reduced reciprocal coordinates; grid '
                          'index n = i + 16 j + 256 k',
            'irreducible_q_reduced': np.round(g.q[g.irr], 6).tolist(),
            'multiplicity': g.mult.tolist(),
            'grid_to_irreducible': g.star_of.tolist(),
            'note': 'extra_rates_ps and all per-mode outputs use this irreducible order; branches are sorted by '
                    'frequency at each q (6 branches)'},
        'requests': requests,
        'budget': BUDGET,
        'units': {'omega': 'rad/ps', 'rates': '1/tau in ps^-1', 'kappa': 'W m^-1 K^-1', 'velocity': 'km/s'},
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    print('data stage written; nirr16 =', len(g.irr))


if __name__ == '__main__':
    stages = sys.argv[1:] or ['data', 'reference']
    if 'data' in stages:
        stage_data()
    if 'reference' in stages:
        import reference
        reference.main()
