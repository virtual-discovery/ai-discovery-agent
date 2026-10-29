"""Plausible wrong paths for dense-hydrogen-llt-crossover; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/dense-hydrogen-llt-crossover/make_task.py from the served (released) data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
NQ = 'nqe_crossover_pressure_shift_GPa'


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


def merged(field, part):
    def mutate(payload, trace):
        payload[field] = {**payload[field], **part}
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'widom_from_whole_runs', setting(widom_rs_N2048=A['widom_whole_run'])),
    ('dp2', 'classification_from_small_cells', setting(classification=A['classification_small_cells'])),
    ('dp3', 'molecular_side_liquid_because_it_diffuses', setting(low_density_phase_900K=A['phase_liquid_from_diffusion'])),
    ('dp4', 'absolute_structure_factor_threshold', setting(classification=A['classification_absolute_structure_factor'])),
    ('dp5', 'nuclear_quantum_effects_ignored', setting(pimd_crossover_rs_1000K=A['nqe_ignored']['rs'], **{NQ: A['nqe_ignored']['shift']})),
    ('dp5', 'nqe_shift_at_fixed_density', setting(**{NQ: A['nqe_fixed_density']})),
    ('dp6', 'nqe_900K_from_mislabelled_pimd_files', setting(**{NQ: A['nqe_900_from_mislabelled_files']})),
    ('dp7', 'small_cells_only_budget', setting(widom_rs_N2048=A['widom_N128'], classification=A['classification_small_cells'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('widom_discard_5ps', setting(widom_rs_N2048=R['widom_t0_5ps'])),
    ('widom_discard_20ps', setting(widom_rs_N2048=R['widom_t0_20ps'])),
    ('widom_three_point_fit', setting(widom_rs_N2048=R['widom_three_point_fit'])),
    ('widom_potential_energy_cv', setting(widom_rs_N2048=R['widom_potential_energy_cv'])),
    ('widom_grid_maximum', setting(widom_rs_N2048=R['widom_grid_maximum'])),
    ('pimd_crossover_from_tanh_fit', setting(pimd_crossover_rs_1000K=R['pimd_tanh_fit']['rs'], **{NQ: R['pimd_tanh_fit']['shift']})),
    ('classical_reference_N512', setting(**{NQ: R['classical_reference_N512']})),
    ('pimd_without_transient_removal', setting(pimd_crossover_rs_1000K=R['pimd_no_transient_removal']['rs'],
                                               **{NQ: R['pimd_no_transient_removal']['shift']})),
]
