"""Plausible wrong paths for wurtzite-ferroelectric-switching; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived QE outputs by
benchmarks/construction/wurtzite-ferroelectric-switching/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def apply(**changes):
    """Set fields; dotted names set one key of a dict field (e.g. barrier_meV_per_fu.AlBN)."""
    def mutate(payload, trace):
        for name, value in changes.items():
            if '__' in name:
                field, key = name.split('__')
                payload[field][key] = value
            else:
                payload[name] = value
        return payload, trace
    return mutate


ABLATIONS = [
    # dp1: first-stage (collective) NEB taken as the (Al,B)N minimum-energy path; its saddle is the half-reversed midpoint
    ('dp1', 'AlBN_first_stage_NEB_as_MEP', apply(barrier_meV_per_fu__AlBN=A['barrier_AlBN_first_stage'],
                                                 AlBN_saddle_polarization_change_C_m2=A['collective_saddle_dP'],
                                                 AlBN_switching_mode='collective')),
    # dp2: homogeneous (collective) barrier of (Al,B)N assumed equal to AlN's
    ('dp2', 'collective_barrier_from_AlN', apply(AlBN_collective_barrier_meV_per_fu=A['collective_from_AlN'])),
    # dp3: normalization per atom (32 atoms) or per supercell instead of per formula unit
    ('dp3', 'per_atom_normalization', apply(barrier_meV_per_fu={'AlN': A['per_atom']['AlN'], 'AlBN': A['per_atom']['AlBN']},
                                            AlBN_collective_barrier_meV_per_fu=A['per_atom']['collective'],
                                            AlBN_intermediate={'energy_meV_per_fu': A['per_atom']['intermediate'],
                                                               'inverted_cations': 8})),
    ('dp3', 'per_supercell_energy', apply(barrier_meV_per_fu={'AlN': A['per_supercell']['AlN'], 'AlBN': A['per_supercell']['AlBN']},
                                          AlBN_collective_barrier_meV_per_fu=A['per_supercell']['collective'])),
    # dp4: raw Berry-phase values (no branch choice) used as polarizations
    ('dp4', 'raw_berry_phase_values', apply(spontaneous_polarization_C_m2={'AlN': A['raw_berry']['AlN'], 'AlBN': A['raw_berry']['AlBN']},
                                            AlBN_saddle_polarization_change_C_m2=A['raw_berry']['saddle'])),
    # dp5: nearest-branch unwrapping over too sparse a set of images (budget-driven stride)
    ('dp5', 'sparse_branch_tracking', apply(spontaneous_polarization_C_m2={'AlN': A['sparse_tracking']['AlN_stride2_half_path'],
                                                                           'AlBN': A['sparse_tracking']['AlBN_stride6']})),
    # dp6: Born-charge estimate over the whole path (end images only) to pick the branch
    ('dp6', 'born_charge_end_to_end_branch', apply(spontaneous_polarization_C_m2={'AlN': A['bec_full_path']['AlN'],
                                                                                  'AlBN': A['bec_full_path']['AlBN']})),
    # dp7: the metastable half-switched (Al,B)N intermediate taken as the nonpolar reference
    ('dp7', 'half_switched_state_as_reference', apply(spontaneous_polarization_C_m2__AlBN=A['half_switched_reference_AlBN'])),
    # dp8: Berry folder NN_image matched to NEB image NN (1-based) instead of NN+1
    ('dp8', 'berry_folder_off_by_one', apply(AlBN_saddle_polarization_change_C_m2=A['off_by_one_saddle']['plus'])),
    ('dp8', 'berry_folder_off_by_one_other_way', apply(AlBN_saddle_polarization_change_C_m2=A['off_by_one_saddle']['minus'])),
    # dp9: first sub-barrier / first intermediate of the multi-hump path taken as the saddle / intermediate
    ('dp9', 'first_subbarrier_as_saddle', apply(barrier_meV_per_fu__AlBN=A['first_subbarrier']['barrier'],
                                                AlBN_saddle_polarization_change_C_m2=A['first_subbarrier']['saddle_dP'])),
    ('dp9', 'first_minimum_as_intermediate', apply(AlBN_intermediate={'energy_meV_per_fu': A['first_subbarrier']['first_minimum_energy'],
                                                                      'inverted_cations': 4})),
    # dp10: mechanism inferred from the barrier alone (collective switching through the layered-hexagonal state)
    ('dp10', 'collective_mechanism', apply(AlBN_switching_mode='collective')),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('AlN_barrier_from_first_stage_run', apply(barrier_meV_per_fu__AlN=R['barrier_AlN_first_stage_run'])),
    ('saddle_located_by_parabola', apply(AlBN_saddle_polarization_change_C_m2=R['saddle_dP_parabolic_position'])),
    ('AlN_Ps_midpoint_reference', apply(spontaneous_polarization_C_m2__AlN=R['Ps_AlN_midpoint_reference'])),
    ('AlBN_barrier_vs_relaxed_endpoint', apply(barrier_meV_per_fu__AlBN=R['barrier_AlBN_vs_relaxed_endpoint'])),
]
