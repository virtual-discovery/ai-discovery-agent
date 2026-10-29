"""Plausible wrong paths for amaranta-2d-anisotropic-ferromagnets; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/amaranta-2d-anisotropic-ferromagnets/make_task.py from the teacher model
(the authors' archived AMaRaNTA parameters) and the reference Monte Carlo.
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def change(ground_states=None, polar=None, anisotropy=None, bond=None, Tc=None, D=None, gap=None, predictions=None):
    def mutate(payload, trace):
        payload = copy.deepcopy(payload)
        if ground_states:
            payload['ground_states'].update(ground_states)
        if polar:
            payload['easy_axis_polar_deg'].update(polar)
        if anisotropy:
            payload['anisotropy_meV_per_Cr'].update(copy.deepcopy(anisotropy))
        if bond:
            payload['CrI3_bond_exchange_meV'] = dict(bond)
        if Tc is not None:
            payload['CrI3_Tc_K'] = Tc
        if D is not None:
            payload['CrSBr_dm_vector_meV'] = list(D)
        if gap is not None:
            payload['CrSBr_spin_wave_gap_meV'] = gap
        if predictions:
            payload['energy_predictions_meV_per_site'].update(predictions)
        return payload, trace
    return mutate


F = A['CrSiSe3_assumed_ferromagnet']
S7 = A['CrSBr_symmetric_offdiagonal']
ABLATIONS = [
    ('dp1', 'CrSiSe3_assumed_ferromagnetic', change(ground_states={'CrSiSe3': F['type']}, polar={'CrSiSe3': F['polar']},
                                                     anisotropy={'CrSiSe3': F['anisotropy']})),
    ('dp2', 'ordered_state_anisotropy_attributed_to_single_ion', change(anisotropy=A['mae_as_single_ion']['anisotropy'],
                                                                        predictions=A['mae_as_single_ion']['predictions'])),
    ('dp3', 'exchange_anisotropy_neglected', change(anisotropy=A['single_ion_only']['anisotropy'], polar=A['single_ion_only']['polar'],
                                                    predictions=A['single_ion_only']['predictions'])),
    ('dp4', 'CrI3_bond_tensor_copied_to_all_bonds', change(polar={'CrI3': A['CrI3_tensor_copied']['polar']},
                                                          anisotropy={'CrI3': A['CrI3_tensor_copied']['anisotropy']},
                                                          predictions=A['CrI3_tensor_copied']['predictions'])),
    ('dp5', 'CrI3_mirror_image_octahedral_frame', change(bond=A['CrI3_mirror_frame'])),
    ('dp6', 'CrI3_collinear_only_mapping', change(bond=A['CrI3_diagonal_only']['bond'], predictions=A['CrI3_diagonal_only']['predictions'])),
    ('dp7', 'CrSBr_offdiagonal_read_as_symmetric', change(D=S7['D'], polar={'CrSBr': S7['polar']}, anisotropy={'CrSBr': S7['anisotropy']},
                                                          gap=S7['gap'], predictions=S7['predictions'])),
    ('dp8', 'CrSBr_dm_copied_to_second_chain', change(predictions=A['CrSBr_uniform_dm']['predictions'])),
    ('dp8', 'CrSBr_mapped_in_2a_periodic_cell', change(D=A['CrSBr_aliased_2a']['D'], predictions=A['CrSBr_aliased_2a']['predictions'])),
    ('dp9', 'CrSBr_uniaxial_gap', change(gap=A['CrSBr_uniaxial_gap'])),
    ('dp10', 'Tc_mean_field', change(Tc=A['Tc_mean_field'])),
    ('dp10', 'Tc_nearest_neighbour_model', change(Tc=A['Tc_nn_only'])),
    ('dp11', 'single_site_four_state_at_noise_level', change(anisotropy=A['single_site_four_state']['anisotropy'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('CrI3_bond_from_actual_plane_normal', change(bond=R['CrI3_bond_actual_plane_normal'])),
    ('Tc_susceptibility_peak_largest_L', change(Tc=R['Tc_susceptibility_peak'])),
    ('Tc_binder_crossing_smallest_sizes', change(Tc=R['Tc_binder_small_pair'])),
    ('anisotropy_from_noisy_collective_states', change(anisotropy=R['anisotropy_from_noisy_collective_states'])),
    ('predictions_from_unsymmetrized_archived_tensor', change(predictions=R['predictions_unsymmetrized_tensor'])),
]
