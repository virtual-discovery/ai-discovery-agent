"""Plausible wrong paths for halfheusler-strain-hinge; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/halfheusler-strain-hinge/make_task.py from the archived Hamiltonians with the lab's code.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


def hinge(**changes):
    def mutate(payload, trace):
        payload['hinge'] = {**payload['hinge'], **changes}
        return payload, trace
    return mutate


def nested(field, key, value):
    def mutate(payload, trace):
        payload[field] = {**payload[field], key: value}
        return payload, trace
    return mutate


def hinge_energy(values):
    return hinge(energy_k0_meV=values)


ABLATIONS = [
    # dp1: the file names say "strained"; the geometry is the equilibrium cubic one
    ('dp1', 'strain_from_file_names', setting(strain_percent=A['strain_from_file_names'])),
    ('dp1', 'strain_compressive_from_file_names', setting(strain_percent=A['strain_compressive_from_file_names'])),
    # dp2: gap read at Gamma (Gamma6 - Gamma8) instead of the full-zone band edges (CBM next to L)
    ('dp2', 'gap_at_gamma', setting(fundamental_gap_meV=A['gap_at_gamma'])),
    # dp3: LiBiZn treated as a zero-gap semiconductor with E_F at the Gamma8 touching point, or E_F from a coarse mesh
    ('dp3', 'fermi_at_gamma8', setting(libizn_fermi_level_meV=A['fermi_at_gamma8'])),
    ('dp3', 'fermi_coarse_mesh_8', setting(libizn_fermi_level_meV=A['fermi_coarse_mesh_8'])),
    ('dp3', 'fermi_mesh_12', setting(libizn_fermi_level_meV=A['fermi_mesh_12'])),
    # dp4: LiSbZn energies referred to the FPLO zero (arbitrary inside the gap, different in the two calculations)
    ('dp4', 'hinge_energy_model_zero', hinge_energy(A['hinge_energy_model_zero'])),
    ('dp4', 'removal_to_model_zero', hinge(removal_potential_eV=A['removal_to_model_zero'])),
    # dp5: hinge identity taken from the paper / the authors' mirrored real-space picture
    ('dp5', 'paper_termination', hinge(termination=A['paper_termination'])),
    ('dp5', 'paper_zn_s_weight', hinge(zn_s_weight_k0=A['paper_zn_s_weight'])),
    # dp6: unconverged prism (slow convergence of the k = 0 level)
    ('dp6', 'hinge_energy_small_prism_N10', hinge_energy(A['hinge_energy_small_prism_N10'])),
    ('dp6', 'hinge_energy_N12', hinge_energy(A['hinge_energy_N12'])),
    # dp7: hinge band window from a coarse k grid, clipped to the global gap, bottom assumed at k = 0, or one branch only
    ('dp7', 'window_grid5', hinge(band_window_meV=A['window_grid5'])),
    ('dp7', 'window_clipped_to_bulk_gap', hinge(band_window_meV=A['window_clipped_to_bulk_gap'])),
    ('dp7', 'window_bottom_at_k0', hinge(band_window_meV=A['window_bottom_at_k0'])),
    ('dp7', 'window_lower_branch_only', hinge(band_window_meV=A['window_lower_branch_only'])),
    # dp8: the strained-setup prism built with the orig integer basis (a different prism with (100) facets)
    ('dp8', 'str_hinge_identity_basis', lambda p, t: hinge_energy({**p['hinge']['energy_k0_meV'],
                                                                    'LiSbZn_strained_setup': A['str_hinge_identity_basis']})(p, t)),
    # dp9: the hinge potential applied at the X/X corner (paper's hinge): the level never leaves the gap
    ('dp9', 'removal_potential_wrong_corner', hinge(removal_potential_eV=A['removal_potential_wrong_corner'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('hinge_energy_N28_direct', hinge_energy(R['hinge_energy_N28_direct'])),
    ('hinge_energy_N20_direct', hinge_energy(R['hinge_energy_N20_direct'])),
    ('gap_from_grid_L_point', setting(fundamental_gap_meV=R['gap_grid_L_point'])),
    ('fermi_mesh_48', setting(libizn_fermi_level_meV=R['fermi_mesh_48'])),
    ('fermi_mesh_64_lab_max', setting(libizn_fermi_level_meV=R['fermi_mesh_64_lab_max'])),
    ('window_N12_refined', hinge(band_window_meV=R['window_N12_refined'])),
    ('window_N12_grid25', hinge(band_window_meV=R['window_N12_grid25'])),
    ('zn_weight_N16', hinge(zn_s_weight_k0=R['zn_weight_N16'])),
    ('removal_N16_linear', hinge(removal_potential_eV=R['removal_N16_linear'])),
    ('removal_N16_brentq', hinge(removal_potential_eV=R['removal_N16_brentq'])),
]
