"""Plausible wrong paths for gan-irradiation-defect-levels; each must fail grading.

Values come from solution/ablation_values.json, computed from the released calculations by
benchmarks/construction/gan-irradiation-defect-levels/make_task.py (reference.py holds the analysis).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELDS = ('band_gap_eV', 'transition_levels', 'divacancy_binding_energy_eV', 'frenkel_pair_energy_eV')


def whole(name):
    """Replace every graded field by the values of one complete wrong analysis."""
    def mutate(payload, trace):
        for f in FIELDS:
            payload[f] = A[name][f] if name in A else R[name][f]
        return payload, trace
    return mutate


def levels_of(name, defects):
    def mutate(payload, trace):
        for d in defects:
            payload['transition_levels'][d] = A[name][d]
        return payload, trace
    return mutate


def robust(name):
    def mutate(payload, trace):
        for f in FIELDS:
            payload[f] = R[name][f]
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'band_edges_from_R_point_supercell', whole('edges_R_point')),
    ('dp1', 'vbm_from_perfect_cell_fermi_energy', whole('edges_efermi')),
    ('dp2', 'host_energy_from_gamma_only_supercell', whole('host_gamma_only')),
    ('dp3', 'electronic_dielectric_tensor', whole('eps_inf')),
    ('dp3', 'no_finite_size_correction', whole('no_correction')),
    ('dp4', 'point_charge_term_without_alignment', whole('pc_only')),
    ('dp5', 'core_potentials_read_as_electrostatic_potentials', whole('alignment_wrong_sign')),
    ('dp6', 'single_interstitial_site_i1', whole('single_site_i1')),
    ('dp7', 'adjacent_charge_levels_without_stability', levels_of('adjacent_levels', ('V_N', 'Ga_i'))),
    ('dp8', 'band_edges_from_Ga_d_paw_run', whole('edges_Ga_d')),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('vbm_aligned_by_mean_core_potential', robust('edges_aligned_all_sites')),
    ('isotropic_static_dielectric_constant', robust('eps_0_isotropic')),
    ('alignment_sites_beyond_5.5_A', robust('region_beyond_5.5A')),
    ('alignment_sites_beyond_4.0_A', robust('region_beyond_4.0A')),
]
