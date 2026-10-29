"""Plausible wrong paths for hydride-coulomb-mustar; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/hydride-coulomb-mustar/make_task.py with the lab's solver on the released data.
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


ABLATIONS = [
    ('dp1', 'mu_star_equals_bare_mu', setting(mu_star_cdos_1eV_Mg2IrH6_20GPa=A['mu_star_bare_mu'])),
    ('dp1', 'mu_star_morel_anderson_rescaling', setting(mu_star_cdos_1eV_Mg2IrH6_20GPa=A['mu_star_morel_anderson'])),
    ('dp2', 'mu_from_smeared_NF', setting(tc_full_W_K=A['tc_smeared_NF_mu_sigma0.3'])),
    ('dp3', 'mu_from_smallest_sigma', setting(tc_full_W_K=A['tc_mu_rescaled_sigma0.01'])),
    ('dp4', 'constant_dos_solver', setting(tc_full_W_K=A['tc_constant_dos_solver'])),
    ('dp5', 'allen_dynes_with_archived_mucAD', setting(tc_full_W_K=A['tc_allen_dynes_archived_mucAD'])),
    ('dp6', 'tc_interpolated_in_pressure', setting(tc_full_W_K=A['tc_interpolated_in_pressure'])),
    ('dp7', 'conventional_mustar_0.16', setting(tc_full_W_K=A['tc_conventional_mustar_0.16'])),
    ('dp7', 'mu_held_at_0GPa_value', setting(tc_full_W_K=A['tc_constant_mu_from_0GPa'])),
    ('dp8', 'NF_per_conventional_cell', setting(N_EF_H3S_180GPa_per_fu_per_spin=A['nef_per_cell'])),
    ('dp9', 'mu_star_equivalent_in_vdos_solver', setting(mu_star_cdos_1eV_Mg2IrH6_20GPa=A['mu_star_vdos_7eV'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('bare_mu_in_vdos_without_kappa', setting(tc_full_W_K=R['tc_bare_mu_route'])),
    ('nearest_pressure_W_EF', setting(tc_full_W_K=R['tc_nearest_W00'])),
    ('mu_star_matched_to_lab_full_W_Tc', setting(mu_star_cdos_1eV_Mg2IrH6_20GPa=R['mu_star_vs_lab_TcW'])),
    ('NF_nearest_grid_point', setting(N_EF_H3S_180GPa_per_fu_per_spin=R['nef_nearest_grid'])),
]
