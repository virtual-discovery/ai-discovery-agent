"""Plausible wrong paths for mnte-altermagnet-ahc-doping; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/mnte-altermagnet-ahc-doping/make_task.py from the archived data and the lab's evaluator.
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


def fields(values):
    return setting(**{k: v for k, v in values.items() if v is not None})


ABLATIONS = [
    ('dp1', 'mesh_100', fields(A['mesh_100'])),
    ('dp1', 'mesh_150', fields(A['mesh_150'])),
    ('dp1', 'mesh_200', fields(A['mesh_200'])),
    ('dp2', 'vbm_from_band_path', setting(hole_density_0p30_below_vbm_cm3=A['vbm_from_band_path'])),
    ('dp2', 'vbm_from_dft_grid', setting(hole_density_0p30_below_vbm_cm3=A['vbm_from_dft_grid'])),
    ('dp3', 'vbm_from_outcar_efermi', setting(hole_density_0p30_below_vbm_cm3=A['vbm_from_outcar_efermi'])),
    ('dp4', 'spin_degeneracy_twice', fields(A['spin_degeneracy_twice'])),
    ('dp5', 'raw_forbidden_components', setting(forbidden_components_max_abs_S_per_cm=A['raw_forbidden_components'])),
    ('dp6', 'symmetries_dat_literal', setting(allowed_hall_components=A['symmetries_dat_literal'])),
    ('dp6', 'inplane_isotropic', setting(allowed_hall_components=A['inplane_isotropic'])),
    ('dp7', 'electron_doped_from_model', setting(electron_doped_sigma_xy_S_per_cm=A['electron_doped_from_model'])),
    ('dp8', 'room_temperature_smearing', fields(A['room_temperature_smearing'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('mesh_300', fields(R['mesh_300'])),
    ('model_values_400', fields(R['mesh_400_model_values'])),
    ('wt_in_efermi_as_vbm', setting(hole_density_0p30_below_vbm_cm3=R['wt_in_efermi_as_vbm'])),
    ('step_occupation_400', fields(R['step_occupation'])),
]
