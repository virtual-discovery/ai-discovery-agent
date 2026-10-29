"""Plausible wrong paths for raman-polymorph-orientation; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/raman-polymorph-orientation/variants.py from the reference solution's own purchases
(replayed with the lab's code), with the wrong analysis choice switched on.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        payload.update(changes)
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'raw_dft_frequencies', setting(A['raw_dft_frequencies'])),
    ('dp2', 'scale_without_axis_offset', setting(A['scale_without_axis_offset'])),
    ('dp3', 'nearest_line_peak_matching', setting(A['nearest_line_peak_matching'])),
    ('dp4', 'tetragonal_near_twin', setting(A['tetragonal_near_twin'])),
    ('dp5', 'lab_axes_polarizations_only', setting(A['lab_axes_polarizations_only'])),
    ('dp6', 'to_only_orientation_fit', setting(A['to_only_orientation_fit'])),
    ('dp6', 'to_only_face_b_phonon', setting(A['to_only_face_b_phonon'])),
    ('dp7', 'born_rows_unrotated', setting(A['born_rows_unrotated'])),
    ('dp7', 'born_rows_isotropic', setting(A['born_rows_isotropic'])),
    ('dp7', 'nac_without_eps_inf', setting(A['nac_without_eps_inf'])),
    ('dp8', 'minority_phase_ignored', setting(A['minority_phase_ignored'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('nonpolar_areas_centroid_calibration', setting(R['nonpolar_areas_centroid_calibration'])),
    ('low_window_spectra_only', setting(R['low_window_spectra_only'])),
    ('constant_linewidth_fit', setting(R['constant_linewidth_fit'])),
]
