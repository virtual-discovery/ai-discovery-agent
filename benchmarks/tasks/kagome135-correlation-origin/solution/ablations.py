"""Plausible wrong paths for kagome135-correlation-origin; each must fail grading.

Values come from solution/ablation_values.json, written by
benchmarks/construction/kagome135-correlation-origin/make_task.py from lab-path SSMF evaluations of the
archived models (the same evaluator the service runs).
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
    # dp1: nominal valence counts (27/30/33), ignoring the Cs 6s electron donated to the d-p manifold
    ('dp1', 'nominal_electron_counts', setting(electron_counts=A['nominal_counts'])),
    ('dp1', 'nominal_count_CsCr3Sb5_Z', setting(CsCr3Sb5_U6=A['CrSb_nominal'])),
    ('dp1', 'nominal_count_CsCr3Bi5', setting(CsCr3Bi5_U4p5=A['CrBi_nominal'])),
    # dp2: CrBi_9.2 taken as the CsCr3Sb5-lattice model (folder name read as the c axis of CsCr3Sb5)
    ('dp2', 'swapped_CsCr3Bi5_models', setting(CsCr3Bi5_U4p5=A['CrBi_swapped'])),
    ('dp2', 'swapped_CsCr3Bi5_sigma', setting(sigma_eV=A['sigma_swapped_CrBi'])),
    # dp3: Wannier90's default d order (dz2, dxz, dyz, dx2-y2, dxy) instead of the archive's projection order
    ('dp3', 'w90_default_orbital_order', setting(most_renormalized_orbital=A['orbital_w90_default_order'])),
    # dp4: the paper's/recalled label (B3g = dyz) instead of the archive's local frame (x radial)
    ('dp4', 'recalled_paper_frame_label', setting(most_renormalized_orbital=A['orbital_paper_frame'])),
    # dp5: decomposition by a single path or by cross-compound regression
    ('dp5', 'one_sided_path_Ti_hamiltonian', setting(filling_share=A['share_one_sided_Ti'])),
    ('dp5', 'one_sided_path_Cr_hamiltonian', setting(filling_share=A['share_one_sided_Cr'])),
    ('dp5', 'cross_compound_regression', setting(filling_share=A['share_regression'])),
    # dp6: kinetic scale from a windowed PDOS, from the d-d block only, or orbital-averaged
    ('dp6', 'sigma_windowed_pdos', setting(sigma_eV=A['sigma_windowed'])),
    ('dp6', 'sigma_d_block_only', setting(sigma_eV=A['sigma_d_only'])),
    ('dp6', 'sigma_orbital_average', setting(sigma_eV=A['sigma_mean_orbital'])),
    # dp7: Hund versus Mott read at the wrong (nominal) filling: cascade into the J dependence
    ('dp7', 'hund_ratio_nominal_count', setting(hund_ratio=A['hund_nominal'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('denser_mesh_12x12x8', setting(CsCr3Sb5_U6=R['CrSb_mesh12'], hund_ratio=R['hund_mesh12'],
                                    CsCr3Bi5_U4p5=R['CrBi_mesh12'], filling_share=R['share_mesh12'])),
    ('sigma_from_pdos_histogram', setting(sigma_eV=R['sigma_mesh_histogram'])),
]
