"""Plausible wrong paths for fresh-inp-polar-phonons; each must fail grading.

Every value comes from solution/ablation_values.json, computed by
benchmarks/construction/fresh-inp-polar-phonons/make_task.py from real calculations made through the
task's own lab code (task_lab.run), not from arbitrary perturbations.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELDS = ['a0_A', 'omega_TO_cm-1', 'omega_LO_cm-1', 'born_charge_In', 'eps_inf', 'eps_0',
          'dlnZ_dlnV', 'dln_eps_inf_dlnV', 'gruneisen_TO', 'LOTO_splitting_under_compression',
          'strained_omega_LO_q_z_cm-1', 'strained_omega_LO_q_x_cm-1']


def use(values, fields=FIELDS):
    """Replace the listed submission fields by the values a wrong (or alternative) path produced."""
    def mutate(payload, trace):
        for field in fields:
            payload[field] = values[field]
        return payload, trace
    return mutate


ABLATIONS = [
    # dp1: Gamma-centred mesh (validated on energies/pressure/TO, which converge there), DFPT on 8x8x8 Gamma-centred.
    ('dp1', 'gamma_centred_k8', use(A['gamma_centred_k8'])),
    # dp2: cheap 4x4x4 shifted mesh for eps_inf and Z* because the LO frequency looks converged (errors cancel).
    ('dp2', 'k4_shifted_tensors', use(A['k4_shifted_tensors'])),
    # dp3: converged tensors at a0, but the second volume on the cheap 4x4x4 mesh (the k error does not cancel).
    ('dp3', 'side_volume_on_cheaper_mesh', use(A['side_volume_on_cheaper_mesh'])),
    ('dp3', 'side_volume_on_cheaper_mesh_expansion', use(A['side_volume_on_cheaper_mesh_expansion'])),
    # dp4: experimental lattice parameter instead of the protocol's own a0.
    ('dp4', 'experimental_lattice', use(A['experimental_lattice'])),
    # dp5: the analytic ph.x frequency reported as the LO, or the nonanalytic term screened by eps_0.
    ('dp5', 'analytic_frequency_as_LO', use(A['analytic_frequency_as_LO'])),
    ('dp5', 'nonanalytic_term_with_eps0', use(A['nonanalytic_term_with_eps0'])),
    # dp6: strained cell treated with isotropic tensors, with the TO modes paired to the wrong q, or as cubic.
    ('dp6', 'strained_scalar_tensors', use(A['strained_scalar_tensors'])),
    ('dp6', 'strained_swapped_TO_pairing', use(A['strained_swapped_TO_pairing'])),
    ('dp6', 'strained_as_cubic', use(A['strained_as_cubic'])),
    # dp7: trend inferred from the Born charge alone (Z* falls under compression -> 'splitting shrinks'),
    #      ignoring the faster fall of eps_inf.
    ('dp7', 'splitting_trend_from_born_charge_only', use(A['splitting_trend_from_born_charge_only'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, use(values)) for name, values in R.items()]
