"""Plausible wrong paths for tb2j-groundstate-cuo-multiferroic; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/tb2j-groundstate-cuo-multiferroic/make_task.py from the archived pair lists.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELD = {'kmin': 'kmin_FM_reference_rlu', 'stability': 'stable_at_own_reference', 'converged_q': 'converged_propagation_vector_rlu',
         'chain': 'chain_exchange_meV', 'mean': 'mean_magnon_energy_meV', 'max': 'magnon_band_maximum_meV'}


def setting(values, keys=None):
    """Replace fields from an ablation/robustness record; k-min components are merged into the kmin dict."""
    def mutate(payload, trace):
        for key, value in values.items():
            if keys is not None and key not in keys:
                continue
            if key in ('ligand_corrected', 'uncorrected'):
                payload[FIELD['kmin']] = {**payload[FIELD['kmin']], key: value}
            elif key in FIELD:
                payload[FIELD[key]] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'fm_reference_exchange_as_final_model', setting(A['fm_reference_as_final'])),
    ('dp2', 'uncorrected_exchange_judged_converged', setting(A['uncorrected_exchange_as_final'])),
    ('dp3', 'exchange_truncated_at_6A', setting(A['exchange_truncated_at_6A'])),
    ('dp4', 'commensurate_approximant_reported_as_kmin', setting(A['commensurate_approximant_as_kmin'])),
    ('dp4', 'paper_SI_kmin_quoted_verbatim', setting(A['paper_SI_kmin_verbatim'])),
    ('dp4', 'experimental_AF2_vector_quoted', setting(A['memorised_experimental_AF2_as_kmin'])),
    ('dp5', 'symmetrized_afm_reference_exchange', setting(A['symmetrized_afm_exchange'], keys=('converged_q', 'mean', 'max'))),
    ('dp5', 'symmetrized_afm_reference_exchange_magnons_only', setting(A['symmetrized_afm_exchange'], keys=('mean', 'max'))),
    ('dp6', 'hamiltonian_eigenvalues_as_magnon_energies', setting(A['hamiltonian_eigenvalues_as_magnons'])),
    ('dp7', 'spin_length_one', setting(A['spin_length_one'])),
    ('dp7', 'bond_counted_twice_in_chain_exchange', setting(A['bond_counted_twice_convention'])),
    ('dp7', 'unit_vector_bond_convention_for_chain_exchange', setting(A['unit_vector_bond_convention'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('kmin_from_exchange_within_8A', setting(R['kmin_r8'])),
    ('kmin_from_exchange_within_10A', setting(R['kmin_r10'])),
    ('magnons_from_exchange_within_6.5A', setting(R['magnons_native_r6.5'])),
    ('magnons_from_exchange_within_10A', setting(R['magnons_native_r10'])),
    ('magnons_4x4x4_grid_unrefined_maximum', setting(R['magnons_grid4_unrefined_max'])),
    ('chain_exchange_from_weakest_single_bond', setting({'chain': R['chain_single_bond_extremes']['chain_low']})),
    ('chain_exchange_from_strongest_single_bond', setting({'chain': R['chain_single_bond_extremes']['chain_high']})),
]
