"""Plausible wrong paths for gwbse-1d-chain-excitons; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived data by
benchmarks/construction/gwbse-1d-chain-excitons/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            if isinstance(value, dict) and isinstance(payload.get(field), dict):
                payload[field] = {**payload[field], **value}
            else:
                payload[field] = value
        return payload, trace
    return mutate


def bi2te3_gap(with_soc, without_soc):
    return {'Bi2Te3': {'with_soc': with_soc, 'without_soc': without_soc}}


ABLATIONS = [
    # dp1: trusting the run labels (BSE_G0W0 = G0W0 energies) for Te3 and As2S3
    ('dp1', 'run_labels_trusted', setting(screening_update_shift_eV=A['labels_trusted_dw'])),
    # dp1: same belief for Bi2Te3, E_b from the G0W0 gap and the BSE_G0W0 peak
    ('dp1', 'run_labels_trusted_bi2te3_binding', setting(bi2te3_binding_energy_eV=A['labels_trusted_eb'])),
    # dp2: over-generalizing: S3's runs assumed to share QP energies as well
    ('dp2', 's3_runs_assumed_same_energies', setting(screening_update_shift_eV={'S3': A['s3_same_energies_dw']})),
    # dp3: archived GW assumed to lack SOC, DFT SOC change added again
    ('dp3', 'soc_added_again', setting(scgw_direct_gap_eV=A['soc_added_again_gaps'])),
    # dp4: Bi2Te3 scGW gap = its G0W0 gap + mean self-consistency increment of the released chains
    ('dp4', 'self_consistency_increment_transferred',
     setting(scgw_direct_gap_eV=bi2te3_gap(A['increment_transfer_gap'],
                                          A['increment_transfer_gap'] - A['soc_added_again_gaps']['Bi2Te3']['with_soc']
                                          + A['soc_added_again_gaps']['Bi2Te3']['without_soc']),
             bi2te3_binding_energy_eV=A['increment_transfer_eb'])),
    # dp4: Bi2Te3 self-consistency assumed negligible (G0W0 gap reported)
    ('dp4', 'g0w0_gap_reported_for_bi2te3',
     setting(scgw_direct_gap_eV=bi2te3_gap(A['g0w0_gap_for_bi2te3'],
                                          A['g0w0_gap_for_bi2te3'] - A['soc_added_again_gaps']['Bi2Te3']['with_soc']
                                          + A['soc_added_again_gaps']['Bi2Te3']['without_soc']))),
    # dp5: first visible spectral feature taken as the optical gap, regardless of strength
    ('dp5', 'first_spectral_feature', setting(scgw_optical_gap_eV={'As2S3': A['first_feature_optical_As2S3']},
                                              screening_update_shift_eV=A['first_feature_dw'])),
    # dp6: transfer by a 1D screened-wire model with the DFPT supercell polarizability
    ('dp6', 'wire_model_transfer', setting(bi2te3_binding_energy_eV=A['wire_model_eb'])),
    ('dp6', 'wire_model_calibrated_on_S3_Te3', setting(bi2te3_binding_energy_eV=A['wire_model_calibrated_eb'])),
    ('dp6', 'mean_released_binding_energy', setting(bi2te3_binding_energy_eV=A['released_mean_eb'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
_dsoc_bi = A['soc_added_again_gaps']['Bi2Te3']['without_soc'] - A['soc_added_again_gaps']['Bi2Te3']['with_soc']
ROBUSTNESS = [
    ('bi2te3_gap_from_ip_onset', setting(scgw_direct_gap_eV=bi2te3_gap(R['bi2te3_gap_ip_onset'],
                                                                     R['bi2te3_gap_ip_onset'] + _dsoc_bi))),
    ('bi2te3_gap_from_first_ip_peak', setting(scgw_direct_gap_eV=bi2te3_gap(R['bi2te3_gap_ip_first_peak'],
                                                                          R['bi2te3_gap_ip_first_peak'] + _dsoc_bi))),
    ('s3_shift_from_ip_peak_shift', setting(screening_update_shift_eV={'S3': R['s3_dw_ip_shift']})),
    ('te3_soc_free_gap_k_resolved', setting(scgw_direct_gap_eV={'Te3': {'with_soc': A['soc_added_again_gaps']['Te3']['without_soc'],
                                                                       'without_soc': R['te3_without_soc_k_resolved']}})),
    ('soc_change_taken_at_the_qp_minimum_k', setting(scgw_direct_gap_eV={
        'Te3': {'with_soc': A['soc_added_again_gaps']['Te3']['without_soc'], 'without_soc': R['te3_without_soc_fixed_k']},
        'Bi2Te3': {'with_soc': A['soc_added_again_gaps']['Bi2Te3']['without_soc'], 'without_soc': R['bi2te3_without_soc_fixed_k']}})),
    ('bi2te3_binding_no_w_correction', setting(bi2te3_binding_energy_eV=R['eb_no_w_correction'])),
    ('bi2te3_binding_additive_mean_shift', setting(bi2te3_binding_energy_eV=R['eb_additive_mean_dw'])),
    ('bi2te3_binding_mean_ratio', setting(bi2te3_binding_energy_eV=R['eb_mean_ratio'])),
    ('bi2te3_binding_scaled_ratio', setting(bi2te3_binding_energy_eV=R['eb_ratio_scaled_by_self_consistency'])),
]
