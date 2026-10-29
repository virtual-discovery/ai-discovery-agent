"""Plausible wrong paths for tb2j-reference-selfconsistency-set; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/tb2j-reference-selfconsistency-set/make_task.py from the archived pair lists.
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
KFM, KCV, EN = 'fm_reference_kmin_FeP_rlu', 'converged_exchange_kmin_FeP_rlu', 'stabilization_energy_meV_per_ion'


def setting(kfm=None, kcv=None, energies=None, scale=None):
    """Replace k-vectors and/or individual stabilisation energies {material: {source: value}}; or scale all energies."""
    def mutate(payload, trace):
        if kfm is not None:
            payload[KFM] = list(kfm)
        if kcv is not None:
            payload[KCV] = list(kcv)
        en = copy.deepcopy(payload[EN])
        for mat, vals in (energies or {}).items():
            en[mat].update(vals)
        if scale is not None:
            en = {m: {s: (None if v is None else v*scale) for s, v in e.items()} for m, e in en.items()}
        payload[EN] = en
        return payload, trace
    return mutate


T = A['insulators_truncated']
RI = R['insulators_pairs_below_14A']
ABLATIONS = [
    ('dp1', 'FeP_exchange_truncated_at_14A', setting(kfm=A['kmin_FeP_from_pairs_below_14A'])),
    ('dp1', 'FeP_exchange_truncated_at_11A', setting(kfm=A['kmin_FeP_from_pairs_below_11A'])),
    ('dp2', 'paper_kmin_quoted', setting(kfm=A['kmin_FeP_paper_verbatim'])),
    ('dp2', 'commensurate_approximant_as_kmin', setting(kfm=A['kmin_FeP_commensurate_approximant'])),
    ('dp3', 'FeP_energies_from_kmin_converged_range_18A', setting(kfm=A['kmin_FeP_pairs_below_18A'],
                                                                  energies={'FeP': {'fm_reference': A['FeP_fm_energy_pairs_below_18A']}})),
    ('dp3', 'insulator_energies_from_short_lists', setting(energies={
        'FeCl2': {'fm_reference': T['FeCl2_fm_below_8A']},
        'FePS3_DFT': {'fm_reference': T['FePS3_DFT_fm_below_11A'], 'converged_reference': T['FePS3_DFT_conv_below_11A']}})),
    ('dp4', 'FeCl2_supercell_long_windows_skipped', setting(energies={'FeCl2': {'converged_reference': A['FeCl2_conv_long_windows_skipped']}})),
    ('dp4', 'FeP_supercell_long_windows_skipped', setting(energies={'FeP': {'converged_reference': A['FeP_conv_long_windows_skipped']}})),
    ('dp5', 'FeP_spiral_folded_as_listed', setting(kcv=A['kmin_conv_FeP_folded_as_listed'])),
    ('dp5', 'FeP_spiral_folded_as_listed_sum_over_five', setting(kcv=A['kmin_conv_FeP_folded_as_listed_sum5'])),
    ('dp6', 'FePS3_DFT+U_FM_duplicate_used', setting(energies={'FePS3_DFT+U': {'fm_reference': A['FePS3_DFT+U_fm_from_archived_FM_file']}})),
    ('dp7', 'FeP_fm_energy_of_its_own_LT_helix', setting(energies={'FeP': {'fm_reference': A['FeP_fm_energy_of_its_own_LT_helix']}})),
    ('dp7', 'FeP_reference_independence_assumed', setting(kcv=A['kmin_conv_FeP_assumed_reference_independent'])),
    ('dp8', 'bond_counted_once', setting(scale=0.5)),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('FeP_kmin_from_pairs_below_18A', setting(kfm=R['kmin_FeP_pairs_below_18A'])),
    ('FeP_kmin_from_pairs_below_23A', setting(kfm=R['kmin_FeP_pairs_below_23A'])),
    ('FeP_energies_on_relaxed_spiral', setting(energies={'FeP': R['FeP_relaxed_spiral']})),
    ('insulator_energies_from_shorter_lists', setting(energies={'FePS3_DFT': RI['FePS3_DFT'], 'NiO': RI['NiO'],
                                                               'FeCl2': {'fm_reference': RI['FeCl2_fm_below_11A']}})),
    ('FeP_spiral_folded_sum_over_five', setting(kcv=R['kmin_conv_FeP_sum5'])),
    ('FeP_spiral_folded_without_half_box_pairs', setting(kcv=R['kmin_conv_FeP_drop_half_box'])),
    ('FeP_spiral_bought_below_30A', setting(kcv=R['kmin_conv_FeP_bought_below_30A'])),
]
