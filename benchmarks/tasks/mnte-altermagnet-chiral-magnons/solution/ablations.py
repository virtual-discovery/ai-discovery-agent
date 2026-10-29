"""Plausible wrong paths for mnte-altermagnet-chiral-magnons; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/mnte-altermagnet-chiral-magnons/calib.py and make_task.py from the archived data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FLIP = {'+': '-', '-': '+'}


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


def chirality(transform):
    def mutate(payload, trace):
        payload['in20_chirality'] = transform(payload['in20_chirality'])
        return payload, trace
    return mutate


def flip_all(ch):
    return {s: {k: FLIP[v] for k, v in ch[s].items()} for s in ch}


def copy_state_b(ch):           # HYSPEC pattern (same NMI sign as state B) used for both IN20 states
    return {'A': dict(ch['B']), 'B': dict(ch['B'])}


ABLATIONS = [
    ('dp1', 'shell_averaged_no_split', setting(split_exchange_difference_meV=A['no_split']['dJ'], in20_peak_qh=A['no_split']['peaks'])),
    ('dp2', 'anisotropy_gap_as_splitting', setting(max_splitting_q_rlu=A['q_gamma'])),
    ('dp3', 'in20_nominal_l', setting(in20_peak_qh=A['nominal_l']['peaks'])),
    ('dp3', 'in20_nominal_l_h_scaled', setting(in20_peak_qh=A['nominal_l_h_scaled']['peaks'])),
    ('dp4', 'ignore_in20_pnd', chirality(copy_state_b)),
    ('dp5', 'flip_for_positive_l_only', chirality(flip_all)),
    ('dp6', 'pnd_z_channel', chirality(flip_all)),
    ('dp7', 'max_on_measured_line', setting(max_splitting_q_rlu=A['q_measured_line'])),
    ('dp9', 'chiral_dipole_full_chirality', setting(split_exchange_difference_meV=A['chiral_dipole_full_chirality']['dJ'],
                                                    in20_peak_qh=A['chiral_dipole_full_chirality']['peaks'])),
    ('dp9', 'chiral_dipole_full_chirality_true_l', setting(split_exchange_difference_meV=A['chiral_dipole_full_chirality']['dJ'],
                                                           in20_peak_qh=A['chiral_dipole_full_chirality']['peaks_true_l'])),
    ('dp8', 'ordered_pair_convention', setting(split_exchange_difference_meV=A['half_dj']['dJ'], in20_peak_qh=A['half_dj']['peaks'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('in20_h_rescaled_to_instrument_a', setting(in20_peak_qh=R['h_scaled']['peaks'])),
    ('dJ_from_cut_g_only', setting(split_exchange_difference_meV=R['dj_g_only']['dJ'], in20_peak_qh=R['dj_g_only']['peaks'])),
    ('dJ_from_cut_h_only', setting(split_exchange_difference_meV=R['dj_h_only']['dJ'], in20_peak_qh=R['dj_h_only']['peaks'])),
    ('dJ_published_value', setting(split_exchange_difference_meV=R['dj_published']['dJ'], in20_peak_qh=R['dj_published']['peaks'])),
    ('equivalent_wave_vector', setting(max_splitting_q_rlu=R['q_equivalent'])),
]
