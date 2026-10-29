"""Plausible wrong paths for zno-hydrogen-shallow-donor; each must fail grading.

Values come from solution/ablation_values.json, computed from the released energies by
benchmarks/construction/zno-hydrogen-shallow-donor/make_task.py with the same pipeline as the reference.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
NUMERIC = ('band_gap_eV', 'delta_mu_eV', 'transition_levels', 'V_O_plus_ever_stable', 'fermi_level_eV',
           'log10_concentration_cm3')


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


def pipeline(name, source=A, **extra):
    """Every numerical field as the named pipeline variant produces it."""
    values = {k: source[name][k] for k in NUMERIC}
    values.update(extra)
    return setting(**values)


def levels_with(defect, entries):
    def mutate(payload, trace):
        payload['transition_levels'] = dict(payload['transition_levels'], **{defect: entries})
        return payload, trace
    return mutate


def chain(*mutators):
    def mutate(payload, trace):
        for m in mutators:
            payload, trace = m(payload, trace)
        return payload, trace
    return mutate


ABLATIONS = [
    # dp1: semilocal energetics (PBE relaxations) instead of the hybrid single points; V_O^+ then looks stable
    ('dp1', 'pbe_pipeline', pipeline('pbe_pipeline', epr_centre={'defect': 'V_O', 'charge': 1})),
    ('dp1', 'pbe_pipeline_keep_oh_assignment', pipeline('pbe_pipeline')),
    # dp2: hydrogen chemical potential at the H2 limit, missing the Zn(OH)2 constraint
    ('dp2', 'hydrogen_reservoir_H2', pipeline('hydrogen_reservoir_H2')),
    # dp3: neutral shallow-donor cells taken at face value (electron filling conduction-band states at several k)
    ('dp3', 'no_band_filling', pipeline('no_band_filling')),
    ('dp3', 'no_band_filling_so_OH_O0_unstable_assign_V_O', pipeline('no_band_filling', epr_centre={'defect': 'V_O', 'charge': 1})),
    # dp4: band filling with equal k-point weights instead of the 1/8, 3/8 weights of the 2x2x2 mesh
    ('dp4', 'band_filling_equal_k_weights', pipeline('band_filling_equal_k_weights')),
    # dp5: V_O levels between every adjacent charge pair rather than between consecutive stable states
    ('dp5', 'V_O_adjacent_levels', chain(levels_with('V_O', V['V_O_adjacent_levels']), setting(V_O_plus_ever_stable=True))),
    # dp6: EPR assignment without the paramagnetism / population criterion
    ('dp6', 'epr_V_O_plus', setting(epr_centre={'defect': 'V_O', 'charge': 1})),
    ('dp6', 'epr_ionized_OH_O_plus', setting(epr_centre={'defect': 'OH_O', 'charge': 1})),
    ('dp6', 'epr_H_O_neutral', setting(epr_centre={'defect': 'H_O', 'charge': 0})),
    # dp7: hydrogen geometry read as antibonding (H on the far side of O from the Zn neighbour)
    ('dp7', 'site_antibonding_parallel', setting(oh_o_site='AB_parallel')),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('band_filling_unaligned_cbm', pipeline('band_filling_unaligned_cbm', source=R)),
    ('shallow_donors_at_cbm', pipeline('shallow_donors_at_cbm', source=R)),
    ('toten_energies', pipeline('toten_energies', source=R)),
    ('bulk_zno_gap', pipeline('bulk_zno_gap', source=R)),
    ('boltzmann_carriers', pipeline('boltzmann_carriers', source=R)),
    ('alignment_radius_4A', pipeline('alignment_radius_4A', source=R)),
]
