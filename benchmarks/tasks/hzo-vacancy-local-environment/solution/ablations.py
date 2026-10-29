"""Plausible wrong paths for hzo-vacancy-local-environment; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/hzo-vacancy-local-environment/make_task.py from the released data, using the
oracle's own purchases (the same site lists solve.py buys) and the free x = 0.50 table.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELDS = ('predicted_formation_energies_eV', 'vacancy_Zr_neighbour_distribution', 'vacancy_site_fraction',
          'lda_common_slope_model')


def use(entry, fields=FIELDS):
    def mutate(payload, trace):
        for field in fields:
            payload[field] = entry[field]
        return payload, trace
    return mutate


ABLATIONS = [
    # dp1: what controls E_f -- Zr count alone (sublattices pooled), or geometry descriptors instead of cation identity
    ('dp1', 'pooled_nZr_regression', use(A['dp1_pooled_nZr_regression'])),
    ('dp1', 'geometry_descriptor_model', use(A['dp1_geometry_descriptor_model'])),
    # dp2: composition-blind local model fitted to all compositions (shift absorbed into the Zr slopes)
    ('dp2', 'composition_blind_pooled_model', use(A['dp2_composition_blind_pooled_model'])),
    # dp3: x = 0.50 local model transferred to x = 0.25 / 0.75 without matched-environment purchases
    ('dp3', 'x050_model_transferred', use(A['dp3_x050_model_transferred'])),
    # dp4: environment statistics of the 96-atom cells instead of the macroscopic random solid solution
    ('dp4', 'supercell_decoration_statistics', use(A['dp4_supercell_decoration_statistics'])),
    ('dp4', 'fixed_composition_cell_statistics', use(A['dp4_fixed_composition_cell_statistics'])),
    # dp5: Boltzmann weights without configurational degeneracy, or the lowest-energy class only
    ('dp5', 'no_configurational_degeneracy', use(A['dp5_no_configurational_degeneracy'])),
    ('dp5', 'lowest_energy_class_only', use(A['dp5_lowest_energy_class_only'])),
    # dp6: rare environment classes never purchased at a composition dropped from the population
    ('dp6', 'unpurchased_classes_dropped', use(A['dp6_unpurchased_classes_dropped'])),
    # dp7: CASTEP sites classified with the VASP x = 0.50 cell's environments instead of the authors' labels
    ('dp7', 'castep_sites_given_vasp_cell_environments', use(A['dp7_castep_sites_given_vasp_cell_environments'])),
    # dp8: PBE description reported for LDA (functional transfer)
    ('dp8', 'pbe_description_reported_for_lda', use(A['dp8_pbe_description_reported_for_lda'])),
    ('dp8', 'castep_pbe_reported_for_lda', use(A['dp8_castep_pbe_reported_for_lda'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, use(entry)) for name, entry in R.items()]
