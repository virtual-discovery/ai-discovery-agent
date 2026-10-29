# Source and transformations

- **Source:** I. C. Skogvoll, E. Fransson, B. Williamson, N. Bristowe, S. M. Selbach, P. Erhart, *Machine-learning-guided exploration of domain walls in the hybrid improper ferroelectric Ca3Ti2O7*, arXiv:2609.35346 (2026). Data: Zenodo record 22801877, DOI 10.5281/zenodo.22801877.
- **License:** CC BY 4.0 (all files of the record).
- **File used:** `Data.zip` (5,983,790 bytes, md5 4007d33af8f3b2b26055b4fbc711abeb), from which:
  - `Potentials/Potential-qNEP-w-BECs.txt`, `Potential-qNEP.txt`, `Potential-NEP.txt` (served by the lab, copied unchanged to `environment/lab/potentials/`);
  - `Bulk-structures/Cmc2_1-qNEP-w-BECs.xyz` and `I4mmm-qNEP-w-BECs.xyz` (agent catalog, coordinates rounded to 1e-6 Å);
  - `Domain-walls/qNEP-w-BECs/*.xyz` and `Training-sets/*.xyz` (construction checks only, not staged).

## Construction

Scripts are in `benchmarks/construction/ca3ti2o7-hybrid-improper-domain-walls/`.

- `environment/lab/nep_numpy.py` is a numpy implementation of NEP4 and NEP4-charge (qNEP: machine-learned charges, Ewald electrostatics, Born effective charges), written for this task from the published NEP formulas. It reproduces the NEP_CPU reference forces and virials of the NEP_CPU qNEP test case to 1e-13 eV/Å, and on the released training set it gives energy/force/BEC errors of 1.8 meV/atom, 0.07 eV/Å and 0.19 e for qNEP-w-BECs. The authors' 11 qNEP-w-BECs domain-wall structures have residual forces ≤ 0.02 eV/Å with it, so it is the authors' forward model.
- `builder.py`/`run_wall.py` build two-wall supercells from the bulk structure with parent-symmetry domain variants and relax them with the lab code (`relax.py`, L-BFGS). Rebuilt walls reproduce the authors' archived structures' energies (100-180-ψ_z 158.3 vs 157.9 mJ/m²; 100-180-φ 185.5 vs 185.6).
- `make_task.py` turns the construction relaxations (`runs/*.json`) into the hidden targets, ablation values and the agent catalog.

## What is hidden

The domain-wall structures of the archive and all wall energies. Graded: the converged neutral ψ and φ wall energies, the lowest layer ψ wall energy, the slope of the charged ψ-wall energy with domain width (all qNEP-w-BECs, computed by the builder with the lab code path), and the ψ/φ wall profile types.

## Scope

The forward model is the authors' potential, not DFT. The targets are builder calculations with that potential (released-start task with a released evaluator), not numbers quoted in the paper.
