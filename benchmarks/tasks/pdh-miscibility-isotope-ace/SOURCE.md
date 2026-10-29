# Source and transformations

- **Source:** M. Qamar, A. Ngoipala, M. Mrovec, M. Vandichel, R. Drautz et al., *Hydride formation and phase separation in palladium nanoparticles from a transferable atomic cluster expansion potential*, arXiv:2606.09341 (2026).
- **Data:** Zenodo record 20591788 (DOI 10.5281/zenodo.20591788), https://zenodo.org/records/20591788.
- **License:** CC BY 4.0 (all files).

## Files used

| File | Size | Checksum |
|---|---|---|
| `Pd-H_ACE_potential_files_and_training_data.zip` | 101,757,605 bytes | md5 0308b4e0ed87805975d37eaeaacf976a |
| `potential_files/ace/PdH_PBE_v03a.yaml` (member) | | sha256 1b4269217bf6948c7cecef926dfb4d27cb9a5e4e9ba64d1d96358e87171d0eec |
| `potential_files/ace/df_train.pckl.gzip` (member) | | sha256 a2abf2528b038f1e65bfc4400aa3e58b762af424c076752bdc223a07efc4ecd1 |
| `potential_files/ace/df_test.pckl.gzip` (member) | | sha256 2031d27e273e5181cf60e4f9e498245d8aaeea2ea71461383351b5cb991e6970 |

## Transformations

The scripts are in `benchmarks/construction/pdh-miscibility-isotope-ace/`.

- **Potential conversion.** `convert_potential.py` converts the B-basis ACE potential into the equivalent C-tilde basis, using the library's own exact transformation (python-ace 0.4.0rc1). It writes `PdH_ace_ctilde.json`, in which no parameter is changed.
- **Evaluator.** `ace_numpy.py` is a pure-numpy port of the ACE C-tilde evaluator: simplified-Bessel radial basis, complex harmonics, Finnis–Sinclair shifted-scaled embedding and ZBL inner core.
- **Validation.** `validate_ace.py` checks it against python-ace on 60 released test structures and close-contact cases. The maximum differences are:
  - energy: 6e-12 eV/atom;
  - force: 4e-8 eV/Å;
  - stress: 1e-9 eV/Å³.
- **DFT subset.** `extract_dft.py` copies a validation subset of the authors' DFT records (VASP PBE) into the lab: PdHₓ 16-Pd cells, Pd fcc and PdH references, H₂ molecules, and Pd₃₁ vacancy cells. The energies are the authors' `energy_corrected`, which shares the ACE energy zero.

## What is hidden

None of the graded numbers is archived by the authors. All are reference calculations with the released potential, run through the lab's own code (`reference.py`) and analysed in `make_task.py`:

- the dilute O-site solution energy at 0 K with ZPE;
- the static and ZPE-corrected T–O differences;
- the dilute lattice expansion;
- the T-site fraction at 300 K;
- the Sieverts ratio K_H/K_D at 300 K;
- the incoherent plateau pressure at 300 K;
- the incoherent and coherent mean-field critical temperatures.

They are stored only under `environment/lab/private/`.

The gas-phase data in the catalog were computed from standard experimental spectroscopic constants of H₂ and D₂: anharmonic zero-point energy, rotor levels with nuclear-spin statistics, and ideal-gas translation.
