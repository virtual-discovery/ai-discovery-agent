# Source and transformations

- **Paper:** S. K. Das, I. C. Fulga, R. Dhawan, H. C. Kandpal, J. van den Brink and J. I. Facio, *Electronic structure and hinge states of strained half-Heusler compounds LiSbZn and LiBiZn*, Phys. Rev. B 112, 125132 (2025), DOI 10.1103/PhysRevB.112.125132, arXiv:2504.20190.
- **Data:** Zenodo record 16927281, *Hinge electronic structure of strained half-Heuslers*, DOI 10.5281/zenodo.16927281, https://zenodo.org/records/16927281.
- **License:** CC BY 4.0 (record metadata). No file carries another license.
- **File used:** `strained_hH_zenodo.zip`, 2,681,691 bytes, md5 f66b568ce17e4055122facfe02e01e24, sha256 0d4ac9e4cb785c7f023050f4660655535f4b22e4963685475879b797f2a2c5bb.

## Members (md5)

| Member | md5 | Use |
|---|---|---|
| `wannier/+hamdata_lisbzn_eq_orig_setup.txt` | e46182ddfa7c6319d9e135a84f65a59c | lab model `LiSbZn_orig`; header staged |
| `wannier/+hamdata_lisbzn_eq_strained_setup.txt` | 1657fad07d1e43e9cc3d6952aa20dc08 | lab model `LiSbZn_strained_setup`; header staged |
| `wannier/+hamdata_libizn_eq_orig_setup.txt` | c387330316fa084fa2383ddec25b87d8 | lab model `LiBiZn_orig`; header staged |
| `wannier/+hamdata_libizn_eq_strained_setup.txt` | c341c104858f00961b8f13af8949fe95 | lab model `LiBiZn_strained_setup`; header staged |
| `code/strained_hH.py` | 1cf64e62695bc7f7b488fd75af947331 | staged unchanged |
| `kwant_ham/*.ham` (4 files) | see `fetch.py` | not used: verified identical to the hamdata hoppings |

## Construction

The scripts are in `benchmarks/construction/halfheusler-strain-hinge/`:

- `fetch.py` downloads the zip and verifies every member.
- `hamdata.py` parses each `+hamdata` file the way the authors' `process_hamdata()` does: hopping vector T = R + w_j − w_i, with an integer check on R. The matrix elements are stored unchanged (eV, complex128) in `environment/lab/data/<model>.npz`, with the lattice vectors and centres converted from Bohr to Å. The kwant pickles of the archive reproduce these blocks exactly.
- `tb.py` (copied into the lab) evaluates the Bloch Hamiltonian, slabs and prisms with numpy and scipy.
- `make_task.py` stages the header blocks (everything before the first hopping block), the code file and the record description, and writes the catalog. It then computes every reference value with the lab's own code path.

## What is hidden

Graded values are derived from the archived Hamiltonians and are never served:

- the strain of each model, computed from the lattice vectors;
- the band inversion;
- the fundamental gaps;
- the LiBiZn charge-neutral Fermi level;
- for the LiSbZn prism: the termination of the hinge that carries the in-gap state, its energy in the large-prism limit for both LiSbZn models, the energy window of the hinge band, its Zn-4s weight, and the hinge potential that moves it into the valence band.

The hopping blocks are not staged, so the agent cannot rebuild the models outside the budgeted service.

## Scope and deviations

- Retrospective analysis of the released Hamiltonians, not new DFT.
- The archive contains no strained calculation. Both `eq_strained_setup` files are the equilibrium cubic structures written in a body-centred tetragonal setting, and the kwant pickles named `*_strain.ham` are the same data. The strain dependence studied in the paper therefore cannot be reproduced, and the task grades what the data do contain.
- The prism geometry of the service is the one consistent with the archived hopping vectors. The authors' kwant code assigns each hopping to the opposite tag displacement, which leaves the spectrum unchanged but mirrors the real-space picture (see the construction README).
