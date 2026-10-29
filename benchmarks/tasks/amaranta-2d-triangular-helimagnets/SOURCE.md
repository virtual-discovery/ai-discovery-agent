# Source and transformations

- **Source:** A. Putatunda, G. Cuono, F. Orlando, S. Picozzi, *AMaRaNTA: automated first-principles exchange parameters in 2D magnets [Results]*, Materials Cloud Archive (2026), record wygeh-jrc57, DOI 10.24435/materialscloud:he-3k. Method paper: F. Orlando et al., arXiv:2509.20226.
- **License:** CC BY 4.0 (all files of the record).
- **Files used:**
  - `data_readable.zip` (53,905 bytes, md5 2020bf7c7813bb99f8af19c2030ef0e1): the `*.J_dict.txt` and `*.POSCAR.vasp` files of the eight compounds, and `other_data/PK_list.txt`.
  - Only the member `db.sqlite3` (the AiiDA provenance database) of three AiiDA exports, read with HTTP range requests, never the full archives:
    - `group2.aiida` (archive md5 2d3b5a776d16741b433ed22e710bd953; member sha256 880eb0f8f54f0aaca77fe02e2d2bda90df9d20f844eb8af9dc322eb07a588455): VCl2, VBr2, VI2, NiBr2;
    - `group3.aiida` (md5 06656233b366976650f3fc420c4ca324; sha256 c09d1cbd039409fd3c77cc05f3a5eb7ab89cb5097f0ab7eefd1a2b16f27c8eb2): MnCl2, MnBr2, MnI2;
    - `group10.aiida` (md5 f33b7f9621f32aa4e199d15078794d9b; sha256 80e067e68f2127bd08529222122781a3b4f3f660ccac238793fbae39a4d1e9ef): NiCl2.

## Construction

The scripts are in `benchmarks/construction/amaranta-2d-triangular-helimagnets/`:

- `fetch.py` downloads the zip and the three `db.sqlite3` members, checking the zip md5, each member's CRC32 from the archive's central directory and its sha256.
- `extract.py` reads, for the eight AmarantaWorkChains listed in `PK_list.txt`, the input structure, the authors' `J_dict`, the three pairs ("couples") AMaRaNTA used, and all 48 four-state VASP calculations (supercell, constrained moment directions, total energies, convergence flags). It asserts that the AiiDA `J_dict` equals the `*.J_dict.txt` files and that each archived parameter equals (E1 − E2 − E3 + E4)/4 of its four stored energies. The NiCl2 structure comes from the AiiDA input, because the zip has no NiCl2 POSCAR.
- `model.py` builds the teacher model, which is the AMaRaNTA Hamiltonian with the archived parameters:
  - the nearest-neighbour tensor sits on the archived 1NN pair and is carried to the other 1NN bonds by the site point group of the magnetic ion (D3d, or C2h for the distorted cells), with J(−d) = J(d)ᵀ;
  - scalar J2 and J3 sit on the bonds of exactly the archived 2NN and 3NN lengths (AMaRaNTA groups distances within 1 %, so in VBr2, VI2 and MnCl2 its "2NN" is the second, longer nearest-neighbour bond);
  - the single-ion coefficient of (eᶻ)² is K = 2A, because AMaRaNTA's formula (E₊z − E₊x − E₋x + E₋z)/4 returns half of it.

  Against the archived spin-orbit four-state energies, the model reproduces all within-set energy differences to 0.003–0.03 meV rms (VI2: 0.26 meV); `calibration.json` has the full numbers.
- `make_task.py` writes the catalog, the lab's model parameters (`environment/lab/private/model.json`), the hidden targets and the ablation values.

The lab adds deterministic Gaussian noise (σ = 0.02 meV per supercell energy, seeded by a hash of the request) and an arbitrary constant per magnetic site.

The exported `other_data/exchanges.json` repeats parameters across materials (for example, the MnBr2 entry carries the MnI2 values). It was not used; the AiiDA `J_dict` nodes and the text files agree with each other.

## What is hidden

- The model parameters: all exchange tensors and the single-ion term.
- The graded targets, computed from the teacher model without noise:
  - each compound's classical ground-state type and Luttinger–Tisza wave vector (anisotropy shifts it by at most 0.006 Å⁻¹);
  - the energies of the 24 `to_predict` configurations, in 7×7 or 4×7 supercells that the lab does not serve;
  - the easy axis of VI2's collinear stripe, which is tilted 30.7° from c in the ac plane.

## Scope

The lab is a surrogate built from the authors' archived parameters, not new DFT. The anisotropies (0.01–0.1 meV) are at or below the precision the authors quote (0.1 meV). The task tests inference within this model.
