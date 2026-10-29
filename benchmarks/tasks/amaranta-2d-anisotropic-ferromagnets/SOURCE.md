# Source and transformations

- **Source:** A. Putatunda, G. Cuono, F. Orlando, S. Picozzi, *AMaRaNTA: automated first-principles exchange parameters in 2D magnets [Results]*, Materials Cloud Archive (2026), record wygeh-jrc57, DOI 10.24435/materialscloud:he-3k. Method paper: F. Orlando et al., arXiv:2509.20226.
- **License:** CC BY 4.0 (all files of the record).
- **Files used:**
  - `data_readable.zip` (53,905 bytes, md5 2020bf7c7813bb99f8af19c2030ef0e1): the `*.J_dict.txt` and `*.POSCAR.vasp` files of the Cr compounds and `other_data/PK_list.txt`.
  - Only the member `db.sqlite3` (the AiiDA provenance database) of five AiiDA exports, read with HTTP range requests, never the full archives:
    - `group1.aiida` (archive md5 12d6c201e55cf12df2279b7d4874870d; member sha256 1bf3450a671b696c889a79dcef80c1260f155b9e4fa42724659e30f70867d1a5): CrCl3, CrBr3, CrGeTe3, CrSiSe3 (and CrPS4, not used in the task);
    - `group4.aiida` (md5 63429af96ec9978edb12a769fa52c5d5; sha256 aa4540abf9304336b0343220e2f8bf5e0f3fc81a1ff054b7625487bd1d23b94f): CrI3 (and CrBrO, not used);
    - `group7.aiida` (md5 6d7a6de21ba04f0263af09fea094c685; sha256 a88112265cc4228f3450a303101c84e116d25952149e0b93170866713a4aabf5): CrSBr;
    - `group5.aiida` (md5 d977634a0b59b0b5b08ed47cadeb5829; sha256 07898d7f148a70d252600e4c8d7b9e4160b5a69a18926808c5387e5545d8ea94) and `group10.aiida` (md5 f33b7f9621f32aa4e199d15078794d9b; sha256 80e067e68f2127bd08529222122781a3b4f3f660ccac238793fbae39a4d1e9ef): CrClO and CrSiTe3, read for the model validation only.

## Construction

The scripts are in `benchmarks/construction/amaranta-2d-anisotropic-ferromagnets/`:

- `fetch.py` downloads the zip and the five `db.sqlite3` members. It checks the zip md5, and each member's CRC32 from the archive's central directory and its sha256.
- `extract.py` reads, for the ten Cr AmarantaWorkChains listed in `PK_list.txt`, the following and writes them to `pool.json`:
  - the input structure;
  - the authors' `J_dict`;
  - the three pairs ("couples") AMaRaNTA used;
  - all 120 four-state sets, 480 VASP calculations (supercell, constrained moment directions, total energies, convergence flags).

  It asserts that the AiiDA `J_dict` equals the `*.J_dict.txt` files, that each archived parameter equals (E1 − E2 − E3 + E4)/4 of its four stored energies, and that every run converged.
- `model.py` builds the teacher model, the AMaRaNTA Hamiltonian with the archived parameters:
  - The nearest-neighbour tensor sits on the archived 1NN pair. It is first averaged over the operations that leave that bond invariant, which removes symmetry-forbidden parts of 0.0003–0.04 meV. It is then carried to every other bond of the class by the layer's space-group operation that maps the pair onto it, found from the full structure with ligands. Spins transform as axial vectors, and a reversed bond carries the transpose. The DM vectors of the two Cr chains of CrSBr are therefore opposite. CrCl3's 0.1 % monoclinic distortion breaks C3, so its tensor is carried by the pseudo-D3d operations.
  - The scalar J2 and J3 sit on the bonds whose lengths are within 1 % of the archived 2NN and 3NN pairs (AMaRaNTA's grouping). In the CrSiSe3 structure the layer straddles the cell boundary in z, so the archived pairs carry a −c offset; positions are unwrapped.
  - The single-ion coefficient of (e^z)² is K = 2A, because AMaRaNTA's formula (E₊z − E₊x − E₋x + E₋z)/4 returns half of it.

  The model reproduces the archived spin–orbit four-state energies within each set to 0.002–0.13 meV rms. Taking K = A fits worse. `calibration.json` has the full numbers.
- `make_task.py` writes the following:
  - the catalog;
  - the lab's model parameters (`environment/lab/private/model.json`);
  - the hidden targets;
  - the ablation values.
- `mc.py` and `mc_reference.py` run the reference classical Monte Carlo for the CrI3 ordering temperature (Metropolis plus over-relaxation, 12×12 to 30×30 cells). `tc_analysis.py` takes its Binder-cumulant crossings.

The lab adds deterministic Gaussian noise (σ = 0.02 meV per supercell energy, seeded by a hash of the request) and an arbitrary constant per magnetic site. Magnetic dipole–dipole interactions are not part of the model (nor of the DFT).

## What is hidden

- The model parameters: all exchange tensors and the single-ion terms.
- The graded targets, computed from the teacher model without noise:
  - the classical ground-state type of each compound (Luttinger–Tisza, confirmed by enumerating collinear patterns);
  - the polar angle of the easy axis;
  - E(x) − E(z) of the rigidly rotated ground-state pattern and its single-ion part;
  - the CrI3 nearest-neighbour bond in the ideal octahedral frame (J, K, Γ; Γ′ is reported but not graded);
  - the classical ordering temperature of the CrI3 model (Binder-cumulant crossing of 24×24 and 30×30 cells);
  - the DM vector of the CrSBr bond from site 0 to site 0 + a1;
  - the zone-centre LSWT gap of CrSBr for S = 3/2;
  - the energies of 10 `to_predict` configurations in 6×6 (honeycomb) and 7×5 (CrSBr) supercells, which the lab does not serve.

## Scope

The lab is a surrogate built from the authors' archived parameters, not new DFT. Anisotropies of 0.01–0.1 meV are at or below the precision the authors quote (0.1 meV). The task tests inference within this model.
