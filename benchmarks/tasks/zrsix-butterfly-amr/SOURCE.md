# Source and transformations

- **Source:** S. Zhang and O. V. Yazyev, *Magnetoresistance in ZrSiX (X = S, Se, Te) nodal-line semimetals*, Phys. Rev. B 112, 165147 (2025), arXiv:2504.12547, DOI 10.1103/914y-4gp9.
- **Data:** Materials Cloud Archive record tyan3-q3j67 (2025-10-16), DOI 10.24435/materialscloud:3q-r9, https://archive.materialscloud.org/records/tyan3-q3j67.
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `ZrSiX data.zip` (8,378,201 bytes, md5 cbc8f3896f69ebf7c1cf68ad9f1a1dc8), 32 members verified against `benchmarks/construction/zrsix-butterfly-amr/members.md5`;
  - `README.txt` (1,089 bytes, md5 199687afb593272f9f73fc873188a519).

## What the agent sees

- `/app/archive/`: every member of the zip except the three Wannier Hamiltonians and the three VASP POTCARs (licensed; never staged), unchanged and under their own folder names, plus the record README (`README_record.txt`).
- `/app/archive/input_MR/<compound>/<hr file>.header`: the header of each Hamiltonian (title line, number of orbitals and of lattice vectors, the degeneracy list) and its first matrix-element line.

## Lab evaluators (new calculations on the archived Hamiltonians)

The scripts are in `benchmarks/construction/zrsix-butterfly-amr/` (see its README):

- `make_model.py` converts the three `wannier90_hr` files to npz (H(R) divided by the Wigner–Seitz degeneracy). The lab stores them as float32 (`data/model_*.npz`) and evaluates band energies, Hellmann–Feynman velocities, Fermi-surface slices and semiclassical orbits exactly.
- `make_counts.py` histograms the exact band energies of every model on the 40, 80, 120 and 160 meshes (0.2 meV bins) for the electron-count requests.
- `make_spline.py` and `make_cache.py` build the Chambers tables: periodic cubic splines of the characteristic-polynomial coefficients of the bands that cross the Fermi level (128 × 128 × 64 grid of exact eigenvalues), RK4 orbits for every start point of the 160 × 160 × 80 mesh inside an energy window, orbit averages of v_i v̄_j for 15 values of Bτ, binned by energy and band. `validate_exact.py` and `test_parabolic.py` check the integrator against exact diagonalization and the analytic Drude result.
- Lattices: the POSCAR / `wannier90.win` cell of each compound (ZrSiS a = 3.5496 Å, c = 8.1637 Å; ZrSiSe 3.623 / 8.365 Å; ZrSiTe 3.6974 / 9.508 Å).

## What is hidden

All graded values are computed with the lab's own code path from the archived Hamiltonians, under the conventions of the instruction, and are stored only in `environment/lab/private/hidden.json`. The paper's figures are not used as targets.
