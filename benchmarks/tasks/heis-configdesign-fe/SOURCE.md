# Source and transformations

- **Source:** B. Li, S. von Malottki, G.-M. Rignanese, *Magnetic-configuration design for reliable Heisenberg exchange parameters*, Materials Cloud Archive (2026), record 3p6fv-8hd59, DOI 10.24435/materialscloud:zg-1w.
- **License:** CC BY 4.0.
- **File used:** `Fe.tar.gz` (206,314,109 bytes, md5 3e76b9b6eaeba27b7873be4327806477).

## Construction

The scripts are in `benchmarks/construction/heis-configdesign-fe/`:

- `extract.py` reads the authors' spin-spiral CSV and, for the random states, each OUTCAR (energy(sigma→0) per atom and final local moment vectors) and POSCAR (cell and positions).
- The random-state energies equal the authors' CSV values (OUTCAR energy per atom minus −8274.659 meV/atom) to within 0.002 meV, and this is asserted during extraction.
- No VASP POTCAR or other licensed input is staged.
- `make_task.py` writes the catalog, the lab data and the hidden data.

## What is hidden

- The 20 [214] spin spirals.
- The 20 random states in the 3×3×4 supercell.

Their energies, and reference values derived from the released data, grade the submission:

- **J1:** 8-shell fit to all [001]/[111]/[113]/[121] spirals.
- **D:** small-q limit of the released spirals.

The tolerances were set so that the defensible modelling choices tested in `ablations.py` pass and the plausible wrong paths fail.

## Scope

Retrospective reconstruction from released DFT data with the authors' protocol, not new calculations. The random states and the spirals use different lattice constants in the source (2.829 Å and 2.87 Å).
