# Construction of heis-configdesign-fe

Reference task for the benchmark (see `benchmarks/tools/BUILDING.md`).

## Rebuild

1. Download `Fe.tar.gz` from https://archive.materialscloud.org/records/3p6fv-8hd59 (md5 3e76b9b6eaeba27b7873be4327806477) into `downloads/heis-configdesign/` and unpack it there.
2. `python3 benchmarks/construction/heis-configdesign-fe/extract.py` writes `pool.json` (101 spirals, 89 random states, DFT energies and moments).
3. `cd benchmarks/construction/heis-configdesign-fe && python3 -B make_task.py` writes the catalog, lab data, private targets and ablation values, and prints the calibration numbers.
4. `.venv/bin/python -m benchmarks.tools.check_task heis-configdesign-fe` writes `check.json`.

`model.py` holds the Heisenberg utilities; `j_reference.json` is the all-data 8-shell fit.

## Build notes

- **Curie temperature dropped as a target.** Across equally valid fits the classical RPA Tc ranges 205–915 K (sometimes negative), so it is not uniquely recoverable. It remains a rubric point (`model.py` keeps the mean-field and RPA estimators used for that check).
- **Two lattice constants.** The random states were computed at a = 2.829 Å and the spirals at 2.87 Å. Shell radii are assigned in each cell's own lattice; using the spiral radii finds no neighbours (ablation dp6).
- **J1 window [18.0, 24.5] meV.** Over budgeted purchase plans fitted with 6–10 shells, J1 spans 19.0–23.5 (5th–95th percentile). The nearest wrong paths are 28.0 (large-angle model) and 5.2 (flat-spiral convention).
- **Oracle.** Buys 9 q-points along each of [111], [113] and [121] plus 2 small-q [001] points, and the 2x3x6 and 3x4x1 random states (58 of 60 units). It chooses the shell count by leave-one-direction-out validation (9 shells).
- **Calibration table:** `design.json` → `calibration`.
- **Gate (2026-10-19):** passed, with all 10 ablations failing. The Harbor oracle scored reward 1.0.
