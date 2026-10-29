# Source and transformations

- **Source:** D. K. J. Lee and A. R. Natarajan, *Diffusion coefficients of multi-principal element alloys from first principles*, arXiv:2603.24228 (2026), doi:10.48550/arXiv.2603.24228. Data: Materials Cloud Archive record bhv3r-a8q53, doi:10.24435/materialscloud:5g-95, https://archive.materialscloud.org/records/bhv3r-a8q53.
- **License:** CC BY 4.0 (both files).
- **Files used:**
  - `all_data.json` (143,013,225 bytes, md5 15df9136ed668320057131c8c59cad2a): 1,945 single-vacancy nearest-neighbour hops in 4×4×4 bcc supercells (127 atoms + 1 vacancy), each with pymatgen structures and PBE total energies of both end states and the transition state (MACE-MPA-0 relaxations and CI-NEB, then DFT single points).
  - `files_description.md` (395 bytes, md5 8849580fa4e1bb9f72bca8d1e74e571f).

## Construction

The scripts are in `benchmarks/construction/mpea-vacancy-diffusion-kmc/` (`fetch.py`, `extract.py`, `reference.py`, `calibrate.py`, `make_task.py`, `kmclib.py`; compute runs in the runtime image via `run_in_image.sh`).

- `extract.py` maps every atom to the nearest ideal site of the 128-site cell (asserted one-to-one, displacements < 0.25 a), checks that the two end states differ by one nearest-neighbour vacancy–atom exchange and that the transition-state atom nearest the hop midpoint is the hopping atom. All 1,945 entries pass. Energies are copied unchanged.
- 922 entries store their end-state structures with a 14.29 Å cell edge while their transition state has 11.4–12.5 Å; the occupations (fractional coordinates) are unaffected and the energies are as archived. Both cell edges are served as archived.
- 146 entries have the transition-state energy below one end state; they are served unchanged.
- `make_task.py` writes the catalog (hop ids, compositions, hopping element; occupations of the hops to predict), the lab data (served hop records), the kMC service tables, and the private targets.

## What is hidden

- **100 hops** (random, seed 2026, among the 1,764 alloy hops whose saddle lies above both end states): never served. Their archived barriers E_TS − E(end state 1) grade `barrier_predictions_eV`.
- **Diffusion targets:** builders' reconstruction with the released data, not values from the paper. A barrier model (per-element ridge fit of the KRA on six orbits around the hop pair plus quadratic nearest-neighbour terms; three-shell pair Hamiltonian fitted to end-state total energies) trained on all 1,799 clean archived hops drives the lab's kMC engine (`environment/lab/kmc.c`) with equilibrated 10×10×10 cells, 600 trajectories × 50,000 hops, giving log10 η of three alloys, the tracer ratio and fastest tracer of B, the activation energy of the fastest tracer of C (1700–2300 K), and log10 η of the thermo-ideal and KRA-simple variants of B.
- The authors used a separately trained DFT energy cluster expansion (not in this archive) for end-state energies and a neural-network eLCE for the KRA; the end-state energetics here come from the released hop energies, so the thermodynamic contribution is larger than the paper reports and the binary classifications in the paper (for example Mo–W and Cr–W) are not reproduced. None of the paper's headline compositions is graded.

## Scope

Retrospective reconstruction from released DFT data with the authors' kinetic framework (KRA decomposition, rule-of-mixtures enhancement factor, thermo-ideal and KRA-simple variants), not new DFT. Tolerances sit between the spread of defensible analyses (subsets of purchasable hops, model forms, pathological entries kept) and the realistic wrong paths in `solution/ablations.py`.
