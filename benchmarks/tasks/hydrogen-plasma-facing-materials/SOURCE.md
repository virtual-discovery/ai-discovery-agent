# Source and transformations

- **Source:** A. Fedrigucci, N. Marzari, P. Ricci, *DFT calculations of surface binding and interstitial hydrogen formation energies for plasma-facing materials*, Materials Cloud Archive record gpxwh-dx645, DOI 10.24435/materialscloud:9k-09. The record accompanies A. Fedrigucci, N. Marzari, P. Ricci, *Comprehensive screening of plasma-facing materials for nuclear fusion*, PRX Energy 3, 043002 (2024), arXiv:2407.00858.
- **Archive URL:** https://archive.materialscloud.org/records/gpxwh-dx645
- **License:** CC BY 4.0 (all files).

## Files used

| File | Bytes | md5 (archive) | sha256 |
|---|---|---|---|
| `PFMs_screening_QE_data.tar.xz` | 46,921,676 | 56aadb3e22b5a0d52ff75d5757798aaa | ec5ff80c8c6b43ebd817e7c6386bb71946bf8dde7359d746e6aee34f932274ef |
| `README.md` | 2,514 | 61db996d231cb32b2fbf19328de17307 | ac01105ab9f2124ee57889a729238e35a4feb29e0e87de4ca48ee31a41400ef4 |
| `files_description.md` | 276 | 21acb466968792b84ff2903d878da8c5 | fc6a621f6bad31bc8044efff1b4530e8378c88abbc461766ccdb1df29375988b |

The tarball holds 490 Quantum ESPRESSO runs (pwscf.in / pwscf.out, 765 MB unpacked) for 21 materials.

## Transformations

The scripts are in `benchmarks/construction/hydrogen-plasma-facing-materials/`.

- `fetch.py` downloads the three files and checks their md5.
- `extract.py` parses every run without recomputation. It reads:
  - the input settings (calculation type, cutoffs, k-points, smearing, spin, thresholds, pseudopotential file names);
  - the cell and starting positions;
  - the run header (code version, start time, MPI processes, wall time);
  - the outputs (last SCF and final energies in Ry, number of SCF cycles, relaxation message, final total force, magnetization, Fermi energy, pressure, final positions and final cell).
- `make_task.py` writes:
  - the catalog, with the archive README quoted verbatim;
  - the lab data: one record per run, exactly as parsed;
  - the private references and the ablation values.
- No pseudopotential file or raw output text is staged.

Runs are served under catalog ids derived from the archive folders, for example `SBE/WC-h-187/vacancy/100_t2/W_0_0.0`.

## What is hidden

Nothing in the archive is withheld from purchase. The lab never serves derived quantities; the graded values are computed from the archived runs with the conventions stated in `instruction.md`:

- the hydrogen formation energies (½H₂ reservoir) at the preferred site, for graphite, HfB₂, Mo, Re, TaC, TiN and WC;
- the distinct-site energy lists and the preferred-site geometry of the same hosts;
- the minimum surface binding energies on the lowest-surface-energy termination, for W, WC and TiB₂;
- the cohesive energies of Re, W and WC.

The paper tabulates formation energies at the preferred site (rounded to 0.01 eV) and SBEs averaged over sites and terminations. It contains no distinct-site lists, site geometries, minimum SBEs or cohesive energies.

## Notes on the source data

- **The published compound SBEs mix cutoffs.** The paper's Table V values for W in WC, W₂C and TaW, and for Ta in TaN, are reproduced only by using the W (30/240 Ry) and Ta (45/360 Ry) atoms of the elemental sets. The archive's own WC, W₂C and TaW folders contain W atoms at 45/360 Ry, and the TaN folder a Ta atom at 60/480 Ry. The references here use each material's own folder.
- **The TiB₂ Ti atom is a copy.** The archived TiB₂ Ti atom is byte-identical to the TiN Ti atom (60/480 Ry), while the TiB₂ slabs and bulk use 35/280 Ry. A construction-time check (`cutoff_check/`) puts the Ti atom energy change between 35 and 60 Ry at the meV level. Both the number and `null` are therefore accepted for that entry.

## Scope

This is a retrospective reconstruction from released DFT data with the authors' protocol, not new calculations. The construction-time cutoff check is the only exception, and it is used only to justify a grading choice.
