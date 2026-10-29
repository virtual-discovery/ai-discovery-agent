# Source and transformations

- **Source:** S. Calderon, J. Hayden, S. Baksa, W. Tzou, S. Trolier-McKinstry, I. Dabo, J.-P. Maria, E. Dickey, *Atomic-scale polarization switching in wurtzite ferroelectrics*, Science 380, 1034 (2023), DOI 10.1126/science.adh7670.
- **Data:** "Atomic-scale polarization switching in wurtzite ferroelectrics - DFT data", Zenodo record 7916565, DOI 10.5281/zenodo.7916565, https://zenodo.org/records/7916565.
- **License:** CC BY 4.0.
- **File used:** `AlBN-system-Science.zip` (16,033,027 bytes, md5 814d04e76cfb237221c961154c5e8716, sha256 a9b32ffc6b7af09de803552d480ce2af7ef478abd25df005d1ac842fef4e652e).

## Construction

The scripts are in `benchmarks/construction/wurtzite-ferroelectric-switching/`:

- `fetch.sh` downloads and unpacks the archive.
- `extract.py` parses, for AlN/ and AlBN/:
  - each neb.x run in `NEB/{NEB,refined_NEB,refined_NEB2}`: the path settings and the last energy/error table of `AlBN.neb.out`, the activation-energy lines, the convergence message, the log tail and `scheduler.out`, the engine settings from `neb.in`, and the final image structures from `AlBN.axsf`;
  - each `polarization/nscf/NN_image` Berry-phase run: the POLARIZATION CALCULATION section of the nscf output verbatim, the structure and settings from `nscf.in` and `scf.in`, and the scf total energy;
  - four vc-relax outputs, and the Born-charge CSV.
- The extraction asserts that:
  - the polarization runs are the refined_NEB2 images (same atoms; scf energies agree within 2 meV);
  - the raw polarizations equal the authors' `*_polarization.txt` column.
- `make_task.py` writes the catalog, the lab data, the private targets and the ablation values.

The authors' analysis files are not served: the notebook, `energies*.txt`, the "Adjusted Polarization" column, the plots and the GIFs. The same holds for the pseudopotential-dependent inputs other than their names.

## What is hidden

The graded values are derived from the archived outputs:

- the barriers per f.u. along the converged minimum-energy paths and the (Al,B)N collective path;
- the spontaneous polarizations relative to the layered-hexagonal reference (continuous branch over all images);
- the polarization change at the (Al,B)N saddle;
- the energy and inverted-cation count of the lowest (Al,B)N intermediate;
- the switching mode.

The lab serves only raw run outputs.

## Scope

This is a retrospective analysis of released calculations, not new DFT. A zinc-blende reference calculation is not part of the archive and is not offered.
