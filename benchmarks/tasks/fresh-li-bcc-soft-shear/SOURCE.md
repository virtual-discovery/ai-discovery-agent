# Source and transformations

This is a **fresh-start** task: the lab runs new Quantum ESPRESSO 6.7 pw.x calculations that the agent requests. No archived calculation results are used or served. The only external input is one pseudopotential file.

## Pseudopotential

- **Source:** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018); SSSP library v1.3.0 (PBE precision), Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym (https://archive.materialscloud.org/records/rcyfm-68h65).
- **Archive files used:**

  | file | bytes | md5 |
  |---|---|---|
  | `SSSP_1.3.0_PBE_precision.tar.gz` | 62,963,841 | fde94756886f32ada7bf597547557eb5 |
  | `SSSP_1.3.0_PBE_precision.json` | 22,485 | 1692c5c9ce89e1c7c783f8f0eee0cbfa |
  | `LICENSE.txt` | 3,673 | c0f9df05fb4cd4921f4caa25606815ef |

- **File vendored unchanged into `environment/lab/`** (the agent image has neither QE nor this file):

  | file | SSSP family | md5 (matches SSSP metadata) | sha256 | license |
  |---|---|---|---|---|
  | `li_pbe_v1.4.uspp.F.UPF` (382,791 bytes) | GBRV 1.4 ultrasoft (K. F. Garrity, J. W. Bennett, K. M. Rabe, D. Vanderbilt, Comput. Mater. Sci. 81, 446 (2014)) | e912e257baa3777c20ea3d68f190483c | 02cc4b3810e28a43b570277840d9631ca24f634d922928200c047b72c0a86bf5 | GNU GPL v3 |

- **License:** the archive record is CC BY 4.0, except the pseudopotentials, which keep their family licenses as listed in the archive's `LICENSE.txt` (vendored as `environment/lab/SSSP_LICENSE.txt`): GBRV files are GPL v3. The UPF file is redistributed unmodified, with its embedded generation information.
- **Recommended cutoffs (SSSP metadata):** 40/320 Ry. The lab uses ecutwfc 40 Ry (optionally 60) with ecutrho = 8 × ecutwfc.

## Construction

Scripts are in `benchmarks/construction/fresh-li-bcc-soft-shear/`:

- `fetch_pseudos.py` downloads the archive files, checks their md5, extracts the Li UPF, checks it against the SSSP metadata md5 and the pinned sha256 (`pseudo_checksums.json`), and installs it (and the archive licence) into the lab.
- `reference_jobs.py` lists every reference calculation as a lab request; `run_batch.sh` runs them through `environment/lab/task_lab.run` (which writes the pw.x input in `qe_driver.py`) in the runtime image on the builder's cores; `reference_results.json` keeps their parsed results.
- `make_task.py` computes the reference answers, the ablation and robustness values and the calibration table from those results, with the oracle's own analysis functions (`solution/solve.py`).

## What is hidden

`environment/lab/private/hidden.json` holds the reference values of every graded field and the tolerances. They come from the construction calculations above, made through the same `task_lab.run` code path the agent uses (some at denser meshes than the lab offers, to confirm convergence): a0 and B from a Birch–Murnaghan fit of k48 / 0.01 Ry free energies; C′ from symmetric ±1 % volume-conserving tetragonal stress differences at a0 averaged over the converged meshes and smearing widths; C44 from symmetric ±1 % trigonal shear stresses; the compressed-cell pressure and stress–strain coefficients from the same strains about the a = 3.30 Å cell. The lab never serves these values; the agent recomputes them.

## Scope

New PBE calculations with one pseudopotential and one code version; protocol-specific values, not experimental or literature numbers.
