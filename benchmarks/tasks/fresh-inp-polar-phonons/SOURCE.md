# Source and transformations

This is a **fresh-start** task: the lab runs new Quantum ESPRESSO 6.7 calculations (pw.x and ph.x) that the agent requests. No archived calculation results are used or served. The only external inputs are two pseudopotential files.

## Pseudopotentials

- **Source:** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018); SSSP library, Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym (https://archive.materialscloud.org/records/rcyfm-68h65).
- **Archive files used:**

  | file | bytes | md5 |
  |---|---|---|
  | `SSSP_1.3.0_PBE_precision.tar.gz` | 62,963,841 | fde94756886f32ada7bf597547557eb5 |
  | `SSSP_1.3.0_PBE_precision.json` | 22,485 | 1692c5c9ce89e1c7c783f8f0eee0cbfa |
  | `LICENSE.txt` | 3,673 | c0f9df05fb4cd4921f4caa25606815ef |

- **Files vendored unchanged into `environment/lab/`** (the agent image has neither QE nor these files):

  | file | SSSP family | md5 (matches SSSP metadata) | sha256 | license |
  |---|---|---|---|---|
  | `In.pbe-dn-rrkjus_psl.0.2.2.UPF` | 031US (PSlibrary, A. Dal Corso) | 25e7c42c3b55b68f4bf926be2e7201a4 | ac76860606654bf63d2950129329ee50267d50ad92296728a7d42e9355bcfb39 | GNU GPL v2 or later |
  | `P.pbe-n-rrkjus_psl.1.0.0.UPF` | 100US (PSlibrary 1.0.0, A. Dal Corso) | 8b930f418f0a4573dca56cd030ffe088 | 2d112dfec2e2d9b75a971574d9116aa92249988177791659bbcd2ba15f23c20c | GNU GPL v2 or later |

- **License:** the archive record is CC BY 4.0, except the pseudopotentials, which keep their family licenses as listed in the archive's `LICENSE.txt` (vendored as `environment/lab/SSSP_LICENSE.txt`): PSlibrary files are GPL v2 or later (reference: A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014)). The UPF files are redistributed unmodified, with their embedded generation information.
- **Recommended cutoffs (SSSP metadata):** In 50/400 Ry, P 30/240 Ry. The lab uses ecutwfc 50 Ry (optionally 65) with ecutrho = 8 × ecutwfc.

## Construction

Scripts are in `benchmarks/construction/fresh-inp-polar-phonons/`:

- `fetch_pseudos.py` downloads the archive files, checks their md5, extracts the two UPF files, checks them against the SSSP metadata md5 and the pinned sha256 (`pseudo_checksums.json`), and installs them into the lab.
- `run_request.sh` runs one lab request through `environment/lab/task_lab.run` (which writes the pw.x/ph.x inputs in `qe_driver.py`) in the runtime image on the builder's cores.
- `reference_runs.sh` lists every reference calculation; `reference_results.json` keeps their parsed results.
- `make_task.py` computes the reference answers, the ablation and robustness values and the calibration table from those results, with the oracle's own analysis functions.

## What is hidden

`environment/lab/private/hidden.json` holds the reference values of every graded field and the tolerances. They come from the construction calculations above, made through the same `task_lab.run` code path the agent uses: the Birch–Murnaghan a0 of 6×6×6-shifted SCF energies; the TO and force constants from 6×6×6-shifted DFPT at a0; ε∞ and Z* from the mean of the 8×8×8- and 10×10×10-shifted electric-field responses; centred ±1 % differences in a for the volume derivatives; and 6×6×6-shifted DFPT of the [001]-strained cell. The lab never serves these values; the agent recomputes them.

## Scope

New PBE calculations with one pseudopotential set and one code version; protocol-specific values, not experimental or literature numbers.
