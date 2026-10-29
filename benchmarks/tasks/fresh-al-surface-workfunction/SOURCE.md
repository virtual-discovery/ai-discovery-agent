# Source and transformations

This is a fresh-start task: the lab runs new Quantum ESPRESSO 6.7 calculations requested by the agent. No archived energies, potentials or reference values are served or used.

## Pseudopotential

- **Citation:** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018). SSSP library v1.3.0, PBE precision.
- **Archive:** Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym, https://archive.materialscloud.org/records/rcyfm-68h65.
- **Files fetched:**
  - `SSSP_1.3.0_PBE_precision.tar.gz` (md5 fde94756886f32ada7bf597547557eb5);
  - `SSSP_1.3.0_PBE_precision.json` (md5 1692c5c9ce89e1c7c783f8f0eee0cbfa);
  - `LICENSE.txt` (md5 c0f9df05fb4cd4921f4caa25606815ef).
- **File used:** `Al.pbe-n-kjpaw_psl.1.0.0.UPF`, PSLibrary 1.0.0 PAW, recommended cutoffs 30/240 Ry.
  - 1,762,723 bytes; md5 cfc449ca30b5f3223ec38ddd88ac046d (matches the SSSP metadata); SHA256 fd7b78921e6d0939095b681c328732668eaefd1dee6882fdbc03be463550cc97.
  - It is the same file as in the v1 task `t1-al-fcc-eos-fresh`.
- **License:** the archive is CC BY 4.0, except for pseudopotential families with their own licenses. Its `LICENSE.txt` gives PSLibrary 1.0.0 files under the GNU General Public License, version 2 or later. The UPF file is vendored unmodified into `environment/lab/` (lab image only, not the agent image), together with the record's license text as `environment/lab/SSSP_LICENSE.txt`.
- **Refetch:** `benchmarks/construction/fresh-al-surface-workfunction/fetch_pseudo.py` downloads the files, checks every checksum above and installs the UPF.

## Construction

The scripts are in `benchmarks/construction/fresh-al-surface-workfunction/`.

- `environment/lab/slab_qe.py` writes every pw.x and pp.x input; `task_lab.run` calls it. The settings are PBE, ecutwfc 30 Ry, ecutrho 240 Ry, an fcc lattice constant of 4.04 Å, Marzari–Vanderbilt smearing, local-TF mixing and Γ-centred meshes.
  - The lattice constant is the construction's own PBE equilibrium value, 4.039 Å, rounded; it was found from a 7-point E(a) fit on a 24³ mesh.
  - Planar averages are taken over the surface cell from the pp.x cube files: `plot_num=11` gives the electrostatic potential and `plot_num=1` the total local potential.
- `qe_batch.py` runs lab requests through that same code, in the runtime image, on the builder's cores. The requests are listed in `reference_batch.json`, and `run_reference.sh` reruns them.
- `make_task.py` reduces the runs to the private reference, the ablation values and the calibration table. `reference_runs.json.gz` keeps the parsed results, so the reduction can be redone without QE.

## What is hidden

The graded values come only from these construction calculations and are stored only in `environment/lab/private/hidden.json`:

- **Surface energies and work functions:** for Al(111) and Al(100), the averages over bulk-terminated slabs of 7–14 and 18 layers. They use the 18×18 mesh, MV smearing of 0.02 Ry, a 15 Å vacuum and enough bands (2.5 N + 4). σ uses an oriented bulk reference sampled like the slab (kz = 24 for (111), 27 for (100)). φ uses the electrostatic vacuum plateau.
- **Thick slabs:** the free energies and work functions of the 18-layer slabs of both faces. These are beyond the lab's 12-layer limit on that mesh and are never served.
- **Relaxation:** Δd12 of Al(111), the mean over symmetric relaxations of the 9–12-layer slabs.
- **Ordering:** the work-function ordering.

## Scope

These are new calculations with one fixed protocol. They do not reproduce published slab studies, and the graded numbers are specific to this protocol, not literature values.
