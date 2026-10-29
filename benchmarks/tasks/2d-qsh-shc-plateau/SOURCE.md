# Source and transformations

- **Study:** J. Zhou, S. Poncé, J.-C. Charlier, *High-throughput calculations of spin Hall conductivity in non-magnetic 2D materials*, npj 2D Mater. Appl. 9, 39 (2025), DOI 10.1038/s41699-025-00562-4; arXiv:2501.02110.
- **Archive:** Materials Cloud Archive record rfcdy-wcy32, DOI 10.24435/materialscloud:yv-h7, https://archive.materialscloud.org/records/rfcdy-wcy32.
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `Data.7z`: 408,780,659 bytes, md5 1d261e780f092b8c7c1f2e051855c798 (one solid LZMA2 block, downloaded whole);
  - `README.txt`: md5 aedc173ccb0587b13fad1cbbcd5f58cc;
  - `files_description.md`: md5 10c58cb5ed5d38de19f445f2fca21c19.
- **Members:** the 73 members used are listed with their md5 in `benchmarks/construction/2d-qsh-shc-plateau/members.md5`:
  - for each of the eight materials `Wannier_Hamiltonians/insulators/<folder>/` (Bi2_P-3m1, Hf2Br2_P-3m1, CuLi2As_P-6m2, Ti2I2N2_Pmm2, Li2Tl2_P3m1, Bi2In2_P4|nmm, Hg4O2_P1, KSnAs_P3m1): `aiida.in`, `aiida.out`, `aiida.win`, `cut_aiida.wout`, `aiida_band.dat` and `aiida_hr.dat`;
  - `Top_materials/{Bi2,Hf2Br2,CuLi2As,Ti2N2I2}/SHC-{xyx,xyy,xyz}/`: `aiida-shc-fermiscan.dat` and `aiida.win`;
  - `Data/README.txt`.

## Construction

The scripts are in `benchmarks/construction/2d-qsh-shc-plateau/`:

- `fetch.py` downloads and verifies the archive and extracts the members.
- `build_models.py` (runtime image) reads each `aiida_hr.dat`, applies Wannier90's `use_ws_distance` image selection with the mp_grid of `aiida.win` and the final Wannier centres of `cut_aiida.wout` (the archived runs set `use_ws_distance = .true.`), merges images that differ only along the out-of-plane lattice vector (k_z = 0), and writes `environment/lab/data/<material>.npz` (R, H(R) as complex64, Wannier centres). Each model reproduces the archived Wannier bands `aiida_band.dat` to 0.07 meV or better (without the image selection the error is 1–22 meV).
- `reference.py` (runtime image) computes the references and the ablation values with the lab's own evaluator (`qsh_tb.py`) on the lab's npz files.
- `make_task.py` (runtime image) stages, for each material, `aiida.in` and `aiida.win` unchanged, `aiida.out.excerpt` (`aiida.out` without the band energies and occupations of the final scf iteration) and `cut_aiida.wout.excerpt` (`cut_aiida.wout` without the disentanglement and minimisation iteration lines: header and final state) in `environment/agent/archive/<material>/` (with `|` replaced by `_` in the folder name), together with the archive README; the excerpts avoid unrelated numbers in the logs that coincide with graded values at 4 significant digits. It also writes `catalog.json`, `environment/lab/data/meta.json`, `environment/lab/data/archived_shc.json` (the twelve archived Top_materials spin Hall Fermi scans and their `aiida.win`, values exactly as in the files), the private targets and `solution/ablation_values.json`.
- Apart from the two excerpts, no file is altered in content. No pseudopotentials are in the archive or in the task.

## What is hidden

Graded values are derived quantities stored only in `environment/lab/private/hidden.json`:

- **Z2** of the occupied Wannier bands of each model, from 201 × 300 Wilson loops along both reciprocal axes (both agree).
- **Plateau σ^z_xy** of each model in (ħ/e)(e²/h): the lab evaluator at T = 0, η = 10⁻⁴ eV, with the Wannier-centre position operator, on a 180 × 180 mesh with 8 × 8 adaptive refinement of points with |Ω| > 10 Å², at the mid-gap of the 300 × 300 band edges of the model.
- **Spin direction** (θ, φ) of the plateau of Bi2_P-3m1 and Hg4O2_P1 from the same calculation.

The lab never returns these values directly. The archived Top_materials files (relaxed structures) and the paper's Table 2 give plateaus of 0.53/0.5 (Bi2), −0.91/−0.9 (Hf2Br2), 1.03/1.0 (CuLi2As) and −0.95/−1.0 (Ti2N2I2); they do not describe the database models, and those for Bi2, Hf2Br2 and Ti2N2I2 fail the plateau tolerance.
