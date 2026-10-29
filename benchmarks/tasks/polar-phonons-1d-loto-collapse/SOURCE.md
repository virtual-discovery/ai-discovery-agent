# Source and transformations

- **Study:** N. Rivano, N. Marzari, T. Sohier, *Infrared-active phonons in one-dimensional materials and their spectroscopic signatures*, npj Comput. Mater. 9, 194 (2023), DOI 10.1038/s41524-023-01140-2.
- **Archive:** Materials Cloud Archive record 08gsg-gpd20, DOI 10.24435/materialscloud:46-wj, https://archive.materialscloud.org/records/08gsg-gpd20.
- **License:** CC BY 4.0 (whole record).
- **Files used:**
  - `materialscloud_data_Rivano.zip` (30,776,903 bytes, md5 f3ee3048a131eb831eecee865d51e705);
  - `files_description.md` (280 bytes, md5 0c399e8707d568e68028464c7cb683b0).

  The sha256 of every unpacked member (270 files) is recorded in `benchmarks/construction/polar-phonons-1d-loto-collapse/members.json`.

## Construction

The scripts are in `benchmarks/construction/polar-phonons-1d-loto-collapse/`:

- `fetch.py` downloads the two files through the archive's per-file links, checks the md5 and unpacks the zip (macOS metadata skipped).
- `make_task.py` (run in the runtime image) writes the catalog, the lab data, the private targets and the ablation values. It first checks that its numpy implementation of the authors' 1D long-range interpolation reproduces the archived 1D interpolated dispersions of BN(5,5), BN(6,6), the BN chain and the GaAs nanowire (LO within 0.03, 0.15, 0.50 and 0.52 cm⁻¹).

## Transformations

- **Dynamical-matrix files** are served verbatim (gzip in the lab). Byte-identical copies in the archive are listed once, under the first path, with the other paths given as `also_archived_as`.
- **Interpolated dispersions:** the archived `matdyn.modes` files (6–14 MB) are reduced to frequencies and, for each mode, the axial share of the normalized displacement pattern. The BN-chain path (2144 points) keeps every point up to 0.3 |ΓZ| and every tenth point beyond.
- **Program outputs** are served verbatim up to their first 400,000 characters. Only the BN-chain 1D `ph.out` (8.3 MB) and `matdyn.out` (1.2 MB) are cut.
- **Input files** are copied into the catalog with the `ATOMIC_POSITIONS` blocks of the pw.x inputs omitted; the positions are in every dyn-file header.
- **Truth:** the LO branch is the mode with the largest longitudinal mode effective charge, computed with the Born charges of the Γ run of the same folder.

## What is hidden

- **BN(5,5), 1D boundary conditions:** the explicit small-q runs `longwave/dyn100`, `dyn25`, `dyn50`, `dyn75` and `dyn125` (q/|ΓZ| = 0.01, 0.025, 0.05, 0.075, 0.125), and `dyn10`, an archived copy of `dyn100`.
- **BN(6,6), 1D boundary conditions:** `longwave/dyn100`, `dyn25`, `dyn50`, `dyn75`, `dyn10` and `dyn125` (0.01–0.125), and their copies in the top-level `1D_code` folder.
- **Both tubes, both boundary conditions:** the archived `matdyn.modes` and `matdyn.out`.

The LO frequencies of the withheld explicit runs grade the 1D fields. The q → 0 point of the archived BN(5,5) 3D-periodic `matdyn.modes` grades the 3D field.

## Archive features found during construction

- The BN(4,4) `1D_code/matdyn.modes`, despite `loto_2d = .true.` in its `matdyn.in`, equals a plain Fourier interpolation without the 1D long-range term (within 1.2 cm⁻¹). It lies 33 cm⁻¹ below the tube's own explicit 1D runs at 0.05 |ΓZ|.
- BN(5,5) `longwave/dyn10` is a copy of `dyn100` (q = 0.01 |ΓZ|), so the archive has no 0.10 |ΓZ| run for this tube.
- `GaAs_bulk/ZB_phase/dyn1` is a copy of the wurtzite file. The zincblende results are only in its `ph.out`.
- The `GaAs_nanowire/3D_code` folder has no `dyn0`.
- Several archived titles name the wrong system: the BN(6,6) files and the 24-atom GaAs nanowire are titled "GaAs nanowire wz 84".

## Scope

This is a retrospective analysis of released DFPT data; no new calculation is run. The truth is the authors' explicit DFPT runs, not their interpolation.
