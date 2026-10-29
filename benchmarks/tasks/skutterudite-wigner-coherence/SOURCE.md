# Source and transformations

- **Source:** E. Di Lucente, M. Simoncelli and N. Marzari, *Crossover from Boltzmann to Wigner thermal transport in thermoelectric skutterudites*, Phys. Rev. Research 5, 033125 (2023), DOI 10.1103/PhysRevResearch.5.033125. Data: Materials Cloud Archive record 5jsfj-cdg70 (2024-01-05), DOI 10.24435/materialscloud:xn-qr, https://archive.materialscloud.org/records/5jsfj-cdg70.
- **License:** CC BY 4.0 (all files).
- **Archive file:** `materials_cloud_archive_submission.zip` (1,213,559,335 bytes, md5 9dc2bd3b8659022e5703862f82adaa40). Only 87 members (0.61 GB) were read, through HTTP range requests on the zip's central directory and single members. Each member was checked against its zip CRC-32. The sha256 of every member is in `benchmarks/construction/skutterudite-wigner-coherence/members.json`.

## Main files used

| Member (under `materials_cloud_archive_submission/`) | Bytes | sha256 |
|---|---|---|
| `Co_based/CoSb3/5_phono3py/tc_calc_no_isotope/kappa-m151515.hdf5` | 2637974 | 92998e76b5ce0885922a7bf410d2a237ad2f8e0342db309235213889b662dacc |
| `Co_based/CoSb3/5_phono3py/tc_calc_with_isotope/kappa-m151515.hdf5` | 2637974 | 92998e76b5ce0885922a7bf410d2a237ad2f8e0342db309235213889b662dacc |
| `Co_based/CoSb3/5_phono3py/v_operator_mean_free_path_isotope/kappa-m151515.hdf5` | 260873300 | 6ef7b6dcb79da01e7f167c2e125f5446904114c21f7117e021fbbdc52ff27062 |
| `Co_based/YbCo4Sb12/3_hiphive_FCS/5_phono3py/kappa-m999.hdf5` | 71603500 | 30b2054d10db0d363e8bcbfd33687e9fc3b18a6effaef2f543a54908c08fd932 |
| `Fe_based/FeSb3/4_hiphive_FCS/5_phono3py/kappa-m999.hdf5` | 66902432 | e1527b5bd07f15ba1fd568c13376679f1993377d26f33e72d92101a5ec682df2 |
| `Fe_based/YbFe4Sb12/NM/6_hiphive/5_phono3py/fc2hip_fc3hip/kappa-m999.hdf5` | 63625867 | 1d80e69626378bc75672f3055112bfad1e0f2aea85adccce98fb8ae65f9a6184 |
| `Fe_based/YbFe4Sb12/NM/3_FCS_thirdorder/5_phono3py/no_isotope/kappa-m999.hdf5` | 745651 | c66da7108c976e7eca433c35dd5d80588f3560fda1d4e78f8dfeb992d37ce4b9 |
| `Fe_based/YbFe4Sb12/NM/3_FCS_thirdorder/5_phono3py/isotope/kappa-m999.hdf5` | 745651 | c66da7108c976e7eca433c35dd5d80588f3560fda1d4e78f8dfeb992d37ce4b9 |
| `Ir_based/IrSb3/3_hiphive_FCS/5_phono3py/kappa-m999.hdf5` | 57678105 | a8dd97748fa48f0e14475a5ea26c69e5cd72f086bfaf659991960ef5a09de77b |
| `Ir_based/YbIr4Sb12/3_hiphive_FCS/5_phono3py/kappa-m999.hdf5` | 75207170 | 0816a8c5b63e5147c34b297dbbb2e522b9bf93c2d220a3a4531b8c7d2df4ca9f |

The other members are each run's `phono3py.yaml`, launch script and log (`tc_*.out`, `nohup.out`), the CoSb3 post-processing files `tc_data.dat` and `2_plot_MFP.py`, the phonopy `projected_dos.dat`, `phonopy.yaml` and `thermal_displacements.yaml` of the six compounds, the MSD `POSCAR`s and commands, the PDOS `cmd.sh` and `settings.conf`, the YbFe4Sb12 NM and AFM DFPT `matdyn__crystal.freq` files, and `README.txt`.

## Transformations

The scripts are in `benchmarks/construction/skutterudite-wigner-coherence/`:

- `fetch.py` refetches and verifies the members.
- `extract.py` (needs h5py, so it runs on the host) reduces each phono3py file to compact arrays: frequency, gamma(T), gamma_isotope, heat_capacity, q-points and weights. It stores velocities in two forms:
  - for the phono3py 2.1.0 runs, (1/3) Σ_a |V_a(q)_{ss'}|² per irreducible q, which is the trace of the archived star-summed `gv_by_gv_operator` divided by the q-point weight;
  - for the 1.19 runs, (1/3)|v_s(q)|² from the archived group velocities.

  It also parses the K_P, K_C and K_T tables and headers from the logs, and sums the PDOS and MSD over the atoms of each element.
- `make_task.py` writes the lab data (per-mode arrays and velocity operators, without the archived κ arrays), the catalog, the private reference and the ablation values. It first asserts two things:
  - the reference evaluator (`model.py`) reproduces every archived `kappa_P_RTA` and `kappa_C` of the 2.1.0 runs to 1e-6 from the archived inputs;
  - the 1.19 per-mode files reproduce their log K_P.

The served arrays are rounded to 7 significant digits (6 for the velocity operator).

## What is hidden

- The archived outputs of `IrSb3/hiphive` and `YbIr4Sb12/hiphive`: the κ arrays of their kappa files and their logs. Their per-mode arrays and velocity operators are served.
- The graded reference values are computed from those archived inputs with the authors' formula, except that the 21 corrupted `gamma_isotope` entries of YbIr4Sb12 are set to zero. These entries are 1.7×10³ to 2.0×10⁷ THz, all at the q-point nearest Γ; every other entry is below 0.01 THz. The graded values are:
  - κ_P of YbIr4Sb12 at 300 and 700 K;
  - κ_C of YbIr4Sb12 at 100, 300 and 700 K;
  - T_half of YbIr4Sb12;
  - κ_P of IrSb3 at 50 K and κ_TOT of IrSb3 at 300 K.
- With the corrupted entries kept, the evaluator gives the authors' archived values, for example κ_P(YbIr4Sb12, 300 K) = 0.816 W/mK. See the construction README.

## Scope

Retrospective analysis of released phono3py data; no new phonon calculations. The archive's own inconsistencies are served as archived. These include:

- byte-identical "isotope"/"no isotope" kappa files of the 1.19 runs;
- tc_data.dat, which is inconsistent with the phono3py log in the same folder;
- the YbFe4Sb12 hiPhive run, which uses hiPhive rather than DFPT second-order constants.
