# Source and transformations

- **Source:** S. Grillo, C. Cignarella, F. Bechstedt, P. Gori, M. Palummo, D. Campi, *Quasiparticle effects and strong excitonic features in exfoliable 1D semiconducting materials*, ACS Nano 20, 2664–2677 (2026), DOI 10.1021/acsnano.5c14061.
- **Archive:** Materials Cloud Archive record 42t0v-jam38, DOI 10.24435/materialscloud:xq-dy (published 2026-02-05), https://archive.materialscloud.org/records/42t0v-jam38.
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `1D_ACS_Nano_Exciton.zip` (3,319,010 bytes, md5 a3cdbff7bbca1465fa0d2a8ec75ccec4);
  - `README.txt` (985 bytes, md5 781fc14a3db2c37c787b60dac65d45cf).

  The md5 checksums of the 45 members read are in `benchmarks/construction/gwbse-1d-chain-excitons/members.md5`. For each chain (S3, Te3, As2S3, Bi2Te3) they are:
  - `GW/<chain>/o-G0W0.qp` and `o-G3W3.qp` (`o-G4W4.qp` for As2S3);
  - `GW/<chain>/o-BSE_G0W0.alpha_q1_haydock_bse` and `o-BSE_scGW.alpha_q1_haydock_bse`;
  - `GW/<chain>/<chain>_soc.bands.dat.gnu` and `<chain>_nosoc.bands.dat.gnu`;
  - `IR/<chain>/ph.epsil.out` and `IR/<chain>/scf.in`;
  - `volume/<chain>/scf.environ.out`, `scf.in` and `environ.in`;
  - plus the archive's `README.txt`.

## Construction

The scripts are in `benchmarks/construction/gwbse-1d-chain-excitons/`.

- **`fetch.py`** downloads both files through the archive's per-file links, checks them and extracts and checks the members.
- **`model.py`** reads the members and holds the estimators.
- **`make_task.py`** writes the catalog, the lab data, the private targets and the ablation and robustness values.
- **`record_oracle.py`** and **`write_design.py`** record the gate's oracle run and write `design.json` with the calibration table.

The lab serves:

- **QP tables:** the numerical rows unchanged (6 decimals), grouped by k-point, and the comment header without the ASCII-art banner.
- **BSE spectra:** all seven archived columns at 6 significant digits, the column names, and the comment header without the banner.
- **Band structures:** the k coordinate converted from units of 2π/alat to the fraction of 2π/c (the path ends at Z), and the 12 highest occupied and 12 lowest empty bands of the SOC files (6 and 6 for the SOC-free files), unchanged at 4 decimals.
  - The index of the highest occupied band was deduced by the builders: for the SOC files it equals the valence-band maximum of the QP tables, whose Eo gaps the SOC files reproduce within 4 meV on all 25 k-points.
  - For the SOC-free files, which use pseudopotentials with different semicore states, it is the band below the gap continuous with the SOC gap.
- **`wire_screening`:** the electronic dielectric tensor (last block of `ph.epsil.out`), the cells of `IR/<chain>/scf.in` and `volume/<chain>/scf.in`, the number of atoms, and the Environ PV energy and pressure of the quantum-volume run.

The k-point coordinates of the QP tables, k = (n − 1)/48 in units of 2π/c, are the irreducible points of the 48-point grid. They were confirmed by matching the Eo gaps with the SOC band structures.

`make_task.py` asserts the facts the task relies on:

- the GW starting point is the SOC band structure;
- each spectrum has 2000 rows;
- the `BSE_G0W0` and `BSE_scGW` runs of Te3, As2S3 and Bi2Te3 have identical independent-particle spectra (relative difference below 1e-5), and their `BSE_G0W0` headers name the scGW QP database;
- the S3 `BSE_G0W0` header names the G0W0 database.

Phonon data and the IR spectra are not staged.

## What is hidden

**Withheld from the agent:**

- the scGW QP table of Bi2Te3 (`o-G3W3.qp`);
- the Bi2Te3 `BSE_scGW` spectrum.

**Graded values** (derived from the archive, never served as numbers):

- the scGW minimum direct gap of Te3 and Bi2Te3 with SOC (as archived) and without SOC (archived gap minus the SOC change of the minimum direct DFT gap);
- the first bright exciton (10% rule) of the `BSE_scGW` spectra of S3, Te3 and As2S3;
- Δ_W, from the first bright excitons of the two runs. For S3, the change of the minimum direct QP gap between G0W0 and scGW is subtracted, because that run used G0W0 energies;
- the Bi2Te3 binding energy: scGW gap minus the first bright exciton of the withheld `BSE_scGW` run, graded by a regime window.

## Scope

This is a retrospective analysis of released calculations, not new GW or BSE runs. The interpretation of the `BSE_G0W0` runs of Te3, As2S3 and Bi2Te3 rests on their headers and independent-particle spectra: scGW QP energies, with a kernel that differs from that of the `BSE_scGW` runs (input files `haydock.in` versus `haydock_W.in`). The archive does not document it otherwise.
