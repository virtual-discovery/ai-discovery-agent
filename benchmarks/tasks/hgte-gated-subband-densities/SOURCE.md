# Source and transformations

- **Study:** M. Hofer et al., *Self-consistent electrostatic modeling of gated narrow-gap topological insulators*, Nano Lett. 26, 3042–3048 (2026), DOI 10.1021/acs.nanolett.5c05762.
- **Calculations and processed data:** Zenodo record 18709991 (DOI 10.5281/zenodo.18709991), `data_repository.zip`, 498,134,382 bytes, md5 a948d21c70f3f911d3bbc8554acccbce. License CC BY 4.0.
- **Raw transport:** Zenodo record 18429111 (DOI 10.5281/zenodo.18429111; concept 10.5281/zenodo.15753790), `data_and_code.zip`, 57,037,973 bytes, md5 e777270ff41037b2072b819a8d4b9961. License CC BY 4.0.
- **Code:** kdotpy (GPL-3.0) produced the archived runs. No kdotpy code is shipped or run.

## Files used

Single members were read from the remote zips with HTTP range requests. The md5 of all 397 members fetched is in `benchmarks/construction/hgte-gated-subband-densities/members.md5`. The files that feed the task are:

- `measurements/{26nm,70nm,110nm}.csv` (raw lock-in records; md5 a0e0a8a0…, a9d458cd…, c47252a7…)
- `FigureS3/calculations/HgTe{26,70,107}nm_FFT_*/…` (`.kpout`, `potential*.csv`, `dos-integrated-byband*.csv` of 80 runs)
- `FigureS3/README.md` (md5 dfc374ad…, staged verbatim in the catalog)
- `FFT_experiment/70nm_FT.csv` (md5 53c7544f…, private: the source of the experimental targets)
- for checking only: `FigureS3/FigureS3.ipynb`, `FFT_experiment/FFT_experiment.ipynb`, `Fig7/dual_carrier_Drude_fit.ipynb`

## Transformations

- **Transport.** The raw records were restricted to B ≤ 2.0 T, with no other change. `hall_scan` is the least-squares slope of Vxy/VI against B over |B| ≤ 0.2 T, plus Vxx/VI at the field closest to 0. The 107 nm sample's record is the archive's `110nm.csv` (nominal 110 nm, modelled by the authors as 107 nm).
- **Calculations.** The integrated DOS per band (bands −3…6 and the total) was cut to −60…80 meV on the archived 0.1 meV grid. The `.kpout` text and the potential profile are served unchanged. All values are rounded to 7 significant digits.
- **Gate calibration.** The convention stated in the instruction reproduces the authors' 70 nm line (9.8197 Vg + 3.2060) exactly.
- **Experimental targets.** The authors' 70 nm SdH chart (`70nm_FT.csv`) is reproduced to 1e-14 from the raw record by `sdh.py`. Each target is the parabolically refined maximum of the chart's log power averaged over Vg ± 0.02 V, on a ridge assigned by tracking it across the 2D chart. The hole-onset voltage is the zero of a linear fit to the majority-pocket ridge at −0.75…−0.55 V.
- **Model targets.** These are the archived runs' band IDOS at E_F(n) of the main run, interpolated linearly in n at n_tot(Vg). Hole pockets use the wedge runs' band −1 at the main run's E_F, divided by 360°/width, following the authors' `FigureS3.ipynb`.

## What is hidden

The agent never sees:

- the authors' processed SdH charts (all samples) and their figure notebooks;
- the 45 nm sample (not used at all).

The graded values are:

- six 70 nm SdH densities and the hole-onset voltage, from `70nm_FT.csv`;
- three model densities from the 70 nm runs;
- the majority-pocket multiplicity, 2.

## Separation from CompMat-Bench

This study is CompMat-Bench source study #5 (tasks 202602_001–004), whose answers are public. Every CompMat-Bench task uses the 45 nm quantum well:

- Hartree inputs and analysis;
- surface dispersion;
- the SdH peak densities at Vg = 0, −0.2 and −0.8 V and the symmetric point of `45nm.csv`.

This task neither stages nor sells any 45 nm file. It grades only 70 nm quantities, plus model values from the 70 nm runs, so no graded value coincides with a CompMat-Bench target. The 70 nm curves appear only graphically in the paper's Fig. S3, with no tabulated numbers.

## Scope

This is a retrospective analysis of released measurements and calculations; no new kdotpy runs were made. The archive holds one band-structure parameter set (E_P = 21 eV, VBO = −620 meV, ε = 20.8, sample-specific boundary potential), so dielectric-constant and VBO calibration were not offered.
