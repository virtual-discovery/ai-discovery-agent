# Source and transformations

- **Source:** F. Grasselli, *Investigating finite-size effects in molecular dynamics simulations of ion diffusion, heat transport, and thermal motion in superionic materials*, J. Chem. Phys. 156, 134705 (2022), doi:10.1063/5.0087382. Data: Materials Cloud Archive record yzhmw-vje16, doi:10.24435/materialscloud:jy-tw (same record as the sibling task fluorite-superionic-finite-size).
- **License:** CC BY 4.0 (all files).
- **Archives read (by HTTP range requests into the remote zips, no full download):**
  - `UO2_heat.zip` (md5 de844afbecd3fd5e72812ccb341644ad): `UO2_heat/NVE_L_T/`, `UO2_heat/NVT_L_T_heat/`
  - `CaF2.zip` (md5 7a88bad620f87627aa4acb18266254d1): `CaF2/NVE_L_T_heat/`, `CaF2/NVT_L_T_heat/`
  - `PbF2.zip` (md5 7c99d63a2e0e0662a2aa634496405e39): `PbF2/NVE_L_T_6p056_heat/`, `PbF2/NVT_L_T_6p056_heat/`
  - Per run: the thermo log (`nve_log`/`nvt_log`) and the input deck. SHA256 of all 511 members: `benchmarks/construction/superionic-heat-transport-size-onsager/members.json`.
- The NVT diffusion folders used by the sibling task are not used here.

## Transformations

Scripts in `benchmarks/construction/superionic-heat-transport-size-onsager/`:

1. `fetch.py` reads the members, checks CRC32 and SHA256, and keeps the log columns Time, TotEng, PotEng, Temp, Press, Lx, `c_myHeatFlux[1-3]` and the anion velocity sums (float32); the six stress columns are dropped.
2. `extract.py` takes the first 800 ps (100000 frames at 8 fs) of each of the 248 runs, band-limits the six flux columns to 6.25 THz by FFT truncation, resamples them to 0.08 ps (10000 frames) and stores each column as int16 with a per-column scale (quantisation noise below 1e-8 of the variance). It records cell edge, mean and standard deviation of the temperature, mean pressure and total-energy drift, and the deck.
3. `analysis.py` computes κ of every run (multicomponent cepstral reference and variants); `cepstral.py` holds the estimator.
4. `make_task.py` writes the catalog, `data.json` and `series.npz` (served runs only), the private targets and the ablation values.

## What is hidden

- All released runs of four states, which are not in the lab: UO₂ 4×4×4 at 1700 K and 3400 K (NVE and NVT), CaF₂ 4×4×4 at 1600 K (NVE and NVT), PbF₂ 4×4×4 at 1200 K (NVT; no NVE run exists). Target: mean of the multicomponent cepstral κ of the withheld runs of each state.
- κ of three served runs (`report_runs`), computed with the same reference estimator.
- The thermostat verdict (NVT/NVE κ ratio over 116 released pairs: 0.998 ± 0.064).

## Scope

Retrospective analysis of released classical MD (empirical force fields at fixed volume); the values are protocol-specific and not experimental. The reference κ is computed by the builder from the released flux series with the multicomponent cepstral method of the thermocepstrum code included in the record (re-implemented in numpy); the paper's own κ values are not used.
