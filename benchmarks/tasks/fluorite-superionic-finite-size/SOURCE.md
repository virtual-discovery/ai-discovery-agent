# Source and transformations

- **Source:** F. Grasselli, *Investigating finite-size effects in molecular dynamics simulations of ion diffusion, heat transport, and thermal motion in superionic materials*, J. Chem. Phys. 156, 134705 (2022), doi:10.1063/5.0087382. Data: Materials Cloud Archive record yzhmw-vje16, DOI 10.24435/materialscloud:jy-tw, https://archive.materialscloud.org/records/yzhmw-vje16.
- **License:** CC BY 4.0 (all files).
- **Archive files read** (by HTTP range requests into the remote zips; the archives were never downloaded whole):
  - `CaF2.zip` (804,978,718 bytes, md5 7a88bad620f87627aa4acb18266254d1), folder `CaF2/NVT_L_T/`;
  - `PbF2.zip` (1,205,555,455 bytes, md5 7c99d63a2e0e0662a2aa634496405e39), folders `PbF2/NVT_L_T_6p056/` and `PbF2/NVT_L_T_6p056_DW/`;
  - `UO2.zip` (317,257,266 bytes, md5 ba43b29a5c7308bcca86b29041c559a6), folder `UO2/NVT_L_T/`;
  - `alpha-AgI.zip` (600,676,721 bytes, md5 6532ccd122eaf453975676e0e8acbbc0), folder `alpha-AgI/NVT_L_T_5.37ang/`.
- **Members used:** for every run, `msd_CM_B4.dat`, `msd_self_B4.dat`, `nvt_log` and the input deck (`*.in`), plus each folder's run scripts: 1207 members in all. Each member's SHA256 is listed in `benchmarks/construction/fluorite-superionic-finite-size/members.json` and checked against the zip CRC32 on refetch.
- **Not used:** none of the `*_heat` folders and not `UO2_heat.zip`. Those hold the heat-flux logs used by the sibling task `superionic-heat-transport-size-onsager`.

## Construction

The scripts are in `benchmarks/construction/fluorite-superionic-finite-size/`:

- `fetch.py` reads the members through `remotezip.py`.
- `extract.py` builds one record per run: 286 runs (CaF2 54, PbF2 71, UO2 111, AgI 50). `PbF2/NVT_L_T_6p056_DW` repeats the main PbF2 set byte for byte at n = 2–4. The script checks this, and only its n = 1 and n = 5–7 runs (and n = 2–4 at 400 K) are added.
- `make_task.py` writes the catalog, the lab data, the private targets and the ablation values.

What a purchase returns:

- the run's input deck, verbatim;
- a subset of the rows of `msd_CM_B4.dat`: 140 lags (220 for UO2), all columns except the type-1 and centre-of-mass variances;
- a thermodynamic summary of `nvt_log`: means, 20-block standard errors, the total-energy variance and the production length.

No values are modified.

## What is hidden

- **Withheld runs:** all 6×6×6 and 7×7×7 PbF2 runs, all 6×6×6 AgI runs, and the 768-atom runs CaF2 1250 K and UO2 2450 and 3150 K.
- **Graded targets:**
  - D of seven withheld runs, from the released MSD (lab-frame fit over 20–100 ps);
  - T_c of the 768-atom cells (central-difference c_V of the released mean energies, parabolic peak);
  - CaF2 and PbF2 activation energies of the 768-atom cells (Arrhenius fits to diffusive runs at least 100 K below or above T_c).

The protocol is in `analysis.py`.

## Scope

This is a retrospective analysis of released empirical-force-field MD (Fumi–Tosi rigid-ion CaF2 and PbF2, CRG-type EAM UO2, Coulomb–Morse AgI), not new simulations. The graded transition temperatures and activation energies are those of these force fields and cells, not experimental values.
