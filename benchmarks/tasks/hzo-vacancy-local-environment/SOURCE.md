# Source and transformations

- **Paper:** S. Taper, B. Cai, N. Strkalj, B. Monserrat, G. Di Martino, *Role of local cation environment on the formation energy of oxygen vacancies in ferroelectric Hf1−xZrxO2*, Phys. Rev. B 113, 134114 (2026), DOI 10.1103/vz64-n5xh.
- **Data:** the same authors, *Research Data supporting "Role of Local Cation Environment on the Formation Energy of Oxygen Vacancies in Ferroelectric Hf1−xZrxO2"*, Apollo – University of Cambridge Repository (2026), DOI 10.17863/CAM.121316, https://www.repository.cam.ac.uk/handle/1810/389433.
- **License:** CC BY 4.0, per DataCite (rightsIdentifier cc-by-4.0) and the repository record. Every staged file derives from this record and carries the same license.
- **Files fetched:**
  - `Data.zip`: 12,791,098 bytes, md5 70721b2c03a34c84e9a2deab61160afb, sha256 e6fd7f7aa4ed5ba2154fcf06dd9808d8c643637576b8b489887780942a0c4ea9.
  - `Readme.md`: 15,116 bytes, md5 77a9d035a0ac4ff76779ff70356b7934.

## Members used (sha256, bytes)

| Member (under `Apollo deposit/`) | sha256 | bytes |
|---|---|---|
| `Main text/Structures/Fig2b_e_0.25_HZO.cif` | 81500cb808458e2bf8f6e8559f56b5daacc4a63a675146fc2767d80c080d26e8 | 5922 |
| `Main text/Structures/Fig2c_f_0.50_HZO.cif` | 402e18c0063fa19ad96a135bdfd61ab213fdec3ad8bb754e9f1c0a84db68bd85 | 5909 |
| `Main text/Structures/FIG2d_g_0.75_HZO.cif` | f74eec66cdde98e487c0212f921cba8d48954ace3179e670df71a16b82f6c2bc | 5890 |
| `Main text/Fig_2/Data_Fig2b_g.csv` | 642a93a8a57605d4a8fbdb40e9a2bed325a7079e5020bad897770c8285b95665 | 47976597 |
| `SM material/FIGS3.csv` | 848557f933a2ba6bd596797726f8b7d52418c7774c65940798b65955916b7f58 | 27604 |
| `SM material/FIGS4.txt` | ece29729a7c957ec312d18707babe9e467ca2bb40ee3ba85349cfe6afedea945 | 972 |
| `SM material/FIGS4 3x3x3_supercell.cif` | abd21fb28670226729bfa05192db7727cf5b1d1bfe03e1efbd0910981887cb93 | 18513 |
| `SM material/FIGS5.txt` | ab7b7ccfb5fddea247394f0a681d864134618a9d9a826d72c33038dea651fca7 | 479 |
| `SM material/FIGS12.csv` | 0dac2a75b8ad7929199aef927300b245a8ac47e9133d7c38feb5720ead1ae0ef | 6944 |

## Construction

The scripts are in `benchmarks/construction/hzo-vacancy-local-environment/`.

- **`fetch.py`** downloads the two files from the repository's bitstream API, verifies their checksums and extracts the members above.
- **`extract.py`** writes `pool.json` and asserts every relation used later against the archive:
  - From `Data_Fig2b_g.csv` it reads the neutral-vacancy and pristine total energies of the three supercells: 64 O sites each, VASP-PBE. These columns do not depend on the Fermi level. The released neutral formation energy equals E_vac − E_bulk + μ_O for μ_O = −10.09, −7.33 and −4.57 eV to 1e-9 eV.
  - From `FIGS12.csv` (x = 0.50 cell) it reads the site table. Its `E_vac` equals the O-rich (μ_O = −4.57 eV) neutral formation energy of the `Fig2c_f_0.50_HZO.cif` cell to 1e-9 eV, although the release does not state this.
  - It also reads the CASTEP enthalpies of `FIGS3.csv` (seven functionals × 64 sites, with the authors' labels) and `FIGS4.txt` (45 CASTEP-PBE sites of the 3×3×3 cell).
- **`make_task.py`** writes the agent catalog, the lab lookup data, the private targets and the ablation and robustness values. It uses `model.py`.

Staged for the agent (`environment/agent/`):

- the three composition CIFs and the 3×3×3 CIF, unchanged apart from line endings normalised from CRLF to LF and the file names;
- `x050_site_table.csv`, which is `FIGS12.csv` unchanged in content, rewritten with LF line endings;
- `catalog.json`: site lists, the CASTEP site coordinates as released, and the conventions.

The authors' CASTEP labels and all energies are served only by the lab, one purchase at a time.

The authors' README says that `FIGS3.csv` belongs to `Fig2a_Fig3_HZO.cif`. Its per-functional labels and coordinates, however, show that each functional's calculation used its own Hf/Zr decoration. For every functional, at most 24 of the 64 released labels match the environment of the same site index in any of the three composition CIFs (7–24 matches, against about 64 expected for a shared decoration). The released labels fit the energies with residuals of 0.02–0.04 eV in a sublattice-plus-Zr-count model. The environments of the x = 0.50 VASP cell, assigned by site index, give 0.10–0.16 eV. The labels are therefore served as released, and `Fig2a_Fig3_HZO.cif` is not staged.

## What is hidden

The lab never serves the following; they grade the submission:

- **Site energies:** the O-poor neutral formation energies of 20 withheld sites: 10 in the x = 0.25 cell, including its only four-Zr site, and 10 in the x = 0.75 cell, including both of its Hf-rich fourfold (one-Zr) sites.
- **Populations:** the equilibrium populations at 700 K in a macroscopic random solid solution at x = 0.25 and 0.75 (distribution over the number of Zr first-shell neighbours, and the site fraction). They are computed from all 64 archived site energies of each composition: for every environment class, the Boltzmann factor is averaged over the class's archived sites and weighted by binomial probabilities. Classes absent from a cell take that cell's per-sublattice additive fit; they carry at most 0.2 % of the population.
- **LDA description:** the least-squares common-slope description of all 64 CASTEP-LDA vacancies, using the authors' labels.

## Separation from CompMat-Bench

CompMat-Bench (github.com/cz2014/CompMat-Bench, tasks 202603_025–030) uses this study and publishes its reference answers. Those targets are:

- environment class counts, relaxation displacements and the shell gap of the `Fig2a_Fig3` (reference x = 0.50) cell;
- the dielectric constant (ε_iso = 22.816) and the FNV-corrected formation energies of nine representative sites of that cell in five charge states, with transition levels;
- COHP integrals.

They are CompMat's own reproductions at x = 0.50. This task grades none of these quantities and stages nothing from CompMat-Bench. It uses neither the `Fig2a_Fig3` cell nor its nine-site, charged-defect, HSE or COHP data. Its graded values are:

- archived per-site energies of the x = 0.25 and x = 0.75 cells;
- populations derived from them;
- a fit to the CASTEP-LDA set.

None of these appears among CompMat-Bench's targets, and the paper reports them only as figure-level trends.

## Scope

This is a retrospective analysis of released DFT results with the authors' protocol, not new calculations. All energies are for neutral vacancies. The three VASP supercells share one lattice (10.46 × 10.00 × 10.10 Å).
