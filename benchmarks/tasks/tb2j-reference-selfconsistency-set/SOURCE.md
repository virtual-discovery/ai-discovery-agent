# Source and transformations

- **Source:** A. Tellez-Mora, X. He, E. Bousquet, L. Wirtz, A. H. Romero, *Systematic determination of a material's magnetic ground state from first principles*, Materials Cloud Archive record rbrcp-0tc34 (2024), DOI 10.24435/materialscloud:5m-2t. Paper: npj Comput. Mater. 10, 20 (2024), DOI 10.1038/s41524-024-01202-z, with its Supplementary Information.
- **License:** CC BY 4.0 (all files of the record).
- **Files used** (md5 from the archive API; sha256 computed on download):

| file | bytes | md5 | sha256 |
|---|---|---|---|
| `NiO.tar.gz` | 327,069 | 02b5fac1c748ee9e612d6d803efaed71 | 994cbf17a204dff1cdfdb9dd2844db2e3c94635ea2c932f49018ee360e8f8a09 |
| `FeCl2.tar.gz` | 1,312,750 | 653cd5334436fb3b69e43352704083bf | ab4f4bb5dd88558325a408e384943a8739841a82972024e23f4889a2c118716f |
| `FePS3.tar.gz` | 5,750,718 | 54ebb117ece8326050225bf6fc141d78 | dca986c02a035c4d213638fb9a72315e09c4ac54fc9dffe1c277b7f9ab2a5cdc |
| `FeP.tar.gz` | 4,700,382 | 9ca0793b00444416ef5a7ea596c741e0 | 27ca2d78df076f0a6a6d0e876bda44f2ec97528b36cbfcdfd8d394c3287405b4 |

- **Members used:**
  - the ten `exchange.out` files: `{NiO,FeCl2,FeP}/{FM,AFM}` and `FePS3/{DFT,DFT+U}/{FM,AFM}`;
  - the `README` files;
  - for validation only (not served): `FeP/AFM/dispersion_AFM`, `FePS3/DFT/FM/dispersion_FM` and `FeCl2/{FM,AFM}/dispersion_*`.

  Member sha256 values are in `benchmarks/construction/tb2j-reference-selfconsistency-set/fetch.py`.
- **Not used:**
  - the `TB2J.yaml` files, which duplicate `exchange.out`;
  - the other `dispersion_*` and `magnons_*` files, which are not served because their extrema are simple transforms of graded quantities;
  - the anisotropic (DMI, J_ani) lines;
  - CuO (used by the sibling task `tb2j-groundstate-cuo-multiferroic`), MnF2, CrSBr and Mn5Si3.

## Transformations

The scripts are in `benchmarks/construction/tb2j-reference-selfconsistency-set/`:

- `fetch.py` downloads and verifies the archives and extracts the members.
- `tb2jmodel.py` parses `exchange.out`: header lines, cell, atoms (charges and moments) and every isotropic pair (i, j, R, J_iso, vector, distance).
- `make_task.py` splits each calculation's pair list into nine windows of the archived distance (the lab data). It computes the graded targets, ablation values and calibration from the pair lists.

No values are altered: pairs are served exactly as archived (J to 1e-4 meV, vectors to 1e-3 Å), with the archived TB2J header lines, cell and atoms of each calculation.

## Validation against the paper and the authors' own output

- **FePS3 (U = 0) FM reference.** The minimum at Z = (0, 0, 1/2) lies 2.00 meV below Γ, as stated in the paper. The archived h(k) at Γ and Z is reproduced to 0.002 meV.
- **FeP FM reference.** k_min = (0, 0.0306, 0.2179) r.l.u., against the paper's (0.00, 0.03, 0.21).
- **Paper tables.**
  - FeP nearest-neighbour coupling −4.642 meV (Table 3);
  - NiO second-neighbour coupling −11.58 (FM) and −11.52 (AFM) meV, the PBE+U rows of Table 1;
  - FeCl2 J1 = 7.330 meV (SI Table 1);
  - FePS3 DFT+U J3 = −3.35 meV (Table 2).
- **FeCl2 archived h(k).** FM and AFM files reproduced at the path nodes to 0.003 and 0.010 meV.
- **FeP AFM (spiral) archived h(k).** Reproduced to 0.1 meV over the whole path from the pair list as listed and the archived moments. Its Γ–Z segment is flat to 1e-8 meV, and its minimum is −1.10 meV near Γ.

## Data properties found during construction

- **FePS3 DFT+U FM file.** It contains the same cell, atoms, moments and 1196 couplings as the U = 0 FM file; only the header (TB2J version, date, working directory) and the pair order differ.
- **R = 0 supercell lists.** The FeCl2 and FeP AFM supercell lists were computed with R = 0 along the long axis only. Bonds that cross the cell boundary are listed at the far image.
- **Non-stationary FeP moments.** The archived FeP spiral moments are not a stationary point of their own exchange (torques 0.3–4.3 meV).

## What is hidden

No graded number is served. The graded values are derived from the archived pair lists by `make_task.py`:

- the k_min of the FeP FM-reference exchange (full list);
- the k_min of the FeP spiral-reference exchange, mapped onto the four-Fe cell with minimum-image pair vectors and averaged over the five supercell copies of each bond;
- for each material, the classical stabilization energy E(FM) − E(converged configuration) per magnetic ion from the FM-reference exchange (carried onto the converged supercell by lattice translation) and from the converged calculation's own exchange, using the archived moment directions.

The FePS3 DFT+U FM-reference entry is null.

The tolerances admit the defensible alternatives listed in `solution/ablations.py` and reject the realistic wrong paths. The alternatives are shorter but converged lists, relaxed FeP moments, and other averaging of the folded bonds.

## Scope

This is a retrospective analysis of released calculations, not new DFT. Each reference calculation used its own relaxed structure: NiO's AFM lattice is about 3 % smaller, and FePS3's U = 0 AFM structure is zigzag-distorted. Each calculation's own cell, atoms and moments are served with its pairs.
