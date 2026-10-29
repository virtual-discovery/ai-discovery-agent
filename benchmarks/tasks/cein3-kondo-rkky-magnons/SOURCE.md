# Source and transformations

- **Source:** W. Simeth, Z. Wang, E. A. Ghioldi, D. M. Fobes, A. Podlesnyak, N. H. Sung, E. D. Bauer, J. Lass, S. Flury, J. Vonka, D. G. Mazzone, C. Niedermayer, Y. Nomura, R. Arita, C. D. Batista, F. Ronning and M. Janoschek, *A microscopic Kondo lattice model for the heavy fermion antiferromagnet CeIn3*, Nat. Commun. 14, 8239 (2023), DOI 10.1038/s41467-023-43947-z.
- **Dataset:** "Dataset package for the manuscript ...", Zenodo record 10146787, DOI 10.5281/zenodo.10146787.
- **License:** CC BY 4.0 (all files).
- **Files:** all 48 `.dat` files of the record (1,005,608 bytes in total). Names, sizes and md5 checksums are listed in `benchmarks/construction/cein3-kondo-rkky-magnons/members.json`.

## Files served by the lab (values unchanged, NaN written as null)

- **12-meV data:**
  - maps `HighEnergy_Int_RG.dat`, `HighEnergy_Int_GX.dat`, `HighEnergy_IntXM.dat` and `HighEnergy_Int_MG.dat`, with their `_Q_` and `_E_` axis files;
  - sold one constant-Q row at a time;
  - the 9 R–Γ rows with no measured pixel are not listed.
- **High-resolution data:**
  - maps `HighResolution_Int_RG.dat` and `HighResolution_Int_RX.dat`, with their axis files, sold one constant-q row at a time;
  - constant-energy cuts `HighResolution_{RG,RX}_ConstEcut_{0p6,1p4}meV.dat`.
- **Band structure:** `fbands.dat` and `allotherbands.dat`, rounded to 0.1 meV or better.

## Withheld and graded

- **High-resolution R–M data:** the map `HighResolution_Int_RM.dat` (with its axes) and the cut `HighResolution_RM_ConstEcut_1p4meV.dat`.
- **Q-integrated χ″:** `Experiment_IntegralOverchi.dat`.
- **Grading of predictions:** each prediction is compared with these archived values by a reduced χ² after the best non-negative scale (plus a constant offset for the R–M map and cut).
  - The cuts and the integrated χ″ use their archived error bars.
  - The map has no archived errors. Its noise (0.28 in the archived units) is estimated from the difference of the +q and −q rows, at energies above 0.25 meV.
- **Scalar windows:** W, η and J1 are graded against windows calibrated from fits to the released data (`calibrate.py`).

## Withheld and not used

- `HighResolution_RM_ConstEcut_0p4meV.dat`: it does not discriminate between the candidate models.
- `Experiment_ConstQ1.dat` and `Experiment_ConstQ2.dat`: their Q positions are not archived.
- The theory files: `Theory_*`, `MOPAM-Dispersion_*`, `J1Model-Dispersion_*` and `InteractionSpectrum.dat`.
- The duplicate `HighResolution_IntRG.dat`.

## Metadata taken from the paper (Methods) and stated in the catalog

- CNCS instrument, T = 1.8 K, a = 4.689 Å.
- Incident energies 12 meV and 3.315 meV; elastic resolution 0.106 meV FWHM for the high-resolution setting.
- The reduction steps behind χ″.
- The momentum integration windows:
  - ±0.17 r.l.u. perpendicular to Γ–X and X–M;
  - ±0.09 r.l.u. perpendicular to R–Γ and M–Γ;
  - "within a distance of 0.08 r.l.u." perpendicular to the high-resolution lines, quoted verbatim;
  - ±0.1 meV energy integration for the cuts.
- The path coordinate convention (R = 0, Γ = 0.866, X = 1.366, M = 1.866) follows from the archived axis files.

## Scope

This is a retrospective analysis of released data. No graded quantity is a literature number. The paper's headline values (W = 2.75 meV, v = 38.6 meV/r.l.u., η = 2.23) lie inside the graded W and η windows, but the three predictions require the data.
