# Source and transformations

- **Source:** H. Sugawara, *A dipole-forbidden conduction-band transition and the mid-infrared absorption of n-type Mg2Si, Mg2Ge and Mg2Sn: all-electron quasiparticle self-consistent GW calculations* (to be submitted to Phys. Rev. B, 2026). Not yet peer reviewed.
- **Archive:** Zenodo record 23319151, DOI 10.5281/zenodo.23319151 (published 2026-10-13), https://zenodo.org/records/23319151.
- **Licenses** (from the archive's `LICENSES.md`):
  - computed outputs, the JSON written by the programs, and the documentation: CC BY 4.0;
  - digitized experimental points (`data/digitized/`): CC0 1.0;
  - programs (`*.jl`) and the Quantum ESPRESSO/BerkeleyGW inputs: MIT;
  - ecalj and Questaal input files: AGPL-3.0-or-later. None is staged or served.
- **Files used:**
  - `20261012aa-paper1-zenodo-archive.zip` (77,333,006 bytes, md5 865ba37d739a8118a8a556cb95c643bb);
  - `README.md` (md5 246ab9331d014b4665cb1629ca218586), `LICENSES.md` (md5 aed582f0574ab7043a3627ea9ac59b82) and `MD5SUMS.txt` (md5 91946a35a975e8bfcda1d31e4a6db173).

  The 776 members read are listed with their md5 (equal to the archive's own `MD5SUMS.txt`) in `benchmarks/construction/mg2x-midir-forbidden-transition/members.md5`:
  - `data/bands/mg2{si,ge,sn}-{lda,qsgw}/bnd00[1-4].spin1`;
  - `data/mass-scan/mg2{si,ge,sn}-{lda,qsgw}/bnd00[12].spin1`;
  - `data/qbox/{qsgw-qbox2,qsgw-qbox-step0.01,qsgw-qbox-step0.005,mg2ge-qbox,mg2sn-qbox}-scf8-cont/bnd*.spin1`;
  - `data/pmat-pbe/box_pmat.csv`;
  - the ecalj TOTE files `ecalj-supplementary/mg2{si,ge,sn}-g0w0-oneshot/TOTE.UP`, `mg2{si,ge}-g0w0-oneshot-pbe/TOTE.UP`, `qsgw-iterations/mg2{si,ge}-666-*-scf8-cont/TOTE.UP.it10`, `mg2{si,ge,sn}-666-*-scf8-cont-ss080/TOTE.UP` and `mg2si-666-a6322-scf8-cont/TOTE.UP`;
  - the Questaal eigenvalue extracts `questaal/{prod-k4-es,basis2-k4-es,bsw-nc8,bsw-nc12,k4-es-lda}/*.eval.txt`;
  - the BerkeleyGW `data/berkeleygw/*/eqp1.dat`;
  - `eph/analysis/eph_renormalization.json`;
  - the digitized curves `data/digitized/{koenig1961-band-curves-th140,20260914-DATA-lott1966-mg2ge-curves,20260916-DATA-koenig1961-edge,20260914-DATA-koenig1961-areas}.json`;
  - the archive's analysis outputs, used only for checks: `analysis/*.json`, `data/mass-scan/effective-masses.json`, `data/quasiparticle/*.json`, `ecalj-supplementary/*.json` and `questaal/questaal_levels.json`.

## Construction

The scripts are in `benchmarks/construction/mg2x-midir-forbidden-transition/`:

- `fetch.py` downloads and verifies the files and members.
- `model.py` holds the readers and a Python port of the archive's Julia estimators (`qbox.jl`, `ir_raman_ratio.jl`, `effective_masses.jl`; MIT, H. Sugawara).
- `make_task.py` writes the catalog, the staged spectra, the lab data, the private targets and the ablation values, and asserts the archive facts the task rests on.

The lab serves:

- **Band files** (`bands`, `x_scan`, `qsgw_box`): the energies of the listed bands (eV, relative to the valence-band maximum, as written by ecalj), rounded to 1e-5 eV (paths) or 1e-6 eV (scans and boxes). Wave vectors are read from columns 5–7 of the band files. For scans and boxes they are converted to the displacement from X.
- **`momentum_matrix`:** the columns `E13_eV`, `E14_eV` and `P2_13_14` of `box_pmat.csv`, unchanged.
- **`tote`:** the rows of the TOTE file at Γ and X for bands N_v − 1 to N_v + 3, as text exactly as written by ecalj, with the file's first line and the column description of the archive's README, quoted verbatim.
- **`questaal`:** the eigenvalues (Ry) of bands 1–12 at Γ and X. The header lines of the extract are dropped, since they quote the derived splitting.
- **`pseudopotential_g0w0`:** the eqp1.dat rows at Γ and X, unchanged.
- **`electron_phonon`:** the archive JSON entries. For AHC, the per-state renormalizations at each broadening; the derived `dE` and `gap` entries are dropped. For ZG, the archive's changes of X3 − X1 and of the gap.
- **`experiment.json`** (staged, CC0): the digitized curves only. Photon energies are converted to meV with 0.1 meV resolution and the ordinates rounded to 4 significant digits. The author's carrier concentrations of the Koenig sample are included. The archive's peak readings and area ratios are omitted.

`make_task.py` asserts:

- that the scans, path and box files agree at X;
- that the X3 − X1 splittings reproduce the archive's `quasiparticle-path.json` and `tote_levels.json`;
- that in every one-shot G0W0 TOTE row column 7 equals ε + Z(column 8 − ε), so that column 7 is the linearized, Z-renormalized solution (the README's label for it reads "eQP (Z = 1)");
- that the QSGW final-iteration TOTE (column 6) reproduces the band files;
- that the longitudinal masses identify X1 and X3 (X3 heavier by a factor above 1.7 in every compound);
- that the Python port reproduces the archive's `ir_raman_ratio.json` (1e-9) and `effective-masses.json`.

## What is hidden

The graded values, derived from the archive and never served as numbers:

- **X3 − X1:** the signed splitting at X (LDA, one-shot G0W0 from the LDA, QSGW) for the three compounds.
- **Ladder correction:** Questaal QSGW + ladder minus QSGW for Mg2Si, production basis and window.
- **Infrared/Raman ratio Q:** ⟨q⊥²⟩(300 K)/⟨q⊥²⟩(77 K) at fixed n with Fermi–Dirac occupation over the QSGW boxes. For Mg2Si these are the archive's values (`analysis/ir_raman_ratio.json`, `points_box`, n = 1e16 and 3e18 cm⁻³); for Mg2Sn (n = 1e16) it is the same estimator applied to the archive's coarse box.
- **Predicted splittings:**
  - Mg2Si at T → 0: the QSGW splitting at the 100 K lattice (6.3218 Å) plus the AHC zero-point change of X3 − X1;
  - Mg2Si and Mg2Ge at 300 K: the production QSGW splitting plus the AHC change at 300 K.

  The AHC values are for the 24³ grid, non-adiabatic, η = 10 meV, Ge 3d in the core.

## Scope

This is a retrospective analysis of released calculations; nothing is recomputed. The archive's paper is not yet peer reviewed. The Mg2Sn ratio and the predicted splittings combine archived results under stated conventions; they are not numbers printed in the archive.
