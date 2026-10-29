# Source and transformations

- **Source:** A. Kundu, Y. Chen, X. Yang, F. Meng, J. Carrete, M. Kabir, G. K. H. Madsen, W. Li, *Electron-induced nonmonotonic pressure dependence of the lattice thermal conductivity of θ-TaN*, Phys. Rev. Lett. 132, 116301 (2024); arXiv:2307.15407.
- **Archive:** Zenodo record 10564835, DOI 10.5281/zenodo.10564835, "Data for Electron-induced non-monotonic pressure dependence of the lattice thermal conductivity of θ-TaN". License: CC BY 4.0.
- **File used:** `IFCs-data-pressure-TaN.tar.gz` (25,502,060 bytes, md5 0e91973a88688ca03ac776d379486f9c). Per-member md5 checksums are in `benchmarks/construction/tan-pressure-kappa-channels/members.json`.

## Contents used

For P = 0, 20, 40, 60, 80, 100 GPa: `POSCAR`, `CONTROL` (ShengBTE/FourPhonon input: 16×16×16 grid, scalebroad 0.1, isotopes, four_phonon), `FORCE_CONSTANTS_2ND` (phonopy format, 5×5×5 supercell, eV/Å²), `FORCE_CONSTANTS_3RD` (ShengBTE format, eV/Å³), `FORCE_CONSTANTS_4TH` (FourPhonon format, eV/Å⁴), `BTE.sr_phel_metal` (phonon–electron rates), and `README.md`.

## Transformations

The construction scripts are in `benchmarks/construction/tan-pressure-kappa-channels/`:

- `tanbte.py` is a numpy implementation of ShengBTE/FourPhonon phonon transport for this two-atom cell: phonopy-format harmonic constants with minimum-image weights, three-phonon rates with adaptive Gaussian broadening, Tamura isotope scattering, four-phonon rates, RTA and iterative BTE. It is served to the lab unchanged. It is validated against the archive by reproducing the ShengBTE irreducible q-point list (270 points) and by the frequency column of `BTE.sr_phel_metal`, which it matches to about 1% (median), the residual being a different phonon calculation in the phonon–electron workflow.
- `fourph.py` computes four-phonon rates on the 16×16×16 grid from `FORCE_CONSTANTS_4TH` by unbiased Monte Carlo sampling of (q′, q″) pairs (8192 per irreducible q, every branch combination evaluated exactly; FourPhonon's own sampling acceleration follows the same idea). Two independent seeds are averaged at 0 and 80 GPa; the seed-to-seed difference of the graded conductivities is reported in the construction README.
- `make_task.py` stages, for every pressure, the harmonic and third-order constants (compact npz), the POSCAR/CONTROL/README files verbatim, the phonon–electron tables at 0, 20, 40, 60 and 100 GPa verbatim, and four-phonon rate lookups at (GPa, K) = (0, 300), (20, 300), (60, 300), (100, 300), (20, 600), (100, 600). `reference.py` computes the references with the lab's own request handler.

## What is hidden

- `FORCE_CONSTANTS_4TH` at 40 and 80 GPa (and any four-phonon rates derived from them), the four-phonon rates at 60 GPa and 600 K, and `BTE.sr_phel_metal` at 80 GPa.
- Graded references (all channels, archived protocol, iterative three-phonon solution with the other channels diagonal): κ_xx and κ_zz at 80 GPa/300 K, κ_xx at 40 GPa/300 K, κ_xx at 60 GPa/600 K, and the phonon–electron share of the 80 GPa resistivity. They are new calculations from the archived inputs; the paper's own κ(P) curve is not used or graded.

## Scope

Retrospective reconstruction from the authors' archived inputs with our own implementation of their protocol; the authors' FourPhonon outputs are not in the archive, so absolute values may differ from the published figures by a few percent.
