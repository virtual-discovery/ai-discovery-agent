# Source and transformations

- **Source:** Z. Liu, H. Kikuchi, Z. Wei, S. Asai, M. Enderle, U. B. Hansen, V. O. Garlea, M. D. Le, G. J. Nilsen, I. A. Zaliznyak, T. Masuda, *Observation of Switchable Chiral Magnons in an Altermagnet*, Phys. Rev. Lett. 136, 236705 (2026), arXiv:2605.14124.
- **Data:** *Switchable Chiral Magnons in an Altermagnet*, Zenodo record 20049187, DOI 10.5281/zenodo.20049187, https://zenodo.org/records/20049187.
- **License:** CC BY 4.0.
- **File used:** `Figures.zip` (5,319,363 bytes, md5 fb87bc34319008f92adee9d7dbe3ada8).
- **Prior Hamiltonian:** the exchange constants in the catalog are those of Z. Liu et al., Phys. Rev. Lett. 133, 156702 (2024). They are the Hamiltonian the authors used for the model curves of this data set.

## Construction

The scripts are in `benchmarks/construction/mnte-altermagnet-chiral-magnons/`:

- `fetch.py` downloads the archive and checks its md5.
- `calib.py` fits the withheld IN20 scans (the graded truth) and runs the reference and alternative analyses.
- `make_task.py` writes the staged, lab, private and ablation files.
- `verify_model.py` checks the spin-wave code and the HYSPEC kinematics against the archived model files.

## What the agent can buy, and what it was made from

- **HYSPEC slices:**
  - symmetric (`Figure2/Fig2c.dat`) and chiral (`Figure2/Fig2e.dat`), 44 × 63 pixels;
  - values and errors are unchanged, rounded to 6 significant digits.
- **HYSPEC cuts:**
  - g, h (symmetric) and i, j (chiral);
  - from `Figure2/Fig2g–j.dat`, unchanged.
- **HYSPEC polarized diffraction:**
  - from `Figure2/Fig2b_gray.dat` (I_y^{+0}) and `Fig2b_green.dat` (I_y^{−0}); the labels follow the legend of paper Fig. 2(b);
  - the rocking angle is re-centred on the peak, and the I^{−0} channel is interpolated onto the I^{+0} angles (the two grids differ by 0.002°).
- **IN20 polarized diffraction:** the raw ILL files `Figure3/Fig3b_red` (state A, file 043731) and `Fig3b_blue` (state B, file 043806), verbatim.
- **Catalog:**
  - **energy grid:** the energy transfer of every HYSPEC pixel, E = E_i − 2.0721 |k_i − Q|². This is derived, not archived. Here E_i = 37 meV and k_i lies 32° from a*, the paper's ω. The reduction lattice is a = 4.100 Å, c = 6.843 Å; these are the only values for which a 32° incident beam reproduces the archived trajectory `Etraj` (Fig1c.mat) exactly, to 1e−8 meV. The grid also matches the energy axis of paper Fig. 2(g) to 0.1 meV.
  - **headers of the withheld IN20 scans:** all lines above `DATA_:`.

The model files `Fig1c.mat`, `Fig2df.mat` and `Fig2ghij.mat` are not staged; they were used only for checks. `FigS2`–`FigS4` are not used.

## What is hidden

The 47 IN20 constant-energy scans at 32 meV (`Figure3/Fig3c/*` = state B, after the PND file 043806; `Figure3/Fig3d/*` = state A, after 043731). They give:

- the two peak centres in the instrument's QH coordinate: a two-Gaussian fit to the summed I_x^{+-} and I_x^{-+} of both states, binned at 0.005;
- the sign of I_x^{-+} − I_x^{+-} at each peak and for each state: Gaussian amplitudes at the fitted centres, all at 9–14 σ.

The graded exchange difference has no archived value. Its window covers the estimates from both HYSPEC cuts, the published value and the value implied by the withheld IN20 separation. The wave vector of maximal splitting is graded through the symmetry form factor of the split bond classes.

## Scope

This is a reconstruction from the authors' archived data, not new measurements.
