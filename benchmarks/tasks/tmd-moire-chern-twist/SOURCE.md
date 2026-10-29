# Source and transformations

- **Study:** X.-W. Zhang, C. Wang, X. Liu, Y. Fan, T. Cao, D. Xiao, *Polarization-driven band topology evolution in twisted MoTe2 and WSe2*, Nat. Commun. 15, 4223 (2024), DOI 10.1038/s41467-024-48511-x; arXiv:2311.12776.
- **Archive:** Zenodo record 10888769 (DOI 10.5281/zenodo.10888769), https://zenodo.org/records/10888769, license CC BY 4.0.
- **File used:** `source-data.tar.gz` (20,803,943 bytes, md5 038b59d589a8f254d76599e6b3e96f74). Per-member md5 sums of the 76 members read by the build are in `benchmarks/construction/tmd-moire-chern-twist/members.md5`. The `ML/` folder (including a licensed VASP POTCAR) is neither extracted nor staged.

## What is staged (unmodified archive files)

`fig1/fig1b` → `archive/fig1/set_P`, `fig1/fig1a` → `archive/fig1/set_Q` (continuum maps and plot scripts); `fig2/fig2c/3.15` and `fig2/fig2d/3.89` (`Es-band.dat`, SIESTA `.fdf`); `fig2/fig2a-b/{3.15,1.25}` (`CONTCAR`, `ILD-int.dat`, `displace-int.dat`); `fig3/fig3a-d/3.15` (`pot_2d_top_sm.dat`, `pos_2d_bot_sm.dat`, `plot-pot.input`, `plot-pot.py`), `fig3/fig3e/3.15/pot_diag-top.dat`, and the scripts `plot-diag-pot.py`, `plot-fft-pot.py`.

## What the lab serves

- `dft_bands`: `Es-band.dat` of the eight `fig2c`/`fig2d` calculations, converted to meV, states sorted from the highest.
- `c3`: the C3z eigenvalue labels of the paper's Fig. 2c/2d (they exist only in the figure). They were transcribed by hand from a 450-dpi rendering of page 2 of arXiv:2311.12776v3 into `c3_labels.json` of the construction folder. Each label gives the spin-up (K) eigenvalue and its spin-down (K′) partner (complex conjugate) for the top three doublets at γ and κ; the caption states that κ′ equals κ. The Chern numbers implied by these labels agree with the C_K↑ values printed in the same figure for every panel. Doublet energies come from `Es-band.dat` (indices 0 and 40 of the path). Within each doublet, the order of the spin-up and spin-down entries is fixed but arbitrary.
- `hartree_map`, `hartree_linecut`, `hartree_fourier`, `structure`: `fig3/fig3a-d/<angle>/pot_2d_{top,bot}_sm.dat` (reshaped to n × n, in eV), `fig3/fig3e/<angle>/pot_diag-top.dat`, `fig3/fig3f/<angle>/pot_top_bot_G.dat` (rows with |g|/g1 ≤ 7) and `fig3/fig3a-d/<angle>/POSCAR`.

## What is hidden and how the targets were obtained

- **DFT Chern pairs** of all eight calculations: from the transcribed Fig. 2 labels with the indicator exp(2πiC/3) = −ξ_γ ξ_κ ξ_κ′ (equal to the C_K↑ values printed in Fig. 2).
- **Continuum Chern pairs** of both parameter sets at 3.89° and of set P at 2.14°, plus the set matching the DFT bands: from our own numpy plane-wave evaluation of the archived maps (`continuum.py`, `reference.py`), converged in the plane-wave cutoff (61–169 plane waves per layer) and the Fukui–Hatsugai–Suzuki mesh (12² to 30²). These calculations reproduce the paper's statements: the local-stacking set gives C1 = −1 and the DFT set gives C1 = +1 for K↑. Set P also reproduces the 3.89° SIESTA bands, with an RMS distance of 1.6 meV against 7.6 meV for set Q.
- **Hartree values:** Δv_H = top − bottom of the archived maps at XM (fractional (1/3, 1/3) of the archived cells, where top-layer Se sits above bottom-layer W in all four relaxed structures), evaluated by Fourier interpolation, and the first-star part at XM.

## Scope

Retrospective analysis of released data. The Chern numbers are the authors' symmetry-indicator results; the continuum values are validated-evaluator values for the archived maps, with m* and a₀ taken from the papers that introduced each parameter set (local stacking: 0.62 mₑ, 3.472 Å; DFT fit: 0.60 mₑ, 3.52 Å).
