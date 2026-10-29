# Source and transformations

- **Study:** K. P. Kluczyk, K. Gas, M. J. Grzybowski et al., *Coexistence of anomalous Hall effect and weak net magnetization in collinear antiferromagnet MnTe*, Phys. Rev. B 110, 155201 (2024), DOI 10.1103/PhysRevB.110.155201; arXiv:2310.09134.
- **Archive:** Zenodo record 18863304, *Data for: Coexistence of anomalous Hall effect and weak magnetization in a nominally collinear antiferromagnet MnTe*, DOI 10.5281/zenodo.18863304, https://zenodo.org/records/18863304.
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `dft_data.zip`: 4,006,544 bytes, md5 cef0196e73cec5deefad7e707137b940;
  - `fig6_script.zip`: 93,091 bytes, md5 15e2f717c3d2c28f11c19025122769ce.
- **Members:** all 34 members are listed with their md5 in `benchmarks/construction/mnte-altermagnet-ahc-doping/members.md5`, for example `dft_data/2_wannier90/wannier90_hr.dat` (md5 db6340b8…), `dft_data/4_wtools/sigma_ahc.txt` (md5 4648d4f1…) and `dft_data/3_wannsymm/symmetries.dat` (md5 79966d13…). `fig.zip` (figures; its `figS10.dat` repeats the withheld Hall curve on a wider window) is not used.

## Construction

The scripts are in `benchmarks/construction/mnte-altermagnet-ahc-doping/`:

- `fetch.py` downloads and verifies the two zip files and their members.
- `make_model.py` converts `wannier90_hr.dat` to npz (exact; the lab copy stores H(R) as integer micro-eV, which is exact for the 6-decimal file).
- `make_cache.py` evaluates the archived Hamiltonian with the lab's own evaluator (`mnte_tb.py`) on Γ-centred N × N × N/2 meshes (N = 48 … 400) and stores Berry-curvature-weighted and state-count energy histograms (0.1 meV bins, 4.80–8.80 eV). The lab's `ahc` and `electron_count` requests evaluate these histograms exactly for any chemical potential and smearing in range; `bands` and `berry_curvature` are evaluated on request. The evaluator follows the archived WannierTools 2.6.1 protocol (atomic-gauge Kubo formula, Fermi–Dirac occupations); the Néel-reversed model is the exact time reverse of the archived Hamiltonian.
- `make_task.py` stages the archive files, writes the catalog, the lab data, the hidden targets, the ablation values and `calibration.json`.
- Staged unchanged: the VASP inputs and outputs (`1_vasp`), the Wannier90 inputs and outputs except `wannier90_hr.dat` (`2_wannier90`), all WannSymm files (`3_wannsymm`), `wt.in` and `sigma_ahc.gnu` (`4_wtools`), and the two Berry-curvature maps with their gnuplot script (`fig6_script`). `4_wtools/WT.out` is staged without its 40,000 k-point progress lines. No pseudopotential is in the archive or the task.

## What is hidden

Graded values are stored only in `environment/lab/private/hidden.json`:

- **σ_xy at hole densities 3.0 × 10²¹ and 6.0 × 10²¹ cm⁻³:** the archived WannierTools output `4_wtools/sigma_ahc.txt` (801 × 801 × 401 mesh, WannSymm-symmetrized Hamiltonian, energies relative to `E_FERMI` of `wt.in`), read at the chemical potentials where the lab evaluator's 400-mesh electron count of the archived Hamiltonian gives these densities.
- **Sign changes of σ_xy between 1.00 and 0.25 eV below the VBM, and the hole density of the sign change:** from the archived curve, the density from the 400-mesh electron count.
- **Hole density 0.30 eV below the VBM:** the lab evaluator (400 mesh); the VBM is the maximum of Wannier band 22 found by minimisation.
- **Forbidden and allowed components:** the magnetic symmetry of the archived state (the archived output has |σ_yz|, |σ_zx| < 10⁻⁷ S/cm); electron doping: not determined by the model (`null`).

The lab never serves `sigma_ahc.txt` or the hidden values. The archived Hamiltonian is not the symmetrized one the authors used (that file is not archived); at the graded hole densities the lab's σ_xy agrees with the archived output within the tolerances (see the construction README).
