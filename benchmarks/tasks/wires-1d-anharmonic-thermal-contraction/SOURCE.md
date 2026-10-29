# Source and transformations

- **Source:** C. Cignarella, L. Bastonero, L. Monacelli and N. Marzari, *Extreme anharmonicity and thermal contraction of 1D wires*, Nano Lett. (2025), DOI 10.1021/acs.nanolett.5c04282; arXiv:2508.07971.
- **Archive:** Materials Cloud Archive record k45dq-ebf35, DOI 10.24435/materialscloud:fj-20, https://archive.materialscloud.org/records/k45dq-ebf35.
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `Extreme-anharmonicity-thermal-contration-1D-wires.zip`: 145,671,625 bytes, md5 1de3fb44caf32b1ed403cbe0e7b687ce. Only its `SSCHA/` and `QHA/` folders are used; the ML potentials are not.
  - `README.txt`: md5 c866d88afdf77c93202d26d1c30e0b38.

## Construction

The scripts are in `benchmarks/construction/wires-1d-anharmonic-thermal-contraction/`.

- `fetch.py` downloads and checks the two files.
- `extract.py` builds `pool.json`:
  - For every SSCHA run (wire, temperature, 10,000 or 100,000 configurations), it takes the auxiliary (`final_sscha_*`) and free-energy Hessian (`hessian_*`) dynamical-matrix files verbatim.
  - It takes a verbatim excerpt of the minimization log: the last minimization step and the final RESULTS block of that temperature section. Configuration-progress and debug lines are removed.
  - It takes the DFPT initial guesses as ph.x wrote them (`*.dyn0-N`) and with the authors' 1D sum rule applied (`*_second_asr_*`).
  - It takes the DFT input/output files in `initial_guess/`.
  - From `QHA/qha_1D.aiida` (read with sqlite3, without aiida-core), it takes the inputs, the pw.x output and the DYN_MAT files of the single archived strained TaSe3 calculation.
- `make_task.py` writes the lab data (the served texts), the catalog, the private targets and the ablation values. It computes them with `calib.py`, `model.py` and `dyn.py`.

## What is hidden

- **Withheld files:** the TaSe3 free-energy Hessians at 650 K and 700 K (10,000 and 100,000 configurations).
- **Graded values,** computed from the archive with the conventions stated in the instruction:
  - the TaSe3 harmonic axial stiffness: the long-wavelength LA branch of the 1D-ASR DFPT force constants, averaged over q = 0.005–0.02 (2π/c);
  - the TaSe3 α_L at 300 K: the SSCHA stress slope from the 200–400 K runs, times the cell volume, divided by the stiffness;
  - the CuC2 and TaSe3 C_V ratios at 300 K: the SSCHA-entropy derivative from the auxiliary frequencies, using mode-by-mode quadratic fits over the converged 100–600 K runs;
  - the TaSe3 instability temperature: the midpoint of the 650–700 K ω² crossings from the withheld 10k and 100k Hessians.

The archive contains no α_L or C_V output files, so no archived number is withheld for those quantities. The graded values do not reproduce the paper's figures, for two reasons:

- the paper's C_V formula (Eq. 2 and the archived notebook) carries the opposite sign of the frequency-derivative term;
- the paper's α_L values are about twice those obtained with k = ∂²E/∂ε².

## Scope

This is a retrospective reconstruction from released files, not a new calculation. The QHA thermal expansion cannot be reconstructed: the QHA archive holds only one strained TaSe3 calculation.
