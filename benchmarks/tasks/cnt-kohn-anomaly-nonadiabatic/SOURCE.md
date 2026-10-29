# Source and transformations

- **Study:** E. Marazzi, S. Poncé, J.-C. Charlier and G.-M. Rignanese, *Importance of nonadiabatic effects in Kohn anomalies in 1D metals*, Phys. Rev. Lett. 135, 256402 (2025), DOI 10.1103/x69r-8fw1 (arXiv:2506.18704).
- **Archive:** Materials Cloud Archive record vs5fd-dbv98, *Importance of non-adiabatic effects on Kohn anomalies in 1D metals*, DOI 10.24435/materialscloud:yx-x1, https://archive.materialscloud.org/records/vs5fd-dbv98.
- **License:** CC BY 4.0 (all files).

## Files used

| File | Bytes | md5 |
|---|---|---|
| `README.md` | 1,237 | 1c6dbfd0eb10564470d215e3eca226fe |
| `dft_inputs.zip` | 4,991 | c85fa7a8af0ffa84597e7fe76e1749b5 |
| `Fig1.zip` | 215,112 | 4c1760d1aaf80add6900cd65233ca571 |
| `Fig2.zip` | 26,803,427 | b17a42c46d5a7856356d71eab691875c |

The md5 of every member used is listed in `benchmarks/construction/cnt-kohn-anomaly-nonadiabatic/members.md5`. The prebuilt binary `Fig1/coupling_calculator` is not used.

## Construction

The scripts are in `benchmarks/construction/cnt-kohn-anomaly-nonadiabatic/`:

- `fetch.py` downloads the four files through the record's per-file links and verifies them.
- `extract.py` reads, without modification:
  - the 24 coupling tables `Fig2/g/g±q.<n m ν>`: 2000 k-points; band energies at k and k±Q; ω_Q; |g|;
  - the 37 self-energy tables `Fig2/PiT/<T>`;
  - `e_band/band.eig` and `ph_band/matdyn.freq`;
  - the DFPT smearing scan, which is the arrays `SIGMA_VALUES` and `OMEGA_EQ1` of `Fig2/fig2_d_f_h.py`;
  - one line of `fig2_a_b.py` (`FERMI_SHIFT`) and one comment of `fig2_c_e_g.py` (summing k+q and k−q), quoted verbatim;
  - the small text files of `dft_inputs/`, `Fig1/` and `README.md`.
- `make_task.py` stages the data and computes the references. It first checks that the lab evaluator `selfenergy.py` reproduces the archived tables (maximum deviation 0.17 meV over modes 1, 20 and 23 at 50–1578.87 K).

## What the lab serves

- **Verbatim text files:** `README.md`, the files of `dft_inputs/`, and `Fig1/coupling_calculator.c`, `Fig1/fig1.py` and `Fig1/backup/output.csv`.
- **Parsed arrays:** the bands, the phonon dispersion, the coupling tables and the DFPT scan.
- **Archived self-energy tables:** for ν = 1, 20 and 23 only, as verbatim lines with the file's own header, at 50, 700, 800, 900, 1000, 1250 and 1578.87 K.
- **New self-energy evaluations:** computed by `selfenergy.py` from the archived coupling tables. The bands and |g|² are interpolated along branches tracked through the band crossings. The chemical potential is fixed at the archived EPW value, −1.0380 eV.

The plotting scripts `fig2_*.py` are not served. They contain the authors' renormalization protocol: the subtraction of the static part of the 1578.87 K table and the Eq. (6) iteration.

## Archive features the task relies on (verified)

- **Self-energy table header.** The tables have 8 columns, but the header names 10, in another order: `... w  Real Sigma(w)  Imag Sigma(w)  Real Sigma(w=0)  Chi_d(w)  Chi_(w=0)`.
  - Column 7 is independent of ω: it is the static Re Π(0).
  - Column 8 is Im Π(ω).
  - The authors' script reads the columns that way.
- **Smearing schemes.** The DFPT run behind ω_Q (`scf.in`) uses `degauss = 0.01` Ry with Marzari–Vanderbilt smearing (ω_Q = 35.28, 130.11 and 143.42 meV). The archived DFPT smearing scan uses Fermi–Dirac smearing: at 10 mRy it gives 48.93 and 137.83 meV for ν = 1 and 20. The authors' protocol subtracts the Fermi–Dirac static self-energy at 1578.87 K from the cold-smearing ω_Q.
- **Threshold program.** In `coupling_calculator.c` the parameter called `mass` in `main()` is `a = 1/m*`, as the function documentation says, and `fig1.py` labels the output rows with the reversed list. The program reproduces `output.csv` exactly. The paper text quotes thresholds of 54.5, 129.0 and 137.5 meV for the CNT modes. These match neither the program (60.33, 151.17, 163.40) nor the paper's Eq. (9) (62.19, 155.91, 167.40).
- **Fermi level.** The band-plot script shifts by −1.0426 eV, but the self-energy input uses −1.0380 eV, which is the Dirac crossing. The self-energies are insensitive to this difference at the target temperatures.

## What is hidden

The 30 self-energy tables at 70–600 K are not served. From the archived data, with the authors' protocol, the hidden references are:

- the complex-pole frequency and half width of ν = 20 at 100 K and 300 K;
- the adiabatic frequency of ν = 20 at 100 K;
- the RBM instability temperature, from the Ω → 0 limit of the pole across the 500–600 K tables;
- the threshold couplings g₀ of the archived program at m* = 2 and ω_Q of ν = 1 and ν = 20.

They are stored only in `environment/lab/private/`.

## Scope

This is a retrospective reconstruction from released data, not a new DFPT or EPW calculation. The lab evaluator is a numpy re-implementation of the EPW phonon self-energy sum on the archived coupling tables; no EPW code was copied.
