# Kohn anomaly in a 4 Å carbon nanotube: nonadiabatic softening and the RBM instability

The (3,3) carbon nanotube is a one-dimensional metal. Its Fermi points are connected by a nesting wavevector Q, where several phonon branches show Kohn anomalies. Using the authors' released first-principles electron–phonon data, determine whether the Kohn anomaly drives a lattice instability, for which mode and below which electronic temperature, and how nonadiabatic (dynamic) phonon self-energy effects change the softening.

## Data and calculation service

`/app/catalog.json` describes the released data (LDA; Quantum ESPRESSO 7.3 and EPW; Materials Cloud record vs5fd-dbv98) and the exact request formats. You can purchase:

- **`archive_file`, `bands`, `phonons`, `coupling`, `dfpt_scan`:** archived input files and model code, the electronic and phonon band structures, the electron–phonon coupling tables at Q, and the authors' DFPT smearing scan.
- **`self_energy_table`:** the archived EPW phonon self-energy Π_ν(ω) at Q for one mode, at one of the released temperatures.
- **`self_energy`:** a numpy evaluator of Π_ν(ω, T) at Q from the archived coupling tables, at any temperature. You choose the k-point sampling, the broadening, the scattering processes and the band pairs; the catalog gives the formula.

```sh
python3 /app/lab_client.py evaluate '{"kind": "self_energy_table", "mode": 20, "temperature_K": 700}'
python3 /app/lab_client.py evaluate '{"kind": "archive_file", "file": "dft_inputs/scf.in"}'
python3 /app/lab_client.py status
```

A `self_energy` evaluation costs 1 unit for nk ≤ 10⁵ k-points, 2 units for nk ≤ 10⁶ and 3 units above; every other request costs 1 unit. The budget is 26 units. A running calculation returns a job id; poll it with `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **Modes** are numbered at Q from the lowest frequency, as in the archive: ν = 1 is the radial breathing mode (RBM) and ν = 20 is an in-plane optical (A1) mode.
- **Self-energy** Π_ν(ω, T) is in meV, with the sign convention of the archive and the lab: Re Π < 0 softens the mode and Im Π(ω > 0) < 0. T is the electronic temperature. Treat any electronic smearing σ used in the archived calculations as a Fermi–Dirac temperature, k_B T = σ.
- **ω_Q** is the adiabatic DFPT frequency at Q from the archived phonon calculation that underlies the coupling and self-energy data (the value listed with them).
- **Renormalized phonon.** Ω − iγ (Ω > 0 the frequency, γ the half width) is the complex pole of the diagonal phonon propagator at Q, with the self-energy evaluated on the real axis at ω = Ω:

  (Ω − iγ)² = ω_b² + 2 ω_Q Π_ν(Ω, T),

  where ω_b is the bare frequency of the mode (without any electronic self-energy), to be obtained from ω_Q. The **adiabatic** frequency is the same pole with Π_ν(Ω, T) replaced by the static Π_ν(0, T), so γ = 0.
- **RBM instability temperature T\*:** the temperature at which, on cooling from high temperature, Ω of the ν = 1 mode first reaches zero (just below T\*, the mode has no real positive Ω).
- **Threshold coupling g₀:** the coupling |g_Q| at which the authors' analytic one-dimensional model (the archived `Fig1/coupling_calculator.c`, at T = 0 with its archived settings) gives Ω → 0. Evaluate it for effective mass m* = 2 (in free-electron masses) and the ω_Q of the mode.

## Deliverables

Save `/app/submission.json`:

```json
{"omega_nonadiabatic_meV": {"100": 0.0, "300": 0.0},
 "linewidth_meV": {"100": 0.0, "300": 0.0},
 "omega_adiabatic_100K_meV": 0.0,
 "rbm_instability_temperature_K": 0.0,
 "threshold_coupling_meV": {"1": 0.0, "20": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`omega_nonadiabatic_meV`, `linewidth_meV`:** Ω and γ of mode ν = 20 at T = 100 K and 300 K.
- **`omega_adiabatic_100K_meV`:** the adiabatic Ω of mode ν = 20 at 100 K.
- **`rbm_instability_temperature_K`:** T\* of the RBM, or `null` if the RBM stays stable at every T ≥ 50 K.
- **`threshold_coupling_meV`:** g₀ for ν = 1 and ν = 20.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how ω_b was obtained from the archived data;
- how the self-energies at the target temperatures were obtained and validated (sampling, broadening, processes, band pairs);
- the adiabatic versus nonadiabatic softening of ν = 20 and its physical origin;
- how T\* was determined, and how much nonadiabatic effects change it;
- how g₀ compares with the computed couplings at the Fermi level, and what that implies for stability;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the eight numbers above; physical reasoning is reviewed separately.
