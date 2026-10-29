# Palladium hydride: why it phase-separates, and what deuterium changes

Below a critical temperature, hydrogen in palladium separates into a dilute α phase and a concentrated β hydride. Using the authors' released machine-learned potential for Pd–H (an atomic cluster expansion fitted to PBE data), establish what drives the separation and what sets the critical temperature and the plateau pressure. Then establish how hydrogen sits in the lattice and how replacing H by D changes solubility.

The competing pictures are:

- short-range attractive H–H chemistry;
- a long-range, volume-mediated (elastic) attraction, with electronic band filling limiting the hydride composition;
- site preference (octahedral vs tetrahedral) together with zero-point motion.

## Calculation service

`/app/catalog.json` describes the potential, the request formats, the gas-phase data you need, and a list of archived DFT records you can retrieve for validation. Structures are periodic, in Cartesian Å, with at most 300 atoms.

```sh
python3 /app/lab_client.py evaluate '{"kind": "energy", "structures": [<structure>, ...]}'
python3 /app/lab_client.py evaluate '{"kind": "relax", "structure": <structure>, "mode": "volume"}'
python3 /app/lab_client.py evaluate '{"kind": "scan", "structure": <structure>, "index": 108, "grid": {"center": [x, y, z], "half_width": 0.5, "spacing": 0.05}}'
python3 /app/lab_client.py evaluate '{"kind": "dft_reference", "ids": ["DFT-001"]}'
python3 /app/lab_client.py job/<id>
python3 /app/lab_client.py status
```

- **`energy`:** up to 8 structures per request, 1 unit.
- **`relax`:** 1 unit. The `mode` is one of:
  - `positions`, with the cell fixed;
  - `volume`, isotropic scaling of the cell to zero pressure;
  - `cell`, full relaxation of the cell to zero stress.
- **`scan`:** 2 units. It returns the total energy with one atom placed at each of up to 10 000 points, given as an explicit list or as a cubic grid, while every other atom stays frozen.
- **`dft_reference`:** up to 4 archived DFT records, 1 unit.

The budget is 30 units. Long requests return a job id; poll it and do not resubmit. Results print as JSON; redirect them to files.

## Conventions

- **Hydrogen reference:** ½ E(H₂), where E(H₂) is the static energy of the isolated relaxed molecule computed with the same potential. Gas-phase chemical potentials are μ°(X₂, 300 K, 1 bar) = E_static(X₂) + Δμ_X₂. The catalog supplies Δμ_X₂ and the molecular zero-point energy (ZPE) for H₂ and D₂; use those values, not your own gas-phase vibrational treatment.
- **Dilute limit:** one interstitial in an infinite Pd crystal at zero pressure. Interstitial sites are octahedral (O, one per Pd) and tetrahedral (T, two per Pd).
- **Nuclear motion:** treat the motion of H and D quantum-mechanically, within the Born–Oppenheimer approximation, in the potential-energy surface of the relaxed host with the Pd atoms frozen. Neglect Pd vibrations.
- **Configurational model** (used for T_c and the plateau):
  - H occupies O sites at random, with ideal configurational entropy per Pd.
  - Its energy Δe(x), per Pd, is the static energy of homogeneous random occupations, each fully relaxed (positions and cell, zero stress), relative to the end members. Here x = H/Pd.
  - The incoherent critical temperature is the top of the spinodal of Δe(x) − T s_ideal(x).
  - The coherent critical temperature is the top of Cahn's coherent spinodal for a cubic crystal. Use the elastic constants of pure Pd at its equilibrium lattice parameter and the dilute-limit Vegard coefficient (1/a)(da/dx).
  - Vibrational contributions are left out of T_c.
- **Plateau pressure:** the equilibrium (incoherent) α/β coexistence at 300 K in the configurational model. Each H additionally carries the vibrational free energy of an isolated O-site interstitial, independent of x.
- **Masses:** as given in the catalog.

## Deliverables

Save `/app/submission.json`:

```json
{"solution_energy_0K_eV": 0.0,
 "OT_difference_static_eV": 0.0,
 "OT_difference_with_zpe_eV": 0.0,
 "lattice_expansion_A3_per_H": 0.0,
 "tetrahedral_fraction_H_300K": 0.0,
 "sieverts_ratio_H_over_D_300K": 0.0,
 "log10_plateau_pressure_H_300K_bar": 0.0,
 "Tc_incoherent_K": 0.0,
 "Tc_coherent_K": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`solution_energy_0K_eV`:** the dilute O-site solution energy of H at 0 K, relative to ½H₂, with the zero-point energies of H in the metal and of the molecule included.
- **`OT_difference_static_eV` and `OT_difference_with_zpe_eV`:** the dilute energy of H on a T site minus that on an O site, static and with the ZPE of H included.
- **`lattice_expansion_A3_per_H`:** the dilute-limit volume increase per dissolved H atom.
- **`tetrahedral_fraction_H_300K`:** the fraction of dissolved H on T sites at 300 K in the dilute limit.
- **`sieverts_ratio_H_over_D_300K`:** the ratio of the dilute equilibrium concentrations of H and D at 300 K and equal gas pressure.
- **`log10_plateau_pressure_H_300K_bar`:** the base-10 logarithm of the H₂ plateau pressure at 300 K, in bar.
- **`evidence_ids`:** the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- what drives the phase separation and what limits it, and the evidence;
- how you obtained each number and established that it is converged, including how you treated the hydrogen nuclei;
- how far you trust the potential, and the evidence;
- the isotope effect and its origin;
- the difference between the coherent and incoherent gaps;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the nine numbers; physical reasoning is reviewed separately.
