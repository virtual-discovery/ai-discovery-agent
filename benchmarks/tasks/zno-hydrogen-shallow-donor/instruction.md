# ZnO particles: which defect is the g = 1.96 shallow donor, and how much of it forms during synthesis?

ZnO particles made by gas-phase synthesis, and annealed in vacuum at 873 K afterwards, show an EPR resonance at g = 1.96 from a paramagnetic shallow donor. The candidates are the oxygen vacancy, the zinc interstitial, hydrogen on or at an oxygen site, and hydrogen-containing complexes. From the authors' released defect calculations, determine:

- which defect, in which charge state, produces the signal;
- where the donor levels lie relative to the band edges;
- the equilibrium Fermi level and defect concentrations at the synthesis temperature.

## Data and calculation service

`/app/catalog.json` describes the released archive (VASP, 72-atom wurtzite supercells) and lists every calculation under its archive folder name:

- the defect-free supercell;
- the defect supercells of seven defects in several charge states;
- the reference phases.

Each folder holds a relaxation, and its `HSE06` subfolder holds a hybrid-functional single point on the relaxed structure.

There are two kinds of purchase:

- **`run`:** returns the run settings, the final energies, the total and site magnetization, the Fermi energy, the initial and final structures, the k-points with their weights, and VASP's per-atom average electrostatic potential at the cores.
- **`eigenvalues`:** returns the final Kohn–Sham eigenvalues and occupations of both spin channels at every k-point.

```sh
python3 /app/lab_client.py evaluate '{"kind": "run", "id": "supercell_3_3_2/OH_O/q_0/HSE06"}'
python3 /app/lab_client.py evaluate '{"kind": "eigenvalues", "id": "supercell_3_3_2/HSE06"}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit; the budget is 35 units. These are retrieved author calculations, not new ones. Results print as JSON; redirect them to files.

## Conventions

- **Level of theory:** use one Hamiltonian throughout, with defect energies, band edges and reference phases all from the same functional.
- **Formation energy:** E_f = E(defect cell) − E(perfect cell) − Σ n_i μ_i + q(E_F + E_VBM) + corrections. Here n_i is the number of atoms of species i added, q is the charge in units of e (the electron count of the neutral defect cell minus that of the cell), and E_F is measured from the VBM. Use the energies extrapolated to σ → 0.
- **Band edges:** the VBM, the CBM and the gap come from the defect-free supercell.
- **Chemical potentials:** use the 0 K total energies of the archived reference phases, with no zero-point or thermal terms.
  - **Zinc and oxygen:** take the Zn-rich limit, μ_Zn = E(Zn metal) per atom and μ_O = E(ZnO) − μ_Zn per formula unit.
  - **Hydrogen:** μ_H is the largest value at which none of the archived hydrogen-containing reference phases is lower in energy than its constituents at these chemical potentials.
  - **Report:** Δμ_O = μ_O − ½E(O₂) and Δμ_H = μ_H − ½E(H₂).
- **Dilute limit:** formation energies and levels refer to an isolated defect in an infinite crystal. Where a cell carries a charge localized at the defect, apply the Lany–Zunger image-charge term (2/3 of the point-charge Madelung energy of the supercell lattice, ε = 10.4) and a potential alignment. The alignment is q·ΔV, where ΔV is the mean change of the core potentials of the atoms farther than 5 Å from the defect.
- **Transition levels:** these are the thermodynamic levels between consecutive thermodynamically stable charge states, in eV above the VBM. Report the levels that lie in the gap. In each entry, `from` is the stable charge just below the level and `to` the one just above it.
- **Equilibrium at T = 1000 K:**
  - **Defects:** every archived defect species has one site per formula unit of ZnO (36 per supercell), and every charge state has degeneracy 1. The concentration of each state is c = N_sites exp(−E_f/k_BT), with k_B = 8.617333 × 10⁻⁵ eV/K.
  - **Free carriers:** use parabolic bands with m_e* = 0.24 and m_h* = 0.59 (in units of m₀), Fermi–Dirac statistics, and the supercell gap.
  - **Fermi level:** E_F is set by charge neutrality over all archived defects and charge states.
- **Hydrogen site of `OH_O`:** classify it from the relaxed structure of the neutral cell as one of `BC_parallel`, `BC_perpendicular`, `AB_parallel` or `AB_perpendicular`. BC is the bond centre of a Zn–O bond and AB the antibonding site on the far side of the O atom; parallel and perpendicular refer to the O–H direction relative to the c axis.

## Deliverables

Save `/app/submission.json`:

```json
{"band_gap_eV": 0.0,
 "delta_mu_eV": {"O": 0.0, "H": 0.0},
 "transition_levels": {"OH_O": [], "H_O": [], "V_O": [{"from": 0, "to": -1, "level_eV": 0.0}]},
 "V_O_plus_ever_stable": true,
 "epr_centre": {"defect": "X", "charge": 0},
 "oh_o_site": "AB_perpendicular",
 "fermi_level_eV": 0.0,
 "log10_concentration_cm3": {"V_O": 0.0, "OH_O^+": 0.0, "OH_O^0": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`V_O_plus_ever_stable`:** whether V_O⁺ is the lowest-energy charge state of the oxygen vacancy anywhere in the gap.
- **`epr_centre`:** the defect (named as in the archive folders) and charge state responsible for the g = 1.96 signal.
- **`log10_concentration_cm3`:** the 1000 K equilibrium concentrations, as log₁₀ of cm⁻³. `V_O` is the total over its charge states; `OH_O^+` and `OH_O^0` are the two charge states of OH_O.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which calculations you used, and why;
- how the chemical potentials were fixed;
- how each charge state was brought to the dilute limit, and whether its carriers are bound to the defect;
- the level ladders and which defect states can be paramagnetic in equilibrium;
- the assignment of the EPR signal, and the competing explanations it rules out;
- the equilibrium Fermi level and concentrations, and what they imply for unintentional hydrogen doping;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks:

- the gap, the chemical potentials, and the OH_O, H_O and V_O level ladders;
- the V_O⁺ stability, the EPR assignment and the hydrogen site;
- the Fermi level and the concentrations.

Physical reasoning is reviewed separately.
