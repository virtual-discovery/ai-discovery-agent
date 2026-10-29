# Irradiated GaN: which native defects set the deep levels, and where are they?

Irradiation of wurtzite GaN creates vacancies, interstitials, antisites and divacancies on both sublattices, and carbon is a common impurity. From author-released hybrid-functional defect calculations, determine the thermodynamic charge-transition levels of the Ga and N vacancies, the Ga and N interstitials and the C_N impurity in the dilute limit, and the energetics of the defect pairs that irradiation produces: the binding of the V_Ga–V_N divacancy and the formation energies of Ga and N Frenkel pairs.

## Data and calculation service

`/app/catalog.json` lists the calculations of the archive (VASP, 96-atom orthorhombic supercells, HSE-type hybrid functional) under the authors' folder names, with the authors' description of the archive. The defect-free supercell runs are `perfect` and its subfolders. Each defect folder is one defect, site label and charge state; for example `Va_Ga1_-2` is a Ga vacancy with q = −2, and `Ga_i1`/`Ga_i2` are two Ga-interstitial configurations. The archive contains no elemental reference calculations.

A purchase returns one released calculation:

- the run settings printed in its OUTCAR, the PAW datasets, the k-point, the electron count and the ionic-step count;
- the final energies and Fermi energy;
- the initial and final structures;
- the final Kohn–Sham eigenvalues and occupations of both spin channels;
- the per-atom "average (electrostatic) potential at core" and spin moments.

A separate purchase returns dielectric tensors of GaN computed by the task builders; use them wherever your analysis needs dielectric constants.

```sh
python3 /app/lab_client.py evaluate '{"kind": "calculation", "id": "Va_Ga1_-2"}'
python3 /app/lab_client.py evaluate '{"kind": "dielectric"}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit; the budget is 43 units. These are retrieved author calculations, not new ones. Results print as JSON; redirect them to files.

## Conventions

- **Charge:** q is in units of e: the electron count of the neutral defect cell minus that of the cell.
- **Formation energies** refer to an isolated defect in an otherwise perfect, infinite crystal (the dilute limit), in the lowest-energy structure available for each charge state. The Fermi level E_F is measured from the valence-band maximum (VBM).
- **Band edges:** the VBM and conduction-band minimum (CBM) of bulk GaN described by the same Hamiltonian as the defect calculations (functional, PAW datasets, cutoff). `band_gap_eV` is this gap.
- **Transition levels:** thermodynamic levels ε(q/q′) between the VBM and the CBM, in eV above the VBM, in increasing order. Each lies between two consecutive thermodynamically stable charge states: `from` is the stable charge just below the level and `to` the one just above it.
- **Divacancy binding:** E_b = E_f(V_Ga) + E_f(V_N) − E_f(V_Ga–V_N) at E_F = E_CBM − 0.5 eV, each defect in its lowest-energy charge state at that E_F. Positive means bound.
- **Frenkel pairs:** E_FP = E_f(vacancy) + E_f(interstitial) of the same element (N: V_N + N_i; Ga: V_Ga + Ga_i), isolated and non-interacting, each in its lowest-energy charge state and structure at E_F = E_CBM − 0.5 eV.

Chemical potentials cancel in E_b and E_FP.

## Deliverables

Save `/app/submission.json`:

```json
{"band_gap_eV": 0.0,
 "transition_levels": {"V_Ga": [{"from": 0, "to": -1, "level_eV": 0.0}],
                       "V_N": [], "Ga_i": [], "N_i": [], "C_N": []},
 "divacancy_binding_energy_eV": 0.0,
 "frenkel_pair_energy_eV": {"N": 0.0, "Ga": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`transition_levels`:** every level of each of the five defects in the gap. Levels within 0.3 eV of a band edge are not graded, and a charge state stable over less than 0.1 eV of E_F is graded as unstable.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which calculations define the host reference and the band edges, and why;
- how the finite-size electrostatics were treated, with the size of the corrections;
- which charge states and structures are physical ground states, and whether any charge state's extra electrons or holes are not bound to the defect;
- the level ladders, including any negative-U behaviour;
- which defects can compensate n-type GaN or pin the Fermi level after irradiation, and the role of C_N;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the band gap, the transition levels, the divacancy binding energy and the Frenkel-pair energies; physical reasoning is reviewed separately.
