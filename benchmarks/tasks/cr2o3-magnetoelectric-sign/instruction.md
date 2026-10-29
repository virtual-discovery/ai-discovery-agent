# Cr2O3: sign and size of the lattice-mediated magnetoelectric response for a given antiferromagnetic domain

Cr2O3 is the prototypical linear magnetoelectric. Its two antiferromagnetic domains have opposite magnetoelectric tensors, and the literature has long disagreed about which domain has which sign. From the authors' archived first-principles calculations (VASP and Elk, with phonopy and Quantum ESPRESSO inputs), determine the zero-kelvin lattice-mediated spin contribution to the transverse coefficient α⊥ for each domain, which polar modes carry it, and which domain is selected by magnetoelectric annealing in parallel electric and magnetic fields.

## Data and calculation service

`/app/catalog.json` describes the archive (record xd06z-0pp29), includes the archive README and lists every item you can purchase. Each purchase costs 1 unit and returns one archived result:

- `born_reference` and `born_displacement`: the finite-displacement Berry-phase runs (LCALCPOL) used for Born effective charges. A `born_displacement` purchase returns the three runs (x, y, z) of one atom at one displacement amplitude.
- `elk_moment` and `vasp_moment`: one self-consistent noncollinear calculation (spin–orbit coupling, LSDA+U) with one polar eigenvector of the force-constant matrix frozen in at one amplitude.
- `force_constant_eigenmodes`, `phonon_gamma` (phonopy runs with real masses, all-Cr or all-H masses) and `structure`.
- `qe_inputs`: input files of the Quantum ESPRESSO calculations (no outputs are archived).

```sh
python3 /app/lab_client.py evaluate '{"kind": "born_displacement", "amplitude": "0_02", "atom": "O3"}'
python3 /app/lab_client.py evaluate '{"kind": "elk_moment", "mode": "Eu_2_x", "amplitude": "-0_01"}'
python3 /app/lab_client.py status
```

The budget is 44 units; the catalog holds far more items than that. These are retrieved author calculations, not new ones. Results print as JSON; redirect them to files.

## Conventions

- **Frame.** Cartesian axes of the Born-charge, phonon and induced-moment calculations (cell vector a₁ = (2.826, 0, 4.498) Å); the trigonal axis is z. α⊥ = α_xx = α_yy.
- **Magnetoelectric tensor.** α_ij = μ₀ ∂M_i/∂E_j in SI (report ps/m), where M is the net magnetic moment per unit volume. Only the lattice-mediated spin contribution is wanted: the field acts on the ions through their Born effective charges, the ions relax according to the force-constant matrix, and the response is the net spin moment induced by that static displacement. Electrons respond only through the frozen-in structure.
- **Moments.** Treat the moments printed by VASP and Elk as magnetic moments (not spins), as the authors do.
- **Domains.** In the *out-pointing* domain the magnetic moments of the two nearest-neighbour Cr ions (the closest Cr–Cr pair) point away from each other; in the *in-pointing* domain they point towards each other.
- **Modes.** `Eu_1_x` … `Eu_4_x` are the archived eigenvectors of the force-constant matrix with polarization along x, in the archive's order and sign. The contribution of mode m to α⊥ is its term in the sum over these modes; its fraction is that term divided by α⊥.
- **Born charges.** Z*_{κ,αβ} = (Ω/e) ∂P_α/∂u_{κβ}, with α the polarization component and β the displacement component. Cr1 is the first Cr atom in POSCAR order.
- **Temperature dependence (experiment, single domain).** α⊥ keeps one sign at all temperatures below T_N = 307 K. At low temperature α∥ has the same sign as α⊥; α∥ changes sign on warming (near 100 K) and keeps the opposite sign up to T_N.

## Deliverables

Save `/app/submission.json`:

```json
{"alpha_perp_ps_per_m": {"out_pointing": 0.0, "in_pointing": 0.0},
 "mode_fractions": {"Eu_1_x": 0.0, "Eu_2_x": 0.0, "Eu_3_x": 0.0, "Eu_4_x": 0.0},
 "dominant_mode": "Eu_1_x",
 "dM_du_out_pointing_muB_per_A": {"Eu_3_x": 0.0, "Eu_4_x": 0.0},
 "born_charge_Cr1": {"xx": 0.0, "xy": 0.0, "zz": 0.0},
 "dataset_domain": {"VASP": "out_pointing", "Elk": "out_pointing"},
 "domain_selected_by_parallel_field_annealing": "out_pointing",
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`alpha_perp_ps_per_m`:** the zero-kelvin lattice-mediated spin α⊥ of each domain.
- **`mode_fractions`** and **`dominant_mode`:** the fraction of α⊥ carried by each mode, and the mode with the largest fraction.
- **`dM_du_out_pointing_muB_per_A`:** the derivative of the net magnetic moment along x per rhombohedral cell (μB) with respect to the amplitude (Å) of the archived eigenvector of that mode, in the out-pointing domain.
- **`born_charge_Cr1`:** Z*_xx, Z*_xy and Z*_zz of Cr1.
- **`dataset_domain`:** the domain modelled in the induced-moment calculations of each code.
- **`domain_selected_by_parallel_field_annealing`:** the domain obtained by cooling through T_N in parallel E and H fields along the trigonal axis.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the formula you used, the mode basis and the units, including the cell volume;
- how the Born charges were constructed and checked;
- how the moment derivatives were obtained and which runs you trusted, and why;
- how the domain of each dataset was established, and how the codes compare after that;
- which modes carry the response and why;
- the temperature argument for the annealed domain;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks α⊥ of both domains, the mode fractions and dominant mode, the two moment derivatives, the Cr1 Born charge, the domain of each dataset and the annealed domain; physical reasoning is reviewed separately.
