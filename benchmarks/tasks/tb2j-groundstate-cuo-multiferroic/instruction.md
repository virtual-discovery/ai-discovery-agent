# CuO: which magnetic ground state do magnetic-force-theorem exchange parameters predict?

Using released TB2J exchange parameters for CuO (tenorite), determine the magnetic ground state predicted by the self-consistent workflow of Tellez-Mora et al. (2024): compute the exchange in a reference magnetic state, test whether that state is stable in linear spin-wave theory with its own exchange, and update the reference until it is. Establish how the reference state and TB2J's spin-ligand correction change the prediction, and characterise the spin waves of the converged state.

## Data and calculation service

`/app/catalog.json` describes four released TB2J calculations: a ferromagnetic (FM) or antiferromagnetic (AFM) reference, each with the spin-ligand correction on or off. It also gives the conventional cell, the exchange convention and what can be purchased. Purchases return archived author output; they are not new calculations.

```sh
python3 /app/lab_client.py evaluate '{"kind": "exchange", "calculation": "FM/ligand", "window": "0-4"}'
python3 /app/lab_client.py evaluate '{"kind": "author_dispersion", "calculation": "AFM/ligand"}'
python3 /app/lab_client.py status
```

- **`exchange`:** every isotropic pair of one calculation within one Cu–Cu distance window, together with that calculation's cell, atoms, charges and moments. A window costs 1–3 units; the catalog lists the windows and their costs.
- **`author_dispersion`:** the authors' archived eigenvalues of their bosonic matrix h(k) for one calculation, along their k-path. Costs 1 unit.

The budget is 24 units. Results are large: redirect them to files.

## Conventions

- **Exchange:** the archived TB2J convention, E = −Σ_{i≠j} J_ij ê_i·ê_j over ordered pairs, with unit vectors along the local moments, J > 0 ferromagnetic, in meV.
- **Wave vectors:** reciprocal-lattice units of the conventional cell given in the catalog. Vectors related by a sign change or by adding integer vectors are treated as equivalent.
- **Stability:** a calculation is stable at its own reference when the linear spin-wave bosonic matrix, built from that calculation's exchange around its own collinear moment configuration, is positive semi-definite at every k.
- **k_min:** for an FM reference, the wave vector at which the lowest eigenvalue of that matrix is smallest (the k_min of the paper).
- **Chain exchange:** the dominant coupling, quoted for H = Σ_bonds J S_i·S_j with each bond counted once, spin-1/2 operators and J > 0 antiferromagnetic, using the classical correspondence ê = S/S.
- **Spin waves:** linear spin-wave theory around the converged collinear state with spin length S = 1/2 (ê_i = S_i/S). Use the archived pairwise exchange of that calculation in its own cell, as the authors' workflow does. The mean magnon energy is the average over the Brillouin zone of all branches; the band maximum is the highest magnon energy anywhere in the zone.

## Deliverables

Save `/app/submission.json`:

```json
{"kmin_FM_reference_rlu": {"ligand_corrected": [0.0, 0.0, 0.0], "uncorrected": [0.0, 0.0, 0.0]},
 "stable_at_own_reference": {"FM/ligand": false, "AFM/ligand": false, "FM/no_ligand": false, "AFM/no_ligand": false},
 "converged_propagation_vector_rlu": [0.0, 0.0, 0.0],
 "chain_exchange_meV": 0.0,
 "mean_magnon_energy_meV": 0.0,
 "magnon_band_maximum_meV": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`kmin_FM_reference_rlu`:** k_min of the FM-reference calculation with and without the spin-ligand correction, to three decimals.
- **`stable_at_own_reference`:** the stability verdict for each of the four calculations.
- **`converged_propagation_vector_rlu`:** the propagation vector of the ground state predicted by the self-consistent (converged) calculation.
- **`chain_exchange_meV`, `mean_magnon_energy_meV`, `magnon_band_maximum_meV`:** for the converged calculation.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the conventions and how the calculations' cells and moment configurations were handled;
- how each calculation was tested, and how the converged one was identified;
- the effect of the reference state and of the spin-ligand correction, with numbers;
- how the spin-wave calculation was set up and validated;
- what an isotropic exchange model can and cannot say about CuO's incommensurate multiferroic phase;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the k_min vectors, the stability verdicts, the converged propagation vector, the chain exchange and the two magnon energies; physical reasoning is reviewed separately.
