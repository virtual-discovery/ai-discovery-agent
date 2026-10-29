# Pb, Nb and Al: thermal expansion versus anharmonicity in the phonon-limited resistivity

How do thermal expansion and phonon anharmonicity separately change the phonon-limited electrical resistivity of lead, niobium and aluminium at high temperature, and what resistivity results when both are included? Answer from the authors' archived first-principles calculations (DFPT, TDEP and EPW). Several of their results are withheld and must be estimated.

## Data and calculation service

`/app/catalog.json` describes the archive (its README is quoted), the methods, the lattice parameters, the request formats and two pools:

- **`transport_runs`:** the archived EPW Boltzmann-transport calculations, under the archive's folder names. Each entry gives its condition and temperatures, and states whether it is served.
- **`force_constant_sets`:** every archived harmonic force-constant set, DFPT at each lattice parameter and TDEP at each sampling temperature and lattice.

You can make two kinds of request:

- **`transport` (2 units):** returns the archived run: its input parameters, fine meshes, and for each temperature the self-energy relaxation-time (SERTA) and iterative (IBTE) conductivity tensors.
- **`phonons` (1 unit):** returns phonons computed by the lab from one force-constant set: frequencies on a uniform grid (with polarisation vectors on request), a DOS histogram, standard moments and high-symmetry dispersions.

```sh
python3 /app/lab_client.py evaluate '{"kind": "transport", "run": "Pb/4-epw/1-DFPT/epw"}'
python3 /app/lab_client.py evaluate '{"kind": "phonons", "ifc": "Pb:TDEP410K:a(0K)", "grid": 16}'
python3 /app/lab_client.py status
```

The budget is 46 units. A running calculation returns a job id; poll it with `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file. Transport results are retrieved author calculations; no new EPW calculation is possible.

## Conventions

- **Resistivity ρ(T):** 1/σ̄ in μΩ cm, where σ̄ is one third of the trace of the converged iterative-BTE conductivity tensor at temperature T.
- **Conditions:**
  - *harmonic:* DFPT force constants at the 0 K lattice parameter;
  - *thermal expansion (TE):* DFPT force constants at the lattice parameter of T;
  - *anharmonic at the 0 K lattice:* TDEP force constants sampled at T, at the 0 K lattice parameter;
  - *TE + anharmonicity:* TDEP force constants sampled at T, at the lattice parameter of T.
- **Change due to a condition X:** Δ_X(T) = 100 [ρ_X(T)/ρ_harmonic(T) − 1] percent. The harmonic reference is at the same T and uses the same numerical settings as the X calculation (fine k and q meshes, Fermi window, smearing).

## Deliverables

Save `/app/submission.json`:

```json
{"al_te_change_percent_481K": 0.0,
 "pb_te_resistivity_610K_uohm_cm": 0.0,
 "pb_anharmonic_change_percent_0K_lattice": {"410": 0.0, "610": 0.0},
 "pb_te_anharm_resistivity_510K_uohm_cm": 0.0,
 "nb_te_anharm_resistivity_2000K_uohm_cm": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`al_te_change_percent_481K`:** Δ_TE of Al at 481 K.
- **`pb_te_resistivity_610K_uohm_cm`:** ρ of Pb with thermal expansion at 610 K.
- **`pb_anharmonic_change_percent_0K_lattice`:** Δ of Pb for anharmonicity at the 0 K lattice, at 410 and 610 K.
- **`pb_te_anharm_resistivity_510K_uohm_cm`:** ρ of Pb with TE + anharmonicity at 510 K.
- **`nb_te_anharm_resistivity_2000K_uohm_cm`:** ρ of Nb with TE + anharmonicity at 2000 K.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how each withheld value was estimated, and how the method was validated on released results;
- what thermal expansion and anharmonicity each do to the phonons and the resistivity of each metal, and why;
- whether the two effects combine additively or multiplicatively;
- how the three metals differ;
- in which temperature range the data and your methods apply;
- the uncertainties, and which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the six values; physical reasoning is reviewed separately.
