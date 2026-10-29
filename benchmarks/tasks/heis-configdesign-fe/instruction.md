# bcc Fe: a transferable Heisenberg exchange model and where it breaks down

Build a Heisenberg exchange model for body-centred-cubic iron from author-released noncollinear DFT energies, test how far it transfers, and establish where the Heisenberg description of this itinerant magnet fails.

## Data and calculation service

`/app/catalog.json` describes two pools of released DFT configurations (PBE, VASP, constrained moments):

- **Spin spirals** in the primitive cell (a = 2.87 Å), listed by direction and q in reciprocal-lattice coordinates of that cell.
- **Random noncollinear states** in several supercells, each given with its cell, atomic positions and spin directions.

Each pool has configurations you may purchase (`*_available`) and configurations whose energies you must predict (`*_to_predict`). Purchases return the released energy and local moments:

```sh
python3 /app/lab_client.py evaluate '{"kind": "spiral", "key": "[111]-03"}'
python3 /app/lab_client.py evaluate '{"kind": "random_state", "label": "236-05"}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit; the budget is 60 units. Energies are meV/atom relative to the collinear ferromagnet. These are retrieved author calculations, not new DFT.

## Conventions

Report exchange constants for H = −Σ_{i≠j} J_ij ê_i·ê_j: unit vectors, every ordered pair counted, so that E/atom − E_FM = (1/N) Σ_i Σ_j J_ij (1 − ê_i·ê_j).

The spin-wave stiffness D is defined by E(q) − E_FM = D q² per atom as q → 0, for a flat spiral (moments perpendicular to the spiral axis), with q in Å⁻¹ including the factor 2π.

## Deliverables

Save `/app/submission.json`:

```json
{"spiral_predictions_meV_per_atom": {"S214-01": 0.0},
 "random_predictions_meV_per_atom": {"R334-01": 0.0},
 "J1_meV": 0.0,
 "spin_wave_stiffness_meV_A2": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **Predictions:** give an energy for every id in `spirals_to_predict` and `random_states_to_predict`.
- **`J1_meV`:** the nearest-neighbour exchange of the model that describes the spin-spiral (small-angle) energy landscape at a = 2.87 Å.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model and conventions;
- how the configurations and the number of exchange shells were chosen and validated;
- how D was obtained;
- where and why the Heisenberg description fails;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the predictions, J1 and D; physical reasoning is reviewed separately.
