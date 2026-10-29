# MnTe: which exchange makes the chiral magnons, and how field cooling switches them

Hexagonal MnTe is a collinear compensated magnet whose two Mn sublattices are related by a screw rotation, not by a translation or inversion. Liu, Masuda and co-workers measured its magnons with polarized neutrons on two crystals: half-polarized time-of-flight data on HYSPEC (crystal S1) and full polarization analysis on IN20 (crystal S2, prepared in two field-cooled states A and B). Using the archived data, determine which exchange interactions split the magnon branches and give them opposite chirality, how large that exchange is, where in the Brillouin zone the splitting is largest, and predict the IN20 inelastic measurements that the lab withholds.

## Data and calculation service

`/app/catalog.json` describes:

- the crystal and magnetic structure;
- a prior Hamiltonian from earlier unpolarized measurements, which leaves one piece undetermined;
- the HYSPEC geometry, including the energy transfer of every slice pixel;
- the IN20 conventions;
- the headers of the withheld IN20 scans (everything above `DATA_:`);
- the request formats and costs.

A purchase returns archived measurements:

```sh
python3 /app/lab_client.py evaluate '{"kind": "hyspec_tile", "quantity": "symmetric", "tile": "u_low_v_low"}' > sym_ll.json
python3 /app/lab_client.py evaluate '{"kind": "hyspec_cut", "cut": "g"}'
python3 /app/lab_client.py evaluate '{"kind": "pnd", "instrument": "IN20", "state": "A"}' > pnd_A.json
python3 /app/lab_client.py status
```

- **HYSPEC slices:** the symmetric and the chiral slice are each sold as four tiles (quadrants in u and v, listed in the catalog), 1 unit per tile.
- **HYSPEC cuts:** each of the four cuts (g, h, i, j) costs 1 unit.
- **Polarized diffraction:** the HYSPEC scan, or the IN20 raw file of state A or state B, costs 1 unit each.

The budget is 8 units. These are retrieved author measurements, not new experiments.

## Conventions

- **Hamiltonian:** H = Σ over Mn pairs (each pair once) J S_i·S_j + D Σ_i (S_i^z)², with S = 5/2, J > 0 antiferromagnetic, z ∥ c, and energies in meV. This is the convention of the prior Hamiltonian.
- **Wave vectors:** r.l.u. of the hexagonal cell given in the catalog (a = 4.148 Å, c = 6.711 Å).
- **Polarization:** I_α^{if} is the intensity with polarization axis α in the Blume–Maleyev frame (x ∥ Q, y in the scattering plane, z = x × y) and incident/scattered spin states i, f.

## Deliverables

Save `/app/submission.json`:

```json
{"split_exchange_difference_meV": 0.0,
 "in20_peak_qh": {"low": 0.0, "high": 0.0},
 "in20_chirality": {"A": {"low": "+", "high": "+"}, "B": {"low": "+", "high": "+"}},
 "max_splitting_q_rlu": [0.0, 0.0, 0.0],
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`split_exchange_difference_meV`:** |J_a − J_b|, where J_a and J_b are the exchange constants of the two symmetry classes into which the 9.83 Å same-sublattice bonds fall; their average is the prior value.
- **`in20_peak_qh`:** the two magnon peaks in the withheld 32 meV constant-energy scans of the IN20 files. Give the peak centres in the QH coordinate exactly as the instrument records it in those files, for the summed I_x^{+-} and I_x^{-+} intensity of both states. "low" is the peak at smaller QH.
- **`in20_chirality`:** for states A and B, the sign ("+" or "−", written "-") of I_x^{-+} − I_x^{+-} at the low-QH and at the high-QH peak.
- **`max_splitting_q_rlu`:** a wave vector (h, k, l) at which the altermagnetic splitting between the two magnon branches of opposite chirality is largest.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which exchange paths are split by symmetry, and why;
- where the splitting vanishes and where it is largest;
- how the exchange difference was obtained from the data, and how robust it is;
- how the IN20 coordinates and energies were related to your model;
- how the chirality of each branch and each field-cooled state was inferred, including the role of polarization channels and magnetic domains;
- why relativistic terms or magnon–phonon mixing do not explain the observations;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks:

- the exchange difference;
- the mean position and the separation of the IN20 peaks;
- the four chirality signs;
- the wave vector of maximal splitting.

Physical reasoning is reviewed separately.
