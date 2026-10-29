# Exfoliable 1D chalcogenide chains: quasiparticle gaps and excitons from released GW-BSE data

Single atomic chains exfoliable from van der Waals crystals (S3, Te3, As2S3 and Bi2Te3) are predicted to have exciton binding energies of order eV. Using the authors' released DFT, GW and Bethe–Salpeter (BSE) results, determine why the binding is so strong, how much eigenvalue self-consistency and spin–orbit coupling (SOC) change the quasiparticle (QP) gaps, how much the screened interaction used in the BSE kernel matters, and whether these results can be transferred to Bi2Te3, whose self-consistent results are withheld.

## Data and calculation service

`/app/catalog.json` describes the four chains, the archived calculations, the request formats and the methods. Purchases return the authors' archived results; nothing is recomputed:

- **`qp`:** rows of a Yambo QP table (G0W0 or scGW level) at chosen k-points, with the run header;
- **`bse_spectrum`:** an archived Haydock BSE spectrum with its run header (two runs per chain, under the authors' labels);
- **`dft_bands`:** archived PBE bands along Γ–Z, computed with or without SOC;
- **`wire_screening`:** the DFPT dielectric tensor of the supercell and the quantum-volume record of the chain.

```sh
python3 /app/lab_client.py evaluate '{"kind": "qp", "chain": "Te3", "level": "scGW", "k_indices": [20, 21, 22]}'
python3 /app/lab_client.py evaluate '{"kind": "bse_spectrum", "chain": "S3", "run": "BSE_scGW"}'
python3 /app/lab_client.py evaluate '{"kind": "dft_bands", "chain": "Bi2Te3", "soc": false}'
python3 /app/lab_client.py evaluate '{"kind": "wire_screening", "chain": "As2S3"}'
python3 /app/lab_client.py status
```

Costs are 1 unit per k-point of a `qp` request, 3 units per `bse_spectrum`, and 1 unit per `dft_bands` or `wire_screening` request. The budget is 45 units. The scGW QP table and the `BSE_scGW` spectrum of Bi2Te3 are not served. Results print as JSON; redirect large ones to a file.

## Conventions

- Energies are in eV. **scGW** is the last archived iteration of the eigenvalue self-consistent GW scheme (catalog `levels`); **G0W0** is the one-shot calculation on PBE.
- **Direct QP gap:** the minimum over the 25 archived k-points of E_c(k) − E_v(k), where v is the band `valence_band_max_index` of the catalog and c the band above it.
- **With and without SOC:** `with_soc` is the value of a calculation that includes SOC; `without_soc` is the same quantity without SOC. Where a gap has to be carried from one case to the other, assume that the QP corrections do not depend on SOC, so that SOC changes the gap exactly as it changes the minimum direct gap of the archived DFT band structures (each band structure's own minimum along Γ–Z).
- **Optical gap:** the energy of the first bright exciton, that is, the lowest-energy local maximum of Im α in an archived BSE spectrum whose height is at least 10% of the largest Im α of that spectrum.
- **scGW optical gap:** the optical gap of a BSE whose QP energies and screened interaction W both come from the scGW calculation.
- **Exciton binding energy:** E_b = (scGW direct QP gap, with SOC) − (scGW optical gap).
- **Screening update shift Δ_W:** the change of the optical gap when the screened interaction of the BSE kernel is switched from the DFT-level W0 of the G0W0 calculation to the scGW-level W, with the QP energies held at their scGW values. Negative means the exciton moves down.

## Deliverables

Save `/app/submission.json`:

```json
{"scgw_direct_gap_eV": {"Te3": {"with_soc": 0.0, "without_soc": 0.0},
                        "Bi2Te3": {"with_soc": 0.0, "without_soc": 0.0}},
 "scgw_optical_gap_eV": {"S3": 0.0, "Te3": 0.0, "As2S3": 0.0},
 "screening_update_shift_eV": {"S3": 0.0, "Te3": 0.0, "As2S3": 0.0},
 "bi2te3_binding_energy_eV": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`scgw_direct_gap_eV`:** the scGW direct QP gap of Te3 and Bi2Te3, with and without SOC.
- **`scgw_optical_gap_eV`:** the scGW optical gaps of S3, Te3 and As2S3.
- **`screening_update_shift_eV`:** Δ_W for S3, Te3 and As2S3.
- **`bi2te3_binding_energy_eV`:** your best estimate of E_b for Bi2Te3.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- why the exciton binding energies of these chains are so large, and how the QP gap opening and the binding compare;
- how much self-consistency and SOC change the QP gaps, chain by chain, and how SOC was treated;
- how each BSE run was interpreted and how Δ_W was obtained;
- which excitons are bright, and how the optical gaps were identified;
- how the Bi2Te3 values were obtained, which transfer models you tested, and their uncertainty;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the gaps, the optical gaps, the screening update shifts and the Bi2Te3 binding energy; physical reasoning is reviewed separately.
