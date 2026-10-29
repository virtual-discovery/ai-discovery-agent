# Construction of gwbse-1d-chain-excitons

This is a released-start task built from Grillo, Cignarella, Bechstedt, Gori, Palummo and Campi, *Quasiparticle effects and strong excitonic features in exfoliable 1D semiconducting materials*, ACS Nano 20, 2664 (2026). It uses Materials Cloud record 42t0v-jam38 (DOI 10.24435/materialscloud:xq-dy, CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. Run `python3 -B fetch.py [DEST]`. It downloads `1D_ACS_Nano_Exciton.zip` (3.3 MB, md5 a3cdbff7…) and `README.txt` through the archive's per-file links, then extracts and verifies the 45 members listed in `members.md5`. The default DEST is `~/.cache/ai-theorist/gwbse-1d-chain-excitons`, which `model.py` reads.
2. Run `python3 -B make_task.py`. It writes:
   - `environment/agent/catalog.json`;
   - `environment/lab/data.json` (1.4 MB);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files.

   It asserts the archive facts the task relies on and prints the calibration numbers. `model.py` holds the readers and the estimators.
3. Run `.venv/bin/python -m benchmarks.tools.check_task gwbse-1d-chain-excitons --cpuset <cpus>`. It writes `check.json`.
4. Run `python3 -B record_oracle.py`, which copies the gate's oracle run into `oracle_run.json`. Then run `python3 -B write_design.py`, which writes `design.json` with the calibration table.

## What the task is

The agent buys archived results for four chains: S3, Te3, As2S3 and Bi2Te3. The requests and costs are:

- Yambo QP-table rows at chosen k-points (G0W0 or scGW): 1 unit per k-point;
- Haydock BSE spectra (`BSE_G0W0` and `BSE_scGW`, the authors' labels): 3 units;
- PBE bands with or without SOC: 1 unit;
- the DFPT supercell dielectric tensor and the quantum volume: 1 unit.

The budget is 45 units; the catalogue is worth 208. The scGW QP table and the `BSE_scGW` spectrum of Bi2Te3 are withheld.

The graded quantities are:

- the scGW direct gap of Te3 and Bi2Te3, with and without SOC;
- the scGW optical gaps of S3, Te3 and As2S3;
- the screening update shift Δ_W of S3, Te3 and As2S3;
- the Bi2Te3 binding energy, graded by a regime window.

## What the archive shows (and the task is built on)

- **The GW and BSE include SOC.** The Eo column of every QP table reproduces the SOC band structure, within 4 meV on all 25 k-points; the SOC-free structure differs by up to 0.52 eV. The paired spinor bands are split away from Γ. The design assumed GW without SOC plus a DFT correction, so that step became the double-counting trap dp3. The BSE headers print "[SOC] Perturbative SOC: no", which reads as if SOC were absent.
- **The `BSE_G0W0` runs of Te3, As2S3 and Bi2Te3 used the scGW QP energies.** Their headers name the scGW database (`KfnQPdb= "E < G3W3/ndb.QP"`, or G4W4), and their independent-particle spectra equal those of the `BSE_scGW` runs to within a relative 5e-6. Only the kernel differs (input files `haydock.in` and `haydock_W.in`). S3's `BSE_G0W0` run alone used G0W0 energies. The archive's README does not mention this.
- **Self-consistency is material dependent.** It opens the direct gap by 0.792, 0.806 and 0.748 eV for S3, Te3 and As2S3, but by 0.111 eV for Bi2Te3. The Bi2Te3 GW used 700 bands, against 1450–1600 for the others.
- **Dark and weak transitions.** The band-edge transitions of S3, Te3 and As2S3 are weak along the chain. As2S3's first bright exciton (2.706 eV) belongs to higher transitions; weak peaks at 1.86–2.66 eV carry 1–2% of the maximum.

## Deviations from the design (candidates.json / DETAILS.md section 32)

- **Withheld chain.** The withheld chain is Bi2Te3, but its G0W0 QP table, its `BSE_G0W0` run, its bands and its screening data are served; only its scGW table and `BSE_scGW` spectrum are withheld. The design withheld everything but DFT-level data. That version is not uniquely recoverable:
  - a 1D screened-wire model with the DFPT polarizability and r0 from the quantum volume overestimates E_b for every chain (S3 3.03 against 2.27; Te3 2.31 against 2.14), so its calibration is chain dependent;
  - Bi2Te3's self-consistency increment (0.11 eV) differs from the others' (0.75–0.81 eV).

  The Bi2Te3 gap is now graded tightly (±0.05 eV) through evidence in Bi2Te3's own data. Its binding energy is graded by a regime window, and the model transfers are the wrong paths.
- **Targets changed.**
  - "scGW QP gap at the exciton k" became the minimum direct gap, with and without SOC, for Te3 and Bi2Te3.
  - "First bright BSE peak" became the scGW optical gaps of S3, Te3 and As2S3, with a stated 10% brightness rule.
  - The Bi2Te3 binding energy became a window.
  - The G0W0-to-scGW increment "predicted before purchase" was dropped, since it cannot be enforced. Δ_W was added: it tests the run-label inconsistency directly.
- **Decision points.** The design's dp1 (direct gap versus indirect fundamental gap) was dropped: the two differ by only 0.006–0.09 eV here, within the gap tolerance. The design's dp2 (mixing levels) is now part of dp1 and dp4, through the mislabelled runs. The design's dp3 (SOC) became double counting.
- **Not staged:** phonons and IR spectra, which are not needed for the question.

## Decision points (6) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | `BSE_G0W0` label versus the QP database in the header and the IP spectrum (Te3, As2S3, Bi2Te3) | metadata inconsistency, found in raw outputs |
| dp2 | S3's runs do differ in QP energies, so the Te3 finding must not be generalized | the same evidence, checked per chain |
| dp3 | The archived GW already includes SOC; adding the DFT SOC change double counts | evidence against a textbook habit, plus a misleading header line |
| dp4 | The Bi2Te3 scGW gap: transferred increment (2.32 eV) versus the IP spectrum of its own scGW-energy run (1.65 eV) | regime judgement; a discriminating purchase |
| dp5 | First bright exciton versus weak features (As2S3) | spectrum reading |
| dp6 | Bi2Te3 E_b: its own W0 binding with a W-update correction, versus a 1D wire or hydrogenic model or the mean of the other chains | model class and transfer |

dp1, dp2, dp3 and dp4 rest on evidence that has to be found in the data. dp4 also cascades into the Bi2Te3 binding energy.

## Calibration

The full table is in `design.json` under `calibration`; values are in eV.

| Field | Reference | Oracle | Defensible variants | Nearest wrong path | Tolerance |
|---|---|---|---|---|---|
| Te3 gap with SOC | 4.310 | 4.310 | 4.274 (fundamental gap) | 4.026 (SOC added again) | ±0.05 |
| Te3 gap without SOC | 4.595 | 4.595 | 4.590 (k-resolved), 4.631 (SOC change at fixed k) | 4.310 (SOC added again) | ±0.05 |
| Bi2Te3 gap with SOC | 1.649 | 1.661 | 1.636 (IP onset at 1%), 1.661 (first IP peak) | 1.561 (SOC added again), 1.539 (G0W0 gap), 2.321 (transferred increment) | ±0.05 |
| Bi2Te3 gap without SOC | 1.737 | 1.749 | 1.724, 1.749, 1.767 (SOC change at fixed k) | 1.649 (SOC added again), 1.627 (G0W0) | ±0.05 |
| optical gap S3 / Te3 / As2S3 | 4.067 / 2.171 / 2.706 | exact | the same peaks for brightness thresholds of 3–30% | As2S3 1.861 (first feature) | ±0.05 |
| Δ_W S3 | −0.347 | −0.347 | −0.335 (IP-peak shift) | +0.445 (runs assumed to share energies) | ±0.10 |
| Δ_W Te3 | −0.545 | −0.545 | none | −1.351 (labels trusted) | ±0.10 |
| Δ_W As2S3 | −0.120 | −0.120 | none | −0.868 (labels trusted), −0.510 (first feature) | ±0.10 |
| Bi2Te3 E_b | 0.854 | 0.844 | 0.755 (no W correction), 0.844 (ratio scaled by self-consistency), 1.026 (mean ratio), 1.093 (additive mean shift) | 0.633 (labels trusted), 1.317 (wire model calibrated on S3 and Te3), 1.573 (wire model), 1.584 (mean of released chains) | window [0.68, 1.15] |

The reference definitions:

- **Gaps:** minimum over the 25 k-points of the archived scGW table. The SOC-free gap subtracts the SOC change of the minimum direct DFT gap: Te3 −0.284, Bi2Te3 −0.088.
- **Optical gaps:** the first local maximum of Im α with at least 10% of the spectrum maximum.
- **Δ_W:** the difference of the two runs' optical gaps. For S3, the G0W0-to-scGW change of the minimum direct QP gap (0.792) is subtracted.
- **E_b:** Bi2Te3 scGW gap minus the first bright exciton of the withheld `BSE_scGW` spectrum (0.795).

The Bi2Te3 window is a regime check. The honest spread of W-update corrections, 0.76–1.09 eV, is wider than any tolerance that would discriminate finer.

## Oracle

The oracle makes 19 evaluations for 39 of 45 units (87%):

- bands with and without SOC for Te3 and Bi2Te3 (4 units);
- QP rows (13 units): Te3 scGW at k 20–22 and G0W0 at k 21; S3 G0W0 and scGW at k 22–23; As2S3 G0W0 and scGW at k 1; Bi2Te3 G0W0 at k 1–3;
- all seven purchasable BSE spectra (21 units);
- the Bi2Te3 screening record (1 unit).

It does the following:

1. verifies that the Eo energies match the SOC bands;
2. reads the KfnQPdb line of every BSE header and checks it against the IP spectra;
3. takes Δ_W directly for Te3 and As2S3, and after removing the QP shift for S3;
4. reads the Bi2Te3 scGW gap from the first IP peak of its `BSE_G0W0` run;
5. estimates the Bi2Te3 E_b by scaling its W0 binding energy with the released chains' W-update ratio, taken proportional to the relative self-consistent gap opening.

The estimated agent length is 40–70 tool calls.

Results:

- **Gate (2026-10-19, cpuset 6,7,22,23):** passed. All 10 ablations fail, the 9 robustness variants pass, and rounding to 4 significant digits passes. The oracle took 9 s in the gate container. There are 6 decision points and 7 rubric items, and the package is 1.5 MB with 13 hidden values.
- **Harbor oracle (`runs/build-gwbse-1d-chain-excitons/oracle-r1-gwbse-1d-chain-excitons`):** reward 1.0, runtime 50 s.

## Open risks

- **Dev-easy risk.** Agents that read the headers will find dp1, dp2 and dp4 at once: the `KfnQPdb` line names `G3W3`, and the catalog lists `G3W3` as the scGW iteration. dp3 needs a cross-check of Eo against the two band structures. Difficulty depends on whether agents read headers and question the "Perturbative SOC: no" line.
- **Kernel interpretation.** That the `BSE_G0W0` runs of Te3, As2S3 and Bi2Te3 use the DFT-level W0 is inferred from the authors' labels and input names; the archive does not document it. The QP-energy finding is certain. If the kernels were other than W0 and W(scGW), the Δ_W definition would still match the archived run pairs, but its physical reading would change.
- **Thin margins.**
  - The Bi2Te3 SOC double-counting path fails by 0.038 eV; the IP-based estimates pass by about 0.036.
  - The Te3 gap read from the first IP peak (4.362, above the weak band-edge transition) fails by 0.002 eV. It is not a defensible estimate of the defined minimum direct gap, since the Te3 QP table is purchasable.
  - The convention deviations "fundamental gap" and "SOC change at fixed k" pass by 0.014 eV.
  - The window edges sit 0.047 eV from the labels-trusted E_b (0.633) and 0.057 eV from the additive-mean variant (1.093).
- **Brightness rule.** The As2S3 W0 run has a 6.4%-high peak at 2.581 eV. A 5% rule would select it, but the rule is stated as 10%.
- **Wire model calibrated on all three released chains.** Including As2S3, whose bright exciton is not its band-edge exciton, lands at about 0.97 eV, inside the window. The ablation uses the S3 and Te3 calibration.
- **Memorisation.** The paper's abstract quotes the optical gaps: 4.07, 2.17, 2.71 and 0.8 eV. The first three can also be read from purchasable data. The 0.8 eV value helps the Bi2Te3 E_b only together with the evidence-based gap.

## Disk

- **Peak:** about 25 MB in `~/.cache/ai-theorist/gwbse-1d-chain-excitons` (the zip plus its extraction) and a 9.5 MB test fetch in the scratchpad. Both were deleted after the gate and the Harbor run.
- **Kept:** the refetch script with checksums.
- **Task package:** 1.5 MB.
