# Co3Sn2S2: one canted kagome ferromagnet, or a second in-plane antiferromagnetic phase?

Co3Sn2S2 is a kagome Weyl ferromagnet (T_C ≈ 177 K, moments along c). Several studies have proposed an additional in-plane antiferromagnetic phase at intermediate temperatures. Tam et al. instead describe a single noncollinear ground state produced by a large Dzyaloshinskii–Moriya interaction (DMI). Use their released diffuse neutron scattering and spin Hamiltonian to decide what the critical magnetic fluctuations reveal, predict the measurements the lab withholds, and quantify what the DMI does in their model.

## Data and calculation service

`/app/catalog.json` describes the crystal, the scattering geometry, the array format and what can be purchased. `/app/archive/` holds the archive's `README.md`, `readcombined.py` and the SpinW script `SpinW_Co3Sn2S2_DMI.m`. The data are single-crystal maps in the (H 0 L) plane from the DMC diffractometer, stored as four channels per 0.015 Å⁻¹ pixel (counts, monitor, normalization, normcount):

- **Narrow scans** around (0 0 3) between 140 and 250 K, sold whole.
- **Wide scans** at 1.5, 165, 180 and 300 K, sold as square regions (±0.48 Å⁻¹) around listed (h 0 l) positions.

```sh
python3 /app/lab_client.py evaluate '{"kind": "narrow_scan", "temperature_K": 192}' > scan192.json
python3 /app/lab_client.py evaluate '{"kind": "wide_region", "temperature_K": 300, "hkl": [1, 0, 4]}' > w300_104.json
python3 /app/lab_client.py status
```

A narrow scan costs 2 units and a wide-scan region 1 unit; the budget is 16 units. Results are large, so save them to files. These are retrieved author measurements, not new experiments. The narrow scans at 170, 175 and 200 K and the (1 0 1) region of the 180 K wide scan are withheld.

## Conventions

**Annular profile** of a scan around a reflection (h 0 l):

- A pixel is valid if normcount > 0, monitor > 0 and normalization is finite and positive. Its normalized intensity is counts·normcount/(monitor·normalization), and pixel centres lie half a bin (0.0075 Å⁻¹) above the lower edges.
- The peak pixel is the valid pixel of largest normalized intensity whose centre lies within 0.05 Å⁻¹ of the nominal (h 0 l) position. The 7 × 7 block of pixels centred on it is excluded.
- Annulus k contains the remaining valid pixels whose integer pixel offsets (Δi, Δj) from the peak pixel satisfy e_k² ≤ Δi² + Δj² < e_{k+1}², with e = (4, 6, 8, 10, 12, 15, 17, 20, 27) pixels: 8 annuli.
- The annulus value is Σ counts / Σ (monitor·normalization/normcount) over its pixels, in units of 10⁻⁶, and its uncertainty is √(Σ counts) / Σ (monitor·normalization/normcount) in the same units.
- A predicted profile is scored by the reduced χ² against the withheld scan: the mean over the 8 annuli of ((prediction − value)/uncertainty)², with no free scale or offset.

**Spin Hamiltonian:** the model that `SpinW_Co3Sn2S2_DMI.m` defines when run in SpinW 3.1.1.

- Energy E = Σ over bonds of S_iᵀ J S_j, each bond counted once, with the spin length given in the script. Cartesian axes are x ∥ a and z ∥ c.
- Bonds are indexed as SpinW's `gencoupling` numbers them for the script's lattice.
- The J1dmi matrix in the script is the coupling of the bond from the Co at (0, ½, ½) to the Co at (½, ½, ½) (hexagonal fractional coordinates); SpinW generates the couplings of the other bonds with the same index from it.
- **`canting_angle_deg`:** the angle in degrees between the Co moments and the c axis in the classical (T = 0) ground state.
- **`spin_wave_gap_meV`:** the lowest linear-spin-wave energy at the zone centre about that ground state, in meV.

## Deliverables

Save `/app/submission.json`:

```json
{"narrow_003_profiles": {"170": [0, 0, 0, 0, 0, 0, 0, 0], "175": [0, 0, 0, 0, 0, 0, 0, 0], "200": [0, 0, 0, 0, 0, 0, 0, 0]},
 "wide_101_profile_180K": [0, 0, 0, 0, 0, 0, 0, 0],
 "canting_angle_deg": 0.0,
 "spin_wave_gap_meV": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`narrow_003_profiles`:** the predicted (0 0 3) annular profiles of the withheld narrow scans.
- **`wide_101_profile_180K`:** the predicted (1 0 1) annular profile of the 180 K wide scan.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the data were reduced and the background treated;
- how the withheld profiles were predicted, and which purchases constrained them;
- which spin components and which type of correlation the diffuse scattering shows near T_C;
- whether these data support a separate intermediate-temperature in-plane antiferromagnetic phase;
- the size of the DMI relative to J1 in the authors' model, and how it produces the canting and the gap;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the four predicted profiles, the canting angle and the gap; physical reasoning is reviewed separately.
