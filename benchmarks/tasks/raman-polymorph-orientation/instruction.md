# SnF2 crystal: polymorph, orientation and unmeasured geometries from polarized Raman

A single crystal of tin(II) fluoride, SnF2, with one polished face (face A) sits in a polarized micro-Raman spectrometer at 295 K. Determine which polymorph it is, whether a second phase is present, how the crystal is oriented, and what a measurement on a second face (face B), which cannot be measured, would show.

## Data and calculation service

`/app/catalog.json` holds:

- the archived records of the three SnF2 structures in the Computational Raman Database (Bagheri & Komsa, Sci. Data 10, 80 (2023); PBEsol phonons): cell, positions, masses, Γ-point frequencies, irreducible representations, eigenvectors, Raman tensors, Born charges, ε∞ and quality flags;
- the database conventions, the lab frame, the description of face B and the instrument.

The lab serves backscattering spectra from face A. They are simulated from a physical forward model of the hidden sample built from these records; treat them as measurements. Each request records one 320 cm⁻¹ window of the spectrum:

```sh
python3 /app/lab_client.py evaluate '{"kind": "raman", "polarizer_deg": 0, "analyzer_deg": 90, "center_cm1": 200, "accumulations": 1}'
python3 /app/lab_client.py status
```

- `polarizer_deg` and `analyzer_deg` are angles in the lab xy plane; `"analyzer_deg": null` removes the analyzer.
- `center_cm1` is the window centre on the spectrometer's Raman-shift axis.
- `accumulations` (1–4) adds 60 s exposures.

Each accumulation costs 1 unit; the budget is 20 units. Results print as JSON; redirect them to a file.

## Conventions

- **Lab frame:** z is the outward normal of face A, so the phonon wavevector in face-A backscattering is along z.
- **Frequency scale factor s:** the sample's phonon frequencies at 295 K are s times the database (PBEsol) harmonic frequencies of the same polymorph, for every phonon (the dynamical matrix scaled by s²).
- **c axis:** the unit vector, in lab coordinates, along the third lattice vector of the main phase's database cell. Its sign is arbitrary.
- **Face B:** perpendicular to face A, with outward normal n_B = (cos 30°, sin 30°, 0) and in-plane direction t_B = (−sin 30°, cos 30°, 0). In face-B backscattering the phonon wavevector is along n_B.
- **Strongest mode:** the Raman-active mode of the main phase with the largest Σ_αβ R_αβ² in its database record, degenerate partners summed.

## Deliverables

Save `/app/submission.json`:

```json
{"main_phase": "mp-0000",
 "minority_phase": "mp-0000",
 "frequency_scale": 0.0,
 "c_axis_lab": [0.0, 0.0, 1.0],
 "faceB_strongest_mode_ratio_tt_over_zz": 0.0,
 "faceB_highest_optical_phonon_cm1": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`main_phase`:** the `mpid` of the crystal's polymorph.
- **`minority_phase`:** the `mpid` of a second polymorph present in the sample, or `"none"`.
- **`frequency_scale`:** s for the main phase.
- **`c_axis_lab`:** the c axis of the crystal.
- **`faceB_strongest_mode_ratio_tt_over_zz`:** for the strongest mode in backscattering from face B, the intensity with incident and analysed polarizations both along t_B divided by that with both along z.
- **`faceB_highest_optical_phonon_cm1`:** the highest zone-centre optical-phonon frequency (cm⁻¹, sample at 295 K) for a phonon wavevector along n_B.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the polymorphs were compared and why the others were rejected;
- how the spectra were calibrated and how s was obtained;
- the scattering model, and how the orientation was determined and shown to be unique;
- how the face-B values were computed;
- the evidence for or against a second phase;
- uncertainties, and which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the phases, s, the c axis and the two face-B values; physical reasoning is reviewed separately.
