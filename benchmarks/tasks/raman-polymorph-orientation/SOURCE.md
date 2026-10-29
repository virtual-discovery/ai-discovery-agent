# Source and transformations

- **Source:** M. Bagheri and H.-P. Komsa, *High-throughput computation of Raman spectra from first principles*, Sci. Data 10, 80 (2023), DOI 10.1038/s41597-023-01988-5.
- **Archive:** Computational Raman Database v1.0.1, Materials Cloud Archive record 0m9r0-g3486, DOI 10.24435/materialscloud:pg-h3, https://archive.materialscloud.org/records/0m9r0-g3486 (published 2024-08-05).
- **License:** CC BY 4.0 (whole record). Phonon data in the database are adapted by the authors from phonondb (A. Togo), also CC BY 4.0.
- **File used:** `Computational_Raman_Database_v1.0.1.zip` (152,753,693 bytes, md5 e2aaea03df70fda95dbd970c60456882). From its member `crd.json` (454,304,470 bytes, an ASE-database JSON document with 5099 records) only three records were used, unchanged:

| record | mpid | space group | sha256 of the record (canonical JSON, sorted keys) |
|---|---|---|---|
| 3216 | mp-2843 | C2/c (α-SnF2) | da176f4d3a76e91f399bd8ca8a1404f0d92c238d31555169d5b961917c860d4c |
| 2564 | mp-7456 | P2_12_12_1 (β-SnF2) | 849d4bfea4de74b2e6022c51bfbd4bf999ee1e1bccea4ccd59f3465c2999c376 |
| 2390 | mp-7457 | P4_12_12 (γ-SnF2) | e5735d7dc160445293f7135b0afd19a52d10710ce85306e1a26d42a6b06f7193 |

## What the agent sees

`environment/agent/catalog.json` holds the three records with their archived content: cell, positions, frequencies, irrep labels and degenerate sets, Raman/IR flags, Raman tensors, eigenvectors, Born charges and ε∞, plus e_above_hull, the negative-frequency flags and the band gap. Two transformations were made:

- floats were rounded to 6 significant digits;
- the archived strings and python-literal lists were converted to JSON arrays.

Per-atom masses (the standard atomic weights phonopy uses) were added, because the database stores only the total mass. The Born charges are staged exactly as archived: one row per symmetry-independent atom (phonopy BORN form). The α record has no Born charges or ε∞ in the database ("Unknown"), and these are staged as null.

## The synthetic measurements (what is constructed, not released)

The spectra the lab serves are not measurements. They are simulated by `environment/lab/task_lab.py` from a forward model of a hidden sample, built by `benchmarks/construction/raman-polymorph-orientation/make_task.py` (with `model.py`) and stored in `private/lab_model.json`. The model has these ingredients:

1. **Main phase:** a β-SnF2 (mp-7456) single crystal. The orientation matrix is Q = Rz(48°) Ry(53°) Rz(62°), mapping the database Cartesian frame to the lab.
2. **Phonons:** the archived frequencies and eigenvectors, with the analytic dynamical matrix rebuilt as D = Σ ω_j² e_j e_jᵀ. The non-analytic (LO) term for a phonon wavevector along the face-A normal is built from the archived Born charges, expanded to all atoms with the space-group operations (acoustic sum rule asserted), and from ε∞. Every frequency is multiplied by a global scale factor s = 1.0237. The Raman tensors of the mixed modes are projections of the archived tensors; the electro-optic contribution is not available and is neglected.
3. **Line intensities:** I ∝ (1 − f)/V_c · ω_s³ (n+1)/ω · |e_s·Q R Qᵀ·e_i|² for Stokes lines, and the corresponding ω_as³ n/ω for anti-Stokes lines. Here ω_s³ is photon counting and n is the Bose factor at 295 K. The analyzer-free setting is the sum over two orthogonal analyzers.
4. **Minority phase:** randomly oriented α-SnF2 (mp-2843) inclusions, volume fraction f = 0.03. They use the orientation-averaged backscattering invariants (45a² + 4γ²)/45 (parallel) and γ²/15 (crossed), mixed as cos²/sin² of the polarizer–analyzer angle, with the same s and lines 1.3 times broader.
5. **Line shapes:** pseudo-Voigt profiles. Phonon FWHM = (1.0 + 0.0065 ω) cm⁻¹, times a fixed random factor of 0.85–1.2 per line. The instrumental profile is a Gaussian of FWHM 1.2 cm⁻¹.
6. **Calibration offset:** the laser emits at 532.06 nm, while the Raman-shift axis is computed from the nominal 532.00 nm. Every reported shift is therefore the true shift plus 2.120 cm⁻¹, on both the Stokes and the anti-Stokes side. The notch filter (±18 cm⁻¹) is centred on the true laser line.
7. **Counts:** an intensity scale of 8×10⁸ and a background of 150 + 0.1·ν counts per accumulation. Laser-power drift is Gaussian with 1.5 % rms per request. The counts follow Poisson statistics, plus an offset of 500 counts and read noise of 5 counts rms per accumulation.

The instrument facts that an experimentalist would know are stated in the catalog: nominal wavelength, window, pixel size, resolution, notch, offset, read noise, power stability and stage temperature. The following are not stated: the sample identity, the orientation, s, the minority phase and its fraction, the true laser wavelength, the line widths and the intensity scale.

Simplifications, stated in the catalog where they define conventions:

- thin-plate polarization (no birefringence or refraction);
- polarization-independent throughput;
- negligible laser heating;
- a single global frequency scale;
- no temperature dependence of the β distortion near its 339 K ferroelastic transition.

## What is hidden and graded

The hidden quantities are:

- the main phase (mp-7456) and the minority phase (mp-2843);
- s;
- the c axis Q·ĉ = (0.5344, 0.5935, 0.6018), graded modulo sign and the 180° turn about the face normal;
- two face-B predictions computed with the same model: I(t_B,t_B)/I(z,z) of the strongest mode, and the highest optical-phonon frequency for q along n_B.

None of these is served. `private/hidden_values.json` lists the graded numbers, the axis offset and the true laser wavelength for the leak scan.

## Scope

This is a released start in an unusual sense. The database records are released author data used unchanged, but the "experiment" is our own construction. The graded truths follow from our hidden sample parameters through the database tensors; they are not author-archived numbers. The construction is deterministic: the lab seeds its noise with the job id and the request.
