# ZrSiX nodal-line semimetals: what shapes the magnetoresistance?

ZrSiS, ZrSiSe and ZrSiTe are tetragonal (P4/nmm) nodal-line semimetals with a large, non-saturating magnetoresistance. In ZrSiS and ZrSiSe the transverse magnetoresistance has a "butterfly" angular pattern: its maximum lies between the c axis and the plane, not where the Lorentz force is largest. It has been attributed to the topology of the nodal line, to a field-induced topological transition, to open orbits and to effective-mass anisotropy. Using the authors' released Wannier models (GGA, no spin–orbit coupling), determine whether semiclassical Boltzmann transport on the actual Fermi surfaces accounts for it, and quantify the carrier balance, the angular pattern, the field dependence and the compound trend.

## Data and calculation service

The released files of Materials Cloud record tyan3-q3j67 are staged in `/app/archive` under their own folder names: the VASP and Wannier90 inputs (`input_vasp+wannier90/<compound>/`; the POTCARs are not staged), the WannierTools inputs (`input_MR/<compound>/input.dat`), the header of each Wannier Hamiltonian (`input_MR/<compound>/*.header`) and the record README (`README.txt`). The Hamiltonians themselves are not staged. `/app/catalog.json` lists the request formats, the menus, the costs and the cells used by the service.

The service evaluates the three archived Hamiltonians (`"compound"`: `"ZrSiS"`, `"ZrSiSe"` or `"ZrSiTe"`; 24 orbitals, no spin index):

- **`bands`** (1 unit): band energies at up to 500 k-points, optionally with band velocities.
- **`fermi_slice`** (1 unit): energies of bands 10–15 on a grid of up to 121 × 121 points in a plane of k-space.
- **`electron_count`** (1 unit): band occupations on a uniform N × N × N/2 mesh (N = 40, 80, 120 or 160) at up to 200 chemical potentials and a temperature of your choice.
- **`orbit`** (1 unit): semiclassical k-space orbits (up to 8 start points) on one band for a field direction, with closure, period and cyclotron mass.
- **`magnetotransport`**: the Chambers (Boltzmann) conductivity and resistivity tensors for one field direction from the menu, the bands you select, up to 4 chemical potentials and a temperature from 10 to 40 K, at 15 values of Bτ from 0 to 20 T ps. Mesh N = 40 or 80 costs 1 unit, N = 160 costs 3.

```sh
python3 /app/lab_client.py evaluate '{"kind": "bands", "compound": "ZrSiSe", "k_frac": [[0, 0, 0], [0.5, 0.5, 0]], "velocities": true}'
python3 /app/lab_client.py evaluate '{"kind": "fermi_slice", "compound": "ZrSiSe", "origin": [-0.5, -0.5, 0], "u": [1, 0, 0], "v": [0, 1, 0], "n": [61, 61]}'
python3 /app/lab_client.py evaluate '{"kind": "electron_count", "compound": "ZrSiTe", "mesh": 80, "mu_eV": [8.80, 8.90], "temperature_K": 30}'
python3 /app/lab_client.py evaluate '{"kind": "orbit", "compound": "ZrSiSe", "band": 12, "k_frac": [[0.2, 0.2, 0.0]], "field_direction": [0, 0, 1], "max_length_invA": 20}'
python3 /app/lab_client.py evaluate '{"kind": "magnetotransport", "compound": "ZrSiSe", "field": {"plane": "bc", "angle_deg": 45}, "mesh": 80, "mu_eV": [8.92], "temperature_K": 30, "bands": [12, 13]}'
python3 /app/lab_client.py status
```

The budget is 50 units. A request body may not exceed 64 kB. A request that takes longer than 20 s returns a job id; poll it with `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **Axes:** x, y and z lie along a, b and c of each compound's tetragonal cell (listed in the catalog). `k_frac` are coordinates in the reciprocal basis of that cell.
- **Energies:** eV on the absolute scale of each Wannier Hamiltonian. Bands are numbered 1–24 from the bottom.
- **Crystal:** the undoped, stoichiometric compound, with rigid bands, at T = 30 K (Fermi–Dirac). Its chemical potential `mu` is the one at which the electron count equals that of the neutral crystal. The carrier density n is the density of electrons in the conduction bands, the bands numbered above N_v/2, where N_v is the number of valence electrons per cell that the 24-orbital model describes; it is given per cm³ of the real crystal, counting both spin directions.
- **Transport:** Boltzmann equation with one constant relaxation time τ for all bands; σ_ij(B) = (2e²/V) Σ_n ∫ d³k/(2π)³ (−∂f/∂ε) v_i v̄_j, v̄ the velocity averaged over the past orbit with weight e^(−t/τ), ħ dk/dt = −e v × B. The resistivity is ρ = σ⁻¹ summed over all bands that cross the chemical potential. Fields are given as Bτ (tesla × picosecond).
- **Geometry:** current along x. "Rotation in the bc plane" means B = B (0, sin θ, cos θ), θ measured from c towards b, so that the field is always perpendicular to the current. "B ∥ b" is θ = 90°.
- **Magnetoresistance:** MR(Bτ) = ρ_xx(Bτ)/ρ_xx(0) − 1.

## Deliverables

Save `/app/submission.json`:

```json
{"neutral_mu_eV": {"ZrSiS": 0.0, "ZrSiSe": 0.0, "ZrSiTe": 0.0},
 "carrier_density_cm3": {"ZrSiSe": 0.0, "ZrSiTe": 0.0},
 "butterfly_peak_angle_deg": 0.0,
 "mr_b_over_c": 0.0,
 "mr_exponent_b": 0.0,
 "zrsite_mr_c_over_b": 0.0,
 "mr_ratio_se_over_te": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`neutral_mu_eV`:** the chemical potential of each undoped crystal at 30 K.
- **`carrier_density_cm3`:** n of undoped ZrSiSe and ZrSiTe at 30 K.
- **`butterfly_peak_angle_deg`:** for undoped ZrSiSe at 30 K and Bτ = 5 T ps, the angle θ in [0°, 90°] of the bc-plane rotation at which ρ_xx is largest.
- **`mr_b_over_c`:** for the same crystal and Bτ, MR(B ∥ b) / MR(B ∥ c).
- **`mr_exponent_b`:** for the same crystal with B ∥ b, the exponent ln[MR(20 T ps)/MR(10 T ps)] / ln 2.
- **`zrsite_mr_c_over_b`:** for undoped ZrSiTe at 30 K and Bτ = 5 T ps, MR(B ∥ c) / MR(B ∥ b).
- **`mr_ratio_se_over_te`:** two samples measured at B = 5 T along c and 30 K: undoped ZrSiSe with τ = 1.0 ps and undoped ZrSiTe with τ = 0.4 ps. Give MR(ZrSiSe) / MR(ZrSiTe).
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the chemical potentials and carrier densities were obtained, and what the archived inputs imply;
- the Fermi-surface pockets and orbits of each compound, and where open orbits occur;
- how the transport results were converged and checked;
- what produces the butterfly pattern in ZrSiSe and the different pattern in ZrSiTe, and whether the magnetoresistance saturates;
- what the Hall resistivity shows and how robust it is;
- whether nodal-line topology is needed to explain these results, and which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the chemical potentials, the carrier densities, the butterfly angle, the two anisotropy ratios, the field exponent and the sample comparison; physical reasoning is reviewed separately.
