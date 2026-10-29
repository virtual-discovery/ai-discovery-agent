# MnTe: anomalous Hall conductivity of a compensated altermagnet under hole doping

Hexagonal MnTe (NiAs type, P6₃/mmc) orders as a compensated collinear antiferromagnet whose two Mn sublattices are related by a screw rotation, not by a translation or by inversion. Experiments find an anomalous Hall effect together with a very weak net magnetization. Using the authors' released HSE06 + spin–orbit calculation and its 32-band spinor Wannier model for in-plane Néel order, determine the intrinsic anomalous Hall conductivity of hole-doped MnTe, which Hall components the magnetic order allows, and what the released model can and cannot predict.

## Data and calculation service

The archive (Zenodo record 18863304) is staged in `/app/archive` under its own folder names: the VASP run (`1_vasp`), the Wannier90 construction (`2_wannier90`), the WannSymm symmetry analysis (`3_wannsymm`), the WannierTools input and log (`4_wtools`), and the data of the paper's Berry-curvature maps (`fig6_script`). Files are unchanged except `4_wtools/WT.out`, from which the k-point progress lines were removed. The authors' Hall-conductivity output and the Wannier Hamiltonian itself are not staged. `/app/catalog.json` lists the staged files, the request formats, the costs and the conventions of the service.

The service evaluates the archived Wannier Hamiltonian (`2_wannier90/wannier90_hr.dat`) with the Néel vector of the archived state (`"neel": "+y"`) or with all moments reversed (`"neel": "-y"`):

- **`bands`** (1 unit): band energies at up to 1000 k-points, or band energies and spin expectation values at up to 500.
- **`berry_curvature`** (1 unit): band energies and band-resolved Berry curvature at up to 500 k-points.
- **`ahc`**: the intrinsic anomalous Hall conductivity (all three components) and the electron count on a uniform N × N × N/2 mesh, at up to 500 chemical potentials and a Fermi smearing of your choice. Meshes N = 48, 72 and 100 cost 1 unit, N = 150 costs 2, N = 200 costs 3, N = 300 costs 5 and N = 400 costs 8.
- **`electron_count`** (1 unit, any mesh): the electron count alone.

```sh
python3 /app/lab_client.py evaluate '{"kind": "bands", "k_frac": [[0, 0, 0], [0, 0, 0.5]], "spin": true}'
python3 /app/lab_client.py evaluate '{"kind": "berry_curvature", "k_frac": [[0.1, 0, 0.4]]}'
python3 /app/lab_client.py evaluate '{"kind": "ahc", "mesh": 100, "mu_eV": [5.6, 5.8], "kT_eV": 0.001}'
python3 /app/lab_client.py evaluate '{"kind": "electron_count", "mesh": 200, "mu_eV": [5.8], "kT_eV": 0.001}'
python3 /app/lab_client.py status
```

The budget is 32 units. A request body may not exceed 64 kB. A request that takes longer than 20 s returns a job id; poll `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **Axes and k-points:** Cartesian axes of the POSCAR cell (a₁ along x, c along z); `k_frac` are coordinates in the reciprocal basis of that cell, `k_cart_invA` Cartesian components in Å⁻¹ including 2π. "Néel vector along y" is the archived state (moments ±y, perpendicular to a₁); "along x" means moments along a₁.
- **Energies:** absolute, in eV, on the scale of the Wannier Hamiltonian (the VASP scale). Wannier bands are numbered 1–32 from the bottom; bands 1–22 are filled in the undoped crystal.
- **VBM:** the maximum of Wannier band 22 over the whole Brillouin zone.
- **Hole density p:** holes per cm³ in the rigid bands of the Wannier model, relative to the 22 filled bands.
- **Hall conductivity:** intrinsic (Berry-curvature) conductivity in S/cm, in the zero-temperature limit, with the sign convention of the archived WannierTools output: σ_ab = (e²/ħ) ∫ d³k/(2π)³ Σₙ fₙ Ω_c,n, (a, b, c) cyclic, Ω = ∇_k × i⟨u|∇_k u⟩. The Hall vector is (σ_yz, σ_zx, σ_xy).

## Deliverables

Save `/app/submission.json`:

```json
{"sigma_xy_S_per_cm": {"2.0e21": 0.0, "6.0e21": 0.0},
 "sigma_xy_sign_changes": 0,
 "sign_change_hole_density_cm3": 0.0,
 "hole_density_0p30_below_vbm_cm3": 0.0,
 "forbidden_components_max_abs_S_per_cm": {"yz": 0.0, "zx": 0.0},
 "allowed_hall_components": {"neel_y": ["xy", "yz", "zx"], "neel_x": ["xy", "yz", "zx"]},
 "electron_doped_sigma_xy_S_per_cm": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`sigma_xy_S_per_cm`:** σ_xy of MnTe with the Néel vector along y at the hole densities given by the keys (cm⁻³).
- **`sigma_xy_sign_changes`:** the number of sign changes of σ_xy(μ) evaluated every 5 meV from μ = VBM − 1.000 eV to μ = VBM − 0.250 eV (151 chemical potentials).
- **`sign_change_hole_density_cm3`:** the hole density at which σ_xy changes sign in that range, with the crossing interpolated linearly between the grid points (if it changes sign more than once, the change at the highest hole density).
- **`hole_density_0p30_below_vbm_cm3`:** the hole density when μ lies 0.30 eV below the VBM.
- **`forbidden_components_max_abs_S_per_cm`:** the largest |σ_yz| and |σ_zx| of MnTe with the Néel vector along y for μ between 1.00 eV below the VBM and the VBM.
- **`allowed_hall_components`:** the Hall components ("xy", "yz", "zx") that the magnetic symmetry of MnTe allows to be nonzero for the Néel vector along y and along x (an empty list if none).
- **`electron_doped_sigma_xy_S_per_cm`:** σ_xy of electron-doped MnTe at an electron density of 1 × 10²⁰ cm⁻³ in the conduction band, or `null` if the released data do not determine it.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the magnetic symmetry of the archived state and which Hall components it allows, for the Néel vector along y and along x, and how you established it;
- how σ_xy depends on hole doping and where its Berry curvature comes from;
- how the hole densities and the energy reference were obtained, and how the results were converged and checked;
- how the Hall signal behaves under reversal of the Néel vector, and why a Hall signal can be observed in a multidomain sample of a compensated magnet;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks σ_xy at the two hole densities, the sign changes and the hole density of the sign change, the hole density 0.30 eV below the VBM, the forbidden components, the allowed components and the electron-doped value; physical reasoning is reviewed separately.
