# Cu2WSe4: orbital versus spin Edelstein response of a gyrotropic, achiral crystal

Cu₂WSe₄ (space group I-42m, point group D2d) is a nonmagnetic semiconductor whose structure is gyrotropic but not chiral. Use the authors' released Wannier tight-binding models (OpenMX, GGA, with and without spin–orbit coupling) to answer four questions:

- Is the magnetization induced by a dc current (the Edelstein effect) carried by orbital or by spin moments, in the valence and in the conduction band?
- Is spin–orbit coupling essential for the orbital response of the valence band?
- Which components of the Edelstein tensor and of the second-order (nonlinear Hall) conductivity does the crystal allow, and how are they related?
- What is the trace of the gyrotropic tensor, which has been proposed as a chirality indicator?

## Data and calculation service

The released files of Zenodo record 20131063 are staged in `/app/archive`:

- the crystal structure (`Fig1a.cif`);
- the DFT band structures with and without spin–orbit coupling (`Fig1b_soc.d.gz`, `Fig1b_wosoc.d.gz`, gzip-compressed);
- the projected densities of states (`Fig1c_*.d`);
- the header lines of the two Wannier Hamiltonians (`Cu2WSe4_soc.HWR.header`, `Cu2WSe4_wosoc.HWR.header`).

The archive describes none of its columns. The archived response curves and the Hamiltonians themselves are not staged. `/app/catalog.json` lists the files, the request formats, the costs and the conventions of the service.

The service evaluates the two archived Wannier Hamiltonians, selected by `"model"`: `"soc"` (`Cu2WSe4_soc.HWR`) or `"nosoc"` (`Cu2WSe4_wosoc.HWR`).

- **`bands`** (1 unit): band energies at up to 1000 k-points.
- **`band_properties`** (1 unit): at up to 200 k-points, the band energies, band velocities and Berry curvature, plus three magnetic moments per band:
  - the spin moment;
  - the Berry-phase orbital moment;
  - the intra-atomic orbital moment.
- **`fermi_surface`**: Fermi-surface sums on a uniform N × N × N/2 mesh, at up to 400 chemical potentials and one temperature. It returns the electron count and, for each moment, the tensor Σ (−∂f/∂ε) m_i ∂ε/∂k_j, plus the Berry-curvature dipole tensor Σ (−∂f/∂ε) (∂ε/∂k_i) Ω_j. The cost depends on the model and the mesh:

  | N | 16 | 24 | 32 | 48 | 64 | 96 | 128 |
  |---|---|---|---|---|---|---|---|
  | `soc` (units) | 1 | 1 | 1 | 2 | 3 | 6 | 10 |
  | `nosoc` (units) | 1 | 1 | 1 | 1 | 2 | 3 | 5 |

```sh
python3 /app/lab_client.py evaluate '{"kind": "bands", "model": "soc", "k_frac": [[0, 0, 0], [0.5, 0.5, 0]]}'
python3 /app/lab_client.py evaluate '{"kind": "band_properties", "model": "nosoc", "k_frac": [[0.1, 0.2, 0.05]]}'
python3 /app/lab_client.py evaluate '{"kind": "fermi_surface", "model": "soc", "mesh": 48, "mu_eV": [-0.7, 0.75], "kT_eV": 0.03}'
python3 /app/lab_client.py status
```

The budget is 30 units. A request body may not exceed 64 kB. A request that takes longer than 20 s returns a job id; poll it with `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **Axes and k-points:** Cartesian x, y and z lie along a, b and c of the conventional cell in `Fig1a.cif`, which is also the cell of both Wannier models. `k_frac` are coordinates in the reciprocal basis of that cell; `k_cart_invA` are Cartesian components in Å⁻¹, including 2π.
- **Energies:** eV, relative to the Fermi level printed in the header of the model used; chemical potentials in requests use the same scale. The valence-band maximum (VBM) and the conduction-band minimum (CBM) are the extrema of the model's bands over the whole Brillouin zone.
- **Moments:**
  - The spin moment of a state is −μ_B⟨σ⟩ (g = 2).
  - The orbital moment of a Bloch state is that of its wave packet, m = (q/2)⟨(r − r_c) × v⟩ with q = −|e|. Here r is the position operator of the tight-binding model, diagonal in the Wannier basis, with each Wannier function at the position of its atom. The orbital magnetization is the response of this moment alone; do not add any other orbital contribution.
- **Edelstein effect:** M_i = Σ_j α_ij E_j, with
  - α_ij = (qτ/V) Σ_{n,k} (−∂f/∂ε) m_i v_j, where q = −|e|, v = ∇_k ε/ħ, V is the crystal volume and f is the Fermi–Dirac distribution;
  - scattering rate γ = ħ/(2τ) = 30 meV and k_BT = 30 meV;
  - the field is E = (+10⁴ V/m, 0, 0);
  - the orbital (spin) magnetization is the part carried by the orbital (spin) moments.
- **Doping:** rigid bands of each model at k_BT = 30 meV; the hole density p is the number of empty valence states of the real crystal per cm³.
- **Nonlinear Hall effect:** dc current j_i = Σ_{j,l} σ_ijl E_j E_l, with σ_ijl = σ_ilj.

## Deliverables

Save `/app/submission.json`:

```json
{"orbital_to_spin_ratio": {"VBM-0.10": 0.0, "CBM+0.10": 0.0},
 "orbital_magnetization_A_per_m": 0.0,
 "soc_to_nosoc_orbital_ratio": 0.0,
 "edelstein_allowed_components": ["xx", "xy"],
 "alpha_yy_over_alpha_xx": 0.0,
 "chirality_indicator": 0.0,
 "nonlinear_hall_allowed_components": ["xyz", "xxx"],
 "sigma_zxy_over_sigma_xyz": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`orbital_to_spin_ratio`:** α^orb_xx / α^spin_xx of the spin–orbit-coupled model at μ = VBM − 0.10 eV and at μ = CBM + 0.10 eV.
- **`orbital_magnetization_A_per_m`:** the orbital magnetization M_x (A/m) of the spin–orbit-coupled model at μ = VBM − 0.10 eV.
- **`soc_to_nosoc_orbital_ratio`:** α^orb_xx of the spin–orbit-coupled model divided by α^orb_xx of the model without spin–orbit coupling, both evaluated for the crystal at a hole density p = 5.0 × 10²⁰ cm⁻³.
- **`edelstein_allowed_components`:** every component "ij" of α_ij (i is the magnetization, j the field) that the crystal symmetry allows to be nonzero in these axes.
- **`alpha_yy_over_alpha_xx`:** α_yy / α_xx of the crystal.
- **`chirality_indicator`:** the trace of α divided by α_xx, (α_xx + α_yy + α_zz) / α_xx.
- **`nonlinear_hall_allowed_components`:** every component "ijl" of σ_ijl that the crystal symmetry allows to be nonzero in these axes.
- **`sigma_zxy_over_sigma_xyz`:** σ_zxy / σ_xyz for the Berry-curvature-dipole contribution of the crystal.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the point-group symmetry of the structure in these axes, how you established it, and which tensor components it allows;
- how the orbital moment was obtained, and why that choice;
- where the band edges lie and how the chemical potentials and hole densities were obtained;
- how the Fermi-surface sums were converged and checked;
- whether orbital or spin moments dominate in each band, and the role of spin–orbit coupling;
- what dominates the nonlinear Hall response, and what the chirality indicator shows;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks:

- the two orbital-to-spin ratios;
- the orbital magnetization;
- the SOC-to-no-SOC ratio;
- the allowed Edelstein components, α_yy/α_xx and the chirality indicator;
- the allowed nonlinear Hall components and σ_zxy/σ_xyz.

Physical reasoning is reviewed separately.
