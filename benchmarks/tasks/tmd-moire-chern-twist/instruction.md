# Twisted MoTe2 and WSe2: how the moiré Chern bands evolve with twist angle

Experiments find opposite Chern numbers at filling ν = −1 in R-stacked twisted bilayer MoTe2 near 3.7° and twisted bilayer WSe2 near 1.2°. A large-scale DFT study (Zhang, Wang, Liu, Fan, Cao and Xiao, Nat. Commun. 2024) computed relaxed moiré superlattices of both materials down to 1.25° and released its source data. Using those archived results, establish whether and how the Chern numbers of the topmost K-valley moiré valence bands change with twist angle, and whether a fixed continuum model or the angle dependence of the moiré potential explains the change.

## Data and calculation service

`/app/catalog.json` describes everything below. Staged in `/app/archive/` (archived files, unmodified):

- `fig1/set_P`, `fig1/set_Q`: the two tMoTe2 continuum-model parameter maps of the paper's Fig. 1. One set comes from the local-stacking approximation and one was fitted to a DFT calculation of 3.89° tMoTe2; the labels P and Q are arbitrary. The catalog gives the effective mass and monolayer lattice constant that belong to each set.
- `fig2/`: archived SIESTA moiré valence bands of 3.15° tWSe2 and 3.89° tMoTe2, and relaxed 3.15° and 1.25° tWSe2 structures.
- `fig3/tWSe2_3.15`: the surface Hartree potentials of 3.15° tWSe2 with the authors' files and scripts.

Purchases return further archived results (no new DFT is run):

| Request | Cost | Returns |
|---|---|---|
| `dft_bands` | 1 | SIESTA moiré valence bands (γ–m–κ–γ) of one calculation |
| `c3` | 1 | spin and C3z eigenvalue of the six highest moiré valence states at γ, κ or κ′ of one calculation, as determined by the authors |
| `hartree_map` | 3 | tWSe2 top- and bottom-surface Hartree potential maps at one angle |
| `hartree_linecut` | 1 | the archived Hartree line cut of one tWSe2 angle |
| `hartree_fourier` | 1 | the archived Hartree Fourier-amplitude file of one tWSe2 angle |
| `structure` | 1 | the relaxed tWSe2 structure of one angle |

Calculations exist for tWSe2 at 3.15°, 1.70°, 1.54°, 1.25° and tMoTe2 at 3.89°, 3.15°, 2.14°, 1.25° (Hartree data: tWSe2 only).

```sh
python3 /app/lab_client.py evaluate '{"kind": "c3", "material": "tWSe2", "angle": "3.15", "point": "gamma"}'
python3 /app/lab_client.py evaluate '{"kind": "hartree_map", "angle": "1.70"}' > map170.json
python3 /app/lab_client.py status
```

The budget is 27 units. Results print as JSON; redirect large ones to a file. If a request returns a job id, poll `python3 /app/lab_client.py job/<id>`.

## Conventions

- **Sector and bands.** Chern numbers refer to the K-valley, spin-up moiré valence bands; band 1 is the topmost band, band 2 the next. Where the evidence fixes a Chern number only modulo 3, report the representative of smallest magnitude.
- **Continuum model** (the paper's Eq. 1), energies of electrons in meV:
  H = [[−ħ²(k − κ₊)²/2m* + Δ_b(r), Δ_T(r)], [Δ_T(r)*, −ħ²(k − κ₋)²/2m* + Δ_t(r)]],
  with κ± = (4π/3a_M)(−√3/2, ∓1/2), a_M = a₀ / (2 sin(θ/2)), and g_j = (4π/√3a_M)(cos((j−1)π/3), sin((j−1)π/3)). The maps tabulate Δ_b (`Delta-bot.dat`), Δ_t (`Delta-top.dat`) and Δ_T (`Delta-T.dat`), which have the forms Δ_{b,t}(r) = 2V Σ_{j=1,3,5} cos(g_j·r ± ψ) and Δ_T(r) = w(1 + e^{−i g₂·r} + e^{−i g₃·r}). Each grid lists values at fractional coordinates (i/(N−1), j/(N−1)) of a right-handed pair of moiré lattice vectors (a₁ × a₂ along +z); the last row and column repeat the first.
- **Chern sign.** C = (1/2π) ∫ Ω d²k with Ω = ∂_{kx}A_y − ∂_{ky}A_x and A = i⟨u|∇_k u⟩, for the band of the Hamiltonian above in these Cartesian axes; the DFT Chern numbers use the same sector and sign.
- **Hartree potential drop.** Δv_H(r) = v_H(top surface) − v_H(bottom surface) in meV, where the top layer is the one at larger z in the archived structures. Local stackings: MM (metal above metal), XM (the chalcogen of the top layer above the metal of the bottom layer) and MX (the metal of the top layer above the chalcogen of the bottom layer).
- **First-star part.** The contribution of the six shortest nonzero moiré reciprocal vectors to the Fourier series of Δv_H, evaluated at XM.

## Deliverables

Save `/app/submission.json`:

```json
{"dft_chern_K_up": {"tWSe2_3.15": [0, 0], "tWSe2_1.70": [0, 0], "tWSe2_1.54": [0, 0], "tWSe2_1.25": [0, 0],
                    "tMoTe2_3.89": [0, 0], "tMoTe2_3.15": [0, 0], "tMoTe2_2.14": [0, 0], "tMoTe2_1.25": [0, 0]},
 "continuum_chern_3p89": {"P": [0, 0], "Q": [0, 0]},
 "continuum_set_matching_dft": "P or Q",
 "continuum_P_chern_2p14": [0, 0],
 "hartree_xm_meV": {"3.15": 0.0, "1.70": 0.0, "1.54": 0.0, "1.25": 0.0},
 "hartree_first_star_xm_meV": {"1.70": 0.0, "1.25": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`dft_chern_K_up`:** (C1, C2) of bands 1 and 2 of each DFT calculation.
- **`continuum_chern_3p89`:** (C1, C2) of each continuum parameter set at θ = 3.89°, with its own m* and a₀.
- **`continuum_set_matching_dft`:** the set that reproduces the archived 3.89° tMoTe2 DFT bands.
- **`continuum_P_chern_2p14`:** (C1, C2) of set P, with all its parameters unchanged, at θ = 2.14°.
- **`hartree_xm_meV`:** Δv_H at XM of tWSe2 at each angle.
- **`hartree_first_star_xm_meV`:** the first-star part of Δv_H at XM at 1.70° and 1.25°.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how each Chern number was established and cross-checked;
- which continuum parameters describe the DFT bands, and why the two sets disagree;
- how and why the band topology changes with twist angle, at which momenta and between which bands;
- what the Hartree potential shows about the mechanism;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the Chern numbers, the continuum-set choice and the Hartree values; physical reasoning is reviewed separately.
