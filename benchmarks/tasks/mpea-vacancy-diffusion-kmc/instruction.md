# Vacancy diffusion in refractory V–Cr–Nb–Mo–Ta–W alloys: sluggish or not, and why?

Multi-principal-element alloys are often said to show "sluggish" diffusion. Using author-released first-principles vacancy-hop energies for bcc V–Cr–Nb–Mo–Ta–W alloys, determine how fast vacancies and atoms diffuse in three alloys relative to a rule-of-mixtures benchmark, which element diffuses fastest, and whether the behaviour is controlled by site-energy (thermodynamic) disorder, by the distribution of local migration barriers, or by neither.

## Data and calculation service

`/app/catalog.json` describes the archived hops, the alloys, the conventions of the hop records and the two purchase types:

- **`hops`** returns the archived records of up to 20 hops: each hop's occupation of the 128-site supercell and the DFT total energies of its two end states and its transition state. These are retrieved author calculations (Lee and Natarajan, 2026), not new DFT. Each purchase costs 1 unit.
- **`kmc`** runs a kinetic Monte Carlo simulation, in the lab's engine, of a vacancy-barrier model you specify, with optional canonical Monte Carlo equilibration of the starting configurations. Its cost depends on the cell, the number of trajectories and the run length (see the catalog).

```sh
python3 /app/lab_client.py evaluate '{"kind": "hops", "ids": ["H0040", "H0041"]}'
python3 /app/lab_client.py job/3
python3 /app/lab_client.py status
```

The budget is 60 units, far short of the catalog. Hop records are long: save results to files. Poll running jobs rather than resubmitting. You may also run your own simulations with the tools in your container.

## Conventions

- **Graded physical model:** one vacancy in the dilute limit, nearest-neighbour vacancy–atom exchanges only, with rate ν exp(−ΔE/k_BT), ν = 1 THz, and ΔE = E_TS − E_initial for the actual local configuration, as given by the DFT energetics of the released calculations. The alloy configuration is the equilibrium state at the temperature of interest under the same energetics. Use one lattice parameter a for all alloys and pure elements; ignore vibrational entropy and thermal expansion.
- **Hop barrier** of an entry in `hops_to_predict`: E_TS − E(end state 1), in eV. End state 1 has the vacancy on `vacancy_site_initial`.
- **Diffusivities:** vacancy D_Va = ⟨ΔR_Va²⟩/(6t). Tracer D*_i = ⟨ΔR²⟩/(6t) averaged over the atoms of element i, divided by the vacancy site fraction (per unit vacancy concentration).
- **Rule-of-mixtures benchmark:** a pure element i has D_Va,i = a²ν exp(−E_m,i/k_BT), where E_m,i is its vacancy migration barrier computed with the same DFT protocol. The enhancement factor is η = D_Va / Π_i D_Va,i^(x_i) for atomic fractions x_i. Diffusion is sluggish if η < 1.
- **Model variants for alloy B** (each at its own equilibrium configuration):
  - *Thermo-ideal:* every configuration has the same energy, so each hop has the barrier E_TS − (E_1 + E_2)/2 of its actual local environment.
  - *KRA-simple:* end-state energies as in the real alloy, but E_TS − (E_1 + E_2)/2 of every hop equals the pure-element value E_m of the hopping element.
- **Activation energy:** E_a = −k_B d ln D*/d(1/T) for the fastest-diffusing element of alloy C, from a least-squares fit of ln D* against 1/T at 1700, 2000 and 2300 K.

Unless stated otherwise, quantities refer to T = 2000 K.

## Deliverables

Save `/app/submission.json`:

```json
{"barrier_predictions_eV": {"<hop id>": 0.0},
 "log10_eta_vacancy": {"A": 0.0, "B": 0.0, "C": 0.0},
 "fastest_tracer_B": "X",
 "slowest_tracer_B": "X",
 "log10_tracer_ratio_B": 0.0,
 "activation_energy_C_eV": 0.0,
 "log10_eta_thermo_ideal_B": 0.0,
 "log10_eta_kra_simple_B": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`barrier_predictions_eV`:** a barrier for every id in `hops_to_predict`.
- **`log10_eta_vacancy`:** log10 η for the alloys A, B and C of the catalog.
- **`log10_tracer_ratio_B`:** log10 of D* of the fastest over D* of the slowest element in B.
- **`log10_eta_thermo_ideal_B`, `log10_eta_kra_simple_B`:** log10 η of alloy B for the two variants.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which records you used, and any problems found in them;
- the barrier model and how it was validated;
- how the simulations were set up, equilibrated and converged;
- how tracer diffusion differs from the hop frequencies of each element, and why;
- whether each alloy is sluggish or anti-sluggish and what controls it, with the evidence;
- the uncertainties, and which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the barrier predictions, the enhancement factors, the fastest tracer and tracer ratio of B, the activation energy and the two variants; physical reasoning is reviewed separately.
