# Construction of mpea-vacancy-diffusion-kmc

Released-start task built on the vacancy-hop dataset of Lee and Natarajan (2026), Materials Cloud record bhv3r-a8q53 (CC BY 4.0).

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` through `./run_in_image.sh` (cpuset from `$CPUSET`, default `4,5,20,21`; it compiles `kmc` from the lab's `kmc.c` first).

1. `python3 fetch.py` downloads `all_data.json` (143 MB, md5 15df9136ed668320057131c8c59cad2a) and `files_description.md` into `~/.cache/ai-theorist/mpea-vacancy-diffusion-kmc/`.
2. `./run_in_image.sh python3 -B extract.py` writes `pool.json` (1,945 hops mapped onto the ideal 128-site cell; all consistency assertions pass).
3. `./run_in_image.sh python3 -B reference.py` writes `reference.json`: the reference barrier model (trained on all 1,799 clean hops) and the hidden diffusion targets (600 trajectories x 50,000 hops, 10x10x10 cells, 150 equilibration sweeps; ~10 min on 4 cores).
4. `./run_in_image.sh python3 -B calibrate.py` writes `calibration.json` and `solution/ablation_values.json` (robustness variants and wrong paths; ~30 min).
5. `./run_in_image.sh python3 -B make_task.py` writes the catalog, lab data, private targets and public budgets; `python3 -B write_design.py` writes `design.json` (decision points and the calibration table).
6. `.venv/bin/python -m benchmarks.tools.check_task mpea-vacancy-diffusion-kmc --cpuset 4,5,20,21` writes `check.json`.

`kmclib.py` holds the geometry, training and kMC driver; `common.py` the alloy definitions, hidden-hop selection (seed 2026) and model variants.

## Design as built

- **Lab:** `hops` (20 archived records per unit) and `kmc` (the lab engine `environment/lab/kmc.c`, compiled at image build, runs the agent's KRA + pair-Hamiltonian model; cost by work, 1-8 units). Budget 60. 1,845 hops are purchasable (93 units); 100 clean alloy hops are hidden and graded by barrier RMSE.
- **Graded:** barrier RMSE on the hidden hops (<= 0.23 eV); log10 eta at 2000 K for A = CrNbMoTaW, B = VCrNbW, C = Nb50Ta25W25 (+-0.3); fastest tracer of B (exact) and log10 tracer ratio (+-0.3); activation energy of the fastest tracer (Nb) of C over 1700-2300 K (+-0.2 eV); log10 eta of the thermo-ideal and KRA-simple variants of B (+-0.3).
- **Reference values:** log10 eta A -1.574, B -1.199, C -0.531; ratio 0.958 (Nb/Cr); E_a 1.107 eV; thermo-ideal B -0.505; KRA-simple B -0.103. The calibration table is in `design.json` (`calibration`).
- **Oracle:** 800 hops (all 35 pure/dilute + 765 random alloy hops; 40 units) and 7 kMC runs (14 units): 54 of 60 units, 47 evaluations, 35 s in the gate container; Harbor oracle reward 1.0 (1 min 19 s). Its deviations from the reference are at most 0.12 dex and 0.04 eV; barrier RMSE 0.164 eV.

## Deviations from the design record

- **Truth.** The design asked for values cross-checked against the paper. The authors' end-state energy model (a separate DFT cluster expansion from Lee et al.) and their neural-network eLCE are not archived. Our reference uses the released hop energies for both the KRA and the end states (pair Hamiltonian on end-state total energies). With these energetics the site-energy disorder is large (median |E2 - E1| = 0.43 eV), the thermo-ideal model raises eta by 0.7 dex in B, and the paper's binary classifications (Mo-W and Cr-W anti-sluggish) are not reproduced (our model gives them sluggish). The paper's mechanism (local barriers dominate) is reproduced for B (KRA-simple raises eta by 1.1 dex). No paper headline composition is graded; all graded alloys are sluggish, so the sluggish/anti-sluggish classification target was replaced by graded log10 eta values.
- **Dropped forks.** Short-range-order equilibration and kMC cell size/run length are not graded. Where they matter (Cr-Nb-Mo-Ta: random vs equilibrated configurations shift log10 eta by 0.4-0.6; L = 4 cells by 0.4), the alloy sits near a miscibility gap at 2000 K (it phase-separates at 1600 K) and its eta is too sensitive to the purchased subset (log10 eta from -2.0 to -3.4 across random 800-hop subsets, and -3.0 for a target-relevant purchase plan). They remain rubric items. In the graded alloys they change log10 eta by < 0.15 (see `calibration.information_only`).
- **Pathological entries** (146 hops with the saddle below an end state, mostly Nb/Ta hops in V/Cr-rich cells) are served unchanged; keeping or dropping them changes no graded value beyond the robustness spread, so this fork is a robustness variant and a rubric item, not an ablation.
- **"Nearest available environment" requests** were replaced by direct purchase of catalog hops (ids with composition and hopping element); the catalog gives no energies.
- 922 records store end-state structures with a 14.29 Å cell (transition state 11.4-12.5 Å). Occupations come from fractional coordinates, so this does not affect the served data; both edges are served as archived.

## Decision points (6) and ablations

dp1 end-state energies / detailed balance (thermo-ideal as real; direct forward-barrier fit), dp2 species-mean barriers, dp3 tracer from jump counts, dp4 rule-of-mixtures pure barriers from the linear model's pure limit, dp5 activation energy from a barrier average (ROM barrier; mean accepted barrier), dp6 mechanism attributed without the variants. Tightest margins: barrier RMSE robust max 0.202 vs tolerance 0.23 vs species-mean ablation 0.268; log10 eta A robust max deviation 0.22 vs species-mean ablation 0.39 and linear-pure-limit ablation 0.36; tracer ratio robust max 0.17 vs species-mean 0.37; E_a robust max 0.11 eV vs ROM barrier 0.33 eV; KRA-simple robust max 0.25 vs tolerance 0.3.

## Risks

- Truth is a builders' surrogate. Fits to fewer than ~800 hops widen the spread (a 600-hop random subset gave E_a 0.29 eV off), so an agent that buys few hops can fail E_a or eta with an otherwise sound workflow.
- Fitting the pair Hamiltonian to swap energies only (E2 - E1, no total energies) leaves the mixing energetics loosely constrained; it shifted log10 eta of the former Cr-Nb-Mo-Ta alloy by 0.32 and was not calibrated for the final alloys.
- Agents that recall the paper (kinetics dominate, Mo-W/Cr-W anti-sluggish) and ignore the end-state energies fail; this is intended, but the graded thermodynamic contribution is a property of the released energies, not of the authors' energy model.
- Statistical noise: an agent's log10 eta from 40 x 20,000-hop runs carries ~0.07 dex of noise.
