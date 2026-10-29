# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Removes the start-up transient of each run before computing averages and energy fluctuations, and shows that whole-run variances inflate c_v and shift the c_v maximum to lower rs, most strongly for the short (20–40 ps) 2048-atom runs at 1200–1800 K.
2. Uses the size dependence of the response (P(rs) loops, the c_v maximum, or max_k S(k)/N) to discriminate: the van der Waals-like loops (dP/drs > 0) and large c_v peaks of N = 128 and 256 at 950–1000 K disappear at N = 512 and 2048, whereas the 900 K loop survives at N = 2048. Concludes crossover above about 950 K and first order at 900 K, with the critical end point bracketed between 900 and 950 K.
3. Identifies the molecular side of the 900 K transition as a (defective) molecular crystal: max_k S(k)/N stays near 0.1 at N = 2048 while liquids give values that fall with N. States that atoms still diffuse, so diffusion alone does not prove a liquid, and that the 900 K first-order transition is melting/freezing next to the melting line rather than a liquid–liquid transition.
4. Explains why small, short simulations suggested a first-order LLT: small periodic cells stabilise defective crystals and show loops, and averages over a few ps near the Widom line fluctuate strongly enough to mimic plateaus and kinks.
5. Reports the Widom line (rs of the N = 2048 c_v maximum, about 1.448, 1.464 and 1.502 at 1100, 1200 and 1400 K) with an uncertainty, and notes that small cells place the c_v maximum differently.
6. Quantifies NQE at 1000 K by comparing crossover pressures (path-integral crossover near rs 1.488, about 40 GPa below the classical one), not pressures at fixed density, and states that the transition remains a smooth crossover with quantum nuclei.
7. Notices that the PIMD files labelled 900K are byte-identical to 1000K files (their temperature column reads 1000 K), so no 900 K path-integral result is available and the 900 K NQE shift is left undetermined.
8. States what the data do not support: the precise critical temperature, experimental transition pressures (one MACE-PBE model), and quantum corrections at temperatures without path-integral data.
