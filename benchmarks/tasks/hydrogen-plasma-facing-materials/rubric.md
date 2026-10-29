# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. States the hydrogen reservoir (½E of the relaxed H₂ molecule) and that the isolated H atom would shift every formation energy by the H₂ binding energy per H (−2.31 eV here). Distinguishes equilibrium dissolution from H₂ gas from the uptake of energetic atomic H arriving from the plasma.
2. Shows from the relaxed outputs that the archive's "inequivalent sites" are starting points, many of which relax into the same configuration: equal energies to below 1 meV and identical coordination shells. Examples are the six WC runs, the two Mo runs, three TaC runs, and graphite run 39, which stopped at the BFGS step limit inside the basin of run 17. Keeps close but genuinely different configurations separate (TiN and HfB₂, about 0.05 eV apart).
3. Reports the site competition with numbers and without textbook assumptions: octahedral over tetrahedral in hcp Re (about 0.56 eV); C-bonded H over the Ta₄ void in TaC (about 0.34 eV); the Ti₄ site over the N–H configuration in TiN (about 0.05 eV). Explains that gaps of order 0.1 eV or less are within zero-point-energy differences, which the archive cannot evaluate.
4. Checks that the runs combined in each energy difference share cutoff, pseudopotentials and supercell. Notes that the isolated-atom references use material-specific cutoffs: the W atom of the W set (30/240 Ry) and of the W-compound sets (45/360 Ry) differ by about 0.13 eV, so each material needs its own atom. Notes that the TiB₂ Ti atom is a copy of the TiN atom (60/480 Ry) while the TiB₂ slabs use 35/280 Ry, and assesses the consequence with evidence rather than ignoring it.
5. Chooses the termination with the lowest surface energy from purchased slab energies, and reports the most weakly bound atom of each element. Explains why this differs from averaging over all sites and terminations.
6. Computes cohesive energies with the bulk cell and atoms at the same cutoff, not with the 60/480 Ry hydrogen supercells. Compares SBE with cohesive energy: for the flat metal surfaces the ratio is about 1.2–1.35, while in the compounds the element- and termination-specific SBEs spread far around the cohesive energy. Concludes that cohesion is not a reliable SBE proxy.
7. Gives a screening judgement that weighs hydrogen uptake against sputtering resistance. Graphite, the carbides, the nitride and the boride dissolve H endothermically, while Re does not. The light element is the weakly bound one in WC (C about 7.8 eV against W about 11.0 eV), so preferential sputtering of C is expected.
8. States what these data cannot support:
   - zero-point and vibrational free energies, hence site order at temperature and solubility rankings;
   - supercell-size or cell-relaxation convergence of the interstitials (one fixed cell per host, pressures up to about 20 kbar);
   - magnetism (all bulk runs are non-magnetic, including Cr);
   - relaxed surface vacancies (the SBEs are unrelaxed).
