# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Writes the lattice-mediated formula α = (μ₀/Ω) Σ_m (∂m/∂u_m)(Z*·e_m)/k_m with eigenvectors and eigenvalues of the force-constant matrix, explains why mass-weighted phonon eigenvectors cannot be paired with the frozen force-constant-mode responses, and states how the units of the archived eigenvalues were established (for example against the all-H phonopy frequencies).
2. Establishes the domain of each dataset from the converged Cr moments and the nearest-neighbour (face-sharing) Cr pair, notes that the Elk input fields point opposite to its converged moments, and concludes that the VASP and Elk induced-moment runs model opposite domains, so their raw signs must not be averaged.
3. Constructs the Born charges with the correct VASP dipole sign (electron units) and the stated index order, checks linearity and the acoustic sum rule (one-sided differences against Disp_0_0 violate it in the z row because that reference run carries a dipole offset), and reports that the Disp_-0_005 runs are inconsistent (their OUTCAR positions and ionic dipoles belong to +0.05 Å or +0.005 Å runs, not to −0.005 Å).
4. Uses the net cell moment (interstitial included) rather than sphere-projected moments, fits symmetric amplitudes in the linear regime, and identifies the anomalous Elk Eu_3 +0.03 Å run (early "Convergence targets achieved", off the linear trend).
5. Gives α⊥ in SI (ps/m) with the rhombohedral cell volume, a sign for each domain, and a magnitude check: about 1 ps/m, and dimensionless cα of order 10⁻⁴.
6. Identifies the stiffer Eu_3/Eu_4 force-constant modes, not the softest mode, as the carriers of the response, with their fractions, and relates this to their large polarity and spin-canting response.
7. Argues the annealed domain from the temperature dependence: α∥ shares the 0 K sign of α⊥ but reverses on warming, and annealing in parallel E and H favours α∥ > 0 near T_N, which gives the out-pointing domain. Distinguishes moments from spins.
8. States what these data do not support: clamped-ion, orbital and finite-temperature contributions, the size of α∥, and the experimental magnitude of α⊥.
