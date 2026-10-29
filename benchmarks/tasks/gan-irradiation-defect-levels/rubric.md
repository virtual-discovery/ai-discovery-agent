# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Identifies from the run metadata that every defect cell and the `perfect` host were computed at the single k-point (½,½,½) with the Ga PAW set (3d in the core), and uses `perfect` (not the Γ-only cell) as the host energy. Explains that the eigenvalues at that k-point are folded states, not the band edges: the supercell gap there is about 6.5 eV. Takes the VBM and CBM from the Γ-point run with the same PAW set (`perfect/kpt_gam`, gap ≈ 3.57 eV). Rejects the perfect-cell Fermi energy (Gaussian smearing places it inside the gap) and the Ga_d runs (`WAV_ref`, `kpt_gam/d10`: 864 electrons, a different eigenvalue scale) as references.
2. Applies a finite-size correction suited to relaxed charged defects in an anisotropic crystal: the point-charge energy with the static (ε0) tensor plus a potential alignment from the core potentials of atoms far from the defect, with the sign convention of the VASP potentials handled. Quantifies the corrections (point-charge term ≈ 0.19 eV × q², i.e. ≈ 1.7 eV for |q| = 3 with ε0 ≈ 10.1/11.1; alignment terms up to a few tenths of an eV). Explains that ε∞ applies to vertical (optical) transitions, not to relaxed geometries.
3. Takes the lowest-energy structure for each charge state and reports the site switching:
   - Ga_i: the i2 structure is lower by about 1 eV for q = +1 and +2, and both runs relax to the same structure for q ≤ 0;
   - N_i: the split-interstitial (i2) structure wins for −2 ≤ q ≤ 0 and for +2 and +3; i1 wins for q = −3 (by more than 2 eV) and, marginally, for q = +1.
4. Identifies the negative-U behaviour of V_N: a single (3+/+) level near 0.66 eV, with the 2+ state never stable (U ≈ −0.07 eV). Notes that for Ga_i the (3+/2+) and (2+/+) levels practically coincide (U ≈ 0), so these data cannot decide whether Ga_i^2+ has a stability window. Notes that ε∞ or an uncorrected treatment changes these conclusions.
5. Checks whether the extra carriers are bound to the defect, using the eigenvalues, occupations, site spin moments and alignment terms. Flags the candidates for delocalized (perturbed-host) states, such as Ga_N +5/+6, V_Ga–V_N +4 and N_i(i2) −2/−3: their spin spreads over many atoms, their states sit in band-like positions, or their alignment terms are anomalously large (up to ~1 eV). Notes that they create no levels in the gap. Notes that sampling at the single special k-point penalizes delocalized carriers: the conduction band there lies about 2.3 eV above the CBM.
6. Interprets the ladders for irradiated material:
   - V_Ga is a triple acceptor with deep levels in the upper half of the gap and compensates n-type GaN;
   - N_i is amphoteric with deep levels;
   - Ga_i stays a donor up to near the CBM;
   - C_N (0/−) is a deep acceptor near 0.9 eV;
   - the V_Ga–V_N divacancy is strongly bound (≈ 3 eV at E_F = E_CBM − 0.5 eV), so vacancy pairs survive;
   - N Frenkel pairs cost less than Ga Frenkel pairs at that E_F.
   States that chemical potentials cancel in these pair energies.
7. States what the data do not support:
   - finite-size convergence (one supercell size, one special k-point);
   - spin states (NUPDOWN fixed to 0 or 1; several even-electron states are broken-symmetry singlets);
   - optical or vertical levels and capture barriers (no configuration-coordinate data);
   - absolute formation energies, concentrations or a self-consistent Fermi level (no elemental references);
   - the sensitivity to the mixing fraction 0.305 and to the PBE dielectric tensors.
