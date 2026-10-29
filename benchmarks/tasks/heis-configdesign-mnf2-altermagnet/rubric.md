# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. States the exchange convention (unit vectors, ordered pairs, per bond) and shows, from the raw VASP total energies, that the archive energies documented as "meV/atom" are total-energy differences of the 6-atom cell (2 Mn). Normalizes J per Mn accordingly.
2. Derives the spiral geometry from the released inputs and moments: the in-cell angle of Mn2 is −2π q·τ(Mn1), so every Mn carries the azimuth 2π q·r (a uniform spiral). Identifies q = 0 as the ferromagnet and the Néel state at integer q with h + k + l odd as the ground state (the energy minimum). Notes that E(q) has the periodicity of the body-centred reciprocal lattice.
3. Identifies the [214] spirals near t = 1 as basis-truncated, from VASP's "choose ENMAX > X" exceeding ENCUT, or from the ~1.5 meV mismatch between symmetry-equivalent Néel energies. Excludes them, and treats the mismatch as a numerical artefact rather than as non-Heisenberg or altermagnetic physics.
4. Shows aliasing with a rank, conditioning or cross-validation argument:
   - single or high-symmetry directions ([001] alone; [001] + [121] beyond 8 shells);
   - commensurate q grids.

   Chooses the shell count by validation on configurations not used in the fit, and finds that shells beyond J1–J3 are required, including the 7.0 Å and 7.8 Å same-sublattice shells.
5. Explains that a uniform spiral gives both sublattices the same q, so only the sublattice sum J(Mn1, d) + J(Mn2, d) of each split shell enters. The altermagnetic difference, and with it the magnon chirality splitting, is therefore undetermined, which is different from zero. Names the first split shell (a⟨110⟩, 7.0 Å) by a site-symmetry argument.
6. Proposes a calculation that would determine the difference: configurations in which the two sublattices see different spin patterns along [110] (for example, collinear supercells with one sublattice modulated along [110], or spirals with sublattice-dependent q).
7. Obtains T_N from classical Monte Carlo of the fitted model, with a finite-size analysis (for example, a Binder-cumulant crossing) and an uncertainty. Contrasts it with the mean-field overestimate. Notes that T_N is insensitive, at second order, to the undetermined altermagnetic difference. Does not present agreement or disagreement with the experimental T_N as validation (PBE, classical spins).
