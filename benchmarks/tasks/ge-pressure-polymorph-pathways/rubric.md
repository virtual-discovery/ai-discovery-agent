# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Identifies that most single-point records are rattled or strained training configurations (forces of order 1–3 eV/Å, shear or anisotropic stress) and excludes them from the equations of state, with a stated criterion.
2. Obtains transition pressures from equal enthalpies (E + PV from E(V) fits, or a common tangent), not from comparing energy minima. Does not take pressures from the printed stresses of intermediate relaxation steps, or from PSTRESS targets that some runs never reached.
3. Notices that the compressed ST12, BC8 and R8 frames carry residual forces, because the ISIF = 6 runs relax the cell but not the internal coordinates. Explains that their energies are therefore upper bounds, and corrects for this (e.g. with ionic relaxations using the released potential), quantifying the shift of the β-Sn/metastable crossings.
4. Assesses the released potential against the DFT data. States that it reproduces the transition pressures to a few tenths of a GPa, but not the few-meV energy differences between the metastable phases, and uses it only where that accuracy suffices.
5. States that thermodynamically CD is the stable phase below the CD/β-Sn pressure for every decompression path. ST12, BC8 and R8 become favoured over β-Sn only at lower pressures and remain metastable relative to CD, so their formation on decompression requires kinetic hindrance of the reconstruction back to diamond.
6. Explains that the near-degeneracy of ST12, BC8 and R8 (a few meV/atom, with an ST12/BC8 crossing near 1 GPa) means static enthalpies cannot decide which metastable phase forms. Rate dependence (slow decompression → ST12, fast → BC8/R8) is attributed to kinetics that these data do not resolve.
7. Discusses that PBE underbinds Ge and gives a CD → β-Sn pressure below the experimental ~10–11 GPa, so PBE pressures need this caveat before comparison with experiment.
