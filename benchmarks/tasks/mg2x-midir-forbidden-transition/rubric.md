# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Assigns the absorption to the vertical transition between the two lowest conduction states at X (X1 → X3 in Mg₂Si and Mg₂Ge, X3 → X1 in Mg₂Sn). States that it is dipole-forbidden at X, because both states are even under inversion and the archived PBE matrix element vanishes at X and along the valley axis, and that it becomes allowed away from X with |p|² ≈ β²q⊥² at leading order. Argues against the alternatives with evidence: the measured features are peaks rather than a free-carrier tail; the Mg₂Si band lies below its interband edge; and the band energies match the QSGW/G0W0 splittings, not the LDA ones.
2. Identifies X1 and X3 from their dispersions, not from their energy order: X1 has a longitudinal mass of about 0.54 and X3 of about 1.03, and their transverse masses differ too. Reports that the order is inverted in Mg₂Sn at all three levels (X3 − X1 ≈ −0.36, −0.165 and −0.131 eV). Notes that the LDA places X3 of Mg₂Sn below the valence-band maximum, so that the G0W0 of Mg₂Sn starts from a semimetal.
3. Compares the levels of theory:
   - the LDA underestimates the splitting by 0.25–0.3 eV (0.16 and 0.33 eV against the measured 0.40 and 0.56 eV), and a scissor correction cannot change a splitting between two conduction states;
   - the all-electron G0W0 and QSGW agree within about 0.02–0.03 eV for Mg₂Si and Mg₂Ge;
   - QSGW overestimates the Mg₂Si band (0.456 against about 0.40 eV) but nearly matches Mg₂Ge (0.576 against about 0.56 eV);
   - the pseudopotential G0W0 gives 0.21–0.39 eV for Mg₂Si, depending on the charge density of the plasmon-pole model, and is not a reliable reference.
4. Explains how the TOTE files were read. For the one-shot G0W0, the linearized quasiparticle energy is the column that equals ε + Z(eQP without Z − ε), with Z ≈ 0.76, although the archive's README calls that column "eQP (Z = 1)". In the QSGW files, the first energy column holds the QSGW levels and the quasiparticle columns do not.
5. Quantifies the ladder correction to the Mg₂Si splitting as about −2 meV (−0.3 to −4.5 meV across the ladder window and the basis), against a lowering of the gap by about 90–125 meV. Contrasts this with QSGW80, which lowers the splitting by about 60 meV: the 0.8 scaling mimics the vertex correction for the gap but overcorrects the splitting, so the remaining QSGW–experiment gap for Mg₂Si is not explained by missing ladder diagrams. States that the Questaal and ecalj QSGW splittings differ by about 28 meV because of the basis, so the correction must be taken within one code.
6. Explains Q:
   - at fixed n, Q equals ⟨q⊥²⟩ at 300 K divided by ⟨q⊥²⟩ at 77 K, which is 300/77 = 3.90 for a parabolic valley;
   - the QSGW boxes give about 4.47 for Mg₂Si because the lower state flattens along the transverse [011] diagonal, where it couples to the upper state; this flattening shows on the box but is missed by X–W scans and fitted masses;
   - degeneracy at 3 × 10¹⁸ cm⁻³ lowers Q to about 3.4;
   - the much smaller splitting of Mg₂Sn makes its lower state far more non-parabolic, so that Q ≈ 6.5;
   - these are leading-order values: the PBE matrix elements show that |p|² saturates at larger q⊥.
7. Assembles the predicted splittings from the QSGW value at the lattice constant of the measurement and the electron–phonon renormalization of X3 − X1:
   - the 100 K lattice raises the Mg₂Si splitting by about 23 meV;
   - the electron–phonon renormalization is about +12 meV at zero point and +27 meV at 300 K for Mg₂Si, and about +20 meV at 300 K for Mg₂Ge;
   - thermal expansion and electron–phonon coupling therefore nearly cancel, and the Mg₂Si splitting decreases by only about 8 meV between T → 0 and 300 K.

   Compares this qualitatively with the Mg₂Ge band of Lott and Lynch, which moves down by about 17 meV between 77 and 300 K.
8. Lists what the data do not support:
   - electron–phonon and lattice effects for Mg₂Sn, and the lattice effect for Mg₂Ge (not archived);
   - absolute band strengths at the QSGW level;
   - convergence with the ladder window;
   - the measured temperature dependence of the band area as a test of Q, because the carrier density of the samples changes with temperature (freeze-out onto donors);
   - spin–orbit effects, which these calculations omit.
