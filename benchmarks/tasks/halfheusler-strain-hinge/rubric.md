# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Shows that both `strained_setup` models have the equilibrium cubic geometry: they are a body-centred tetragonal description with c = √2·a_t = c₀ and the same cell volume, so δ = 0. Supports this with the intact fourfold Γ8 level and with LiBiZn bands that agree with the `orig` model to a few meV. Concludes that the archive contains no strained Hamiltonian, so the paper's strain effects (Dirac or Weyl nodes, strain shifts of the hinge energy) cannot be tested with these data.
2. Notes that the `strained_setup` files use a different Cartesian frame (rotated by 45° about z) and a different lattice basis, and builds the same physical prism there with an explicit integer basis (for example [[1,0,0],[1,1,1],[1,0,1]]), checked through the Zn site's reduced coordinates. Explains that reusing the `orig` integers gives a different prism, with (100) facets and many facet states in the gap.
3. States that the FPLO energy zero lies at an arbitrary point of the LiSbZn gap and differs between the two LiSbZn calculations (about 47 and 289 meV above the VBM). Notes that these two calculations of the same geometry also differ in the fundamental gap (about 344 and 252 meV) and in the depth of the L valley. Attributes this to a difference between the calculations, not to strain.
4. Establishes the band order from the Zn-4s weight at Γ: LiSbZn is normal (Γ6 above Γ8) and LiBiZn is inverted (Γ6 about 1.7 eV below Γ8). Uses a search of the full Brillouin zone to find the LiSbZn conduction-band minimum at or next to L rather than at Γ, which makes the gap indirect.
5. Shows that LiBiZn is a compensated semimetal: band 7 dips about 0.87 eV below the Γ8 maximum near L. Obtains the six-electron count (Li 1, Zn 2 and X-p 3; Li has no Wannier function). From it, derives E_F about 0.40 eV below Γ8 (close to the FPLO zero) on a converged mesh, and notes that coarse meshes are off by tens of meV.
6. Locates the LiSbZn hinge state at the hinge between the two Zn-terminated {111} facets, from the corner weights and the cell geometry (Zn at +1/4 along both sides), and reports its roughly two-thirds Zn-4s weight (cation dangling-bond character). Notes the conflict with the paper's "Bi/Sb hinge, about 95 % on Bi/Sb". Explains it, for instance by noting that the authors' kwant `HoppingKind` assignment reverses the hopping vectors relative to the site positions, which mirrors their real-space picture, or argues why the service geometry is self-consistent.
7. Converges the k = 0 hinge level with prism size and notes that the convergence is slow (about 8 meV above the limit at N = 12 and 3 meV at N = 16). Separates hinge from facet and bulk states by region weights and size scaling. Describes the dispersion:
   - the band is Kramers-degenerate at k = 0 and π, with its minimum at π;
   - two spin-split branches peak near k ≈ 1;
   - the band stays in the projected gap but rises far above the global conduction-band minimum;
   - it never connects to the valence band.
8. Argues that the hinge states are trivial, not higher-order topological:
   - the band is detached from the valence band;
   - it appears at only one of the four hinges, so it depends on the termination;
   - an attractive potential of about −0.38 eV on the hinge cells pushes it into the valence band while pulling new bound states out of the conduction band;
   - the bulk is normal-ordered and the facets are gapped.
9. Explains the visibility contrast. In LiSbZn the bulk and the facets are gapped, so the hinge state is isolated. In LiBiZn, compensated bulk pockets and surface states lie at the same energies, so hinge-weighted states hybridize with them. Also notes what is not established: the strain dependence, surface reconstructions, and the fidelity of the 8-band Wannier models (no DFT bands are archived).
