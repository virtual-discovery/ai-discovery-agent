# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Explains why three-phonon scattering is weak in θ-TaN (large acoustic–optical gap freezing out aao processes, acoustic bunching restricting aaa processes) and notes that, with the archived narrow broadening, three-phonon-only solutions leave modes without allowed processes and do not converge, so a three-phonon-only conductivity is not a meaningful reference.
2. Quantifies the four-phonon contribution (it roughly halves κ at 0–100 GPa) and its pressure trend, distinguishing it from the three-phonon trend.
3. Links the growth of phonon–electron scattering above ~20 GPa to the electronic structure of the semimetal under pressure (states at the Fermi level growing), and explains the non-monotonic κ(P) as stiffening/weakening phonon–phonon scattering competing with rising phonon–electron scattering; states the phonon–electron share at 80 GPa.
4. Reports how the calculation settings were fixed (archived CONTROL: 16×16×16 grid, scalebroad 0.1) and checked (mesh or broadening sensitivity), and the size of the RTA-vs-iterative difference (iterative ≈ 2× RTA).
5. Describes the archived BTE.sr_phel_metal layout: branch-major rows over the ShengBTE irreducible list (third grid index fastest), a frequency column in rad/ps that comes from a different phonon calculation (agreement only ~1%), and rates 1/τ in ps⁻¹; and validates the reading against the lab's own phonon–electron channel at a released pressure.
6. Obtains the withheld 40/80 GPa ingredients by transferring scattering rates (mode by mode, in pressure) and solving one combined BTE, rather than interpolating or Matthiessen-combining conductivities; transfers four-phonon rates to 600 K with a calibrated temperature factor rather than reusing 300 K rates; gives an uncertainty for these transfers.
7. States what the data do not support: e.g. the exact pressure of the κ maximum (20 GPa spacing), phonon–electron rates at other temperatures, and convergence beyond the archived 16×16×16 grid for the four-phonon and phonon–electron channels.
