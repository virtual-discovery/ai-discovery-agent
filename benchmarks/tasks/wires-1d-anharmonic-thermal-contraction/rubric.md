# Physical-reasoning rubric

Each item is a checkable statement that the report must make explicitly.

1. Establishes the sign convention of the archived SSCHA stress from evidence in the data (for example the positive kinetic difference between the SSCHA stress and the ab initio average stress that grows with T, or the pressure line of the pw.x output). Concludes that the wires contract on heating (α_L < 0).
2. Uses the SSCHA free-energy stress, not the ab initio ensemble-average stress, for (∂P/∂T)_L. Explains that the average lacks the kinetic (quantum/thermal) contribution.
3. Explains the stiffness route: the long-wavelength longitudinal-acoustic branch of the 1D-sum-rule DFPT force constants. Notes that the plain ph.x matrices violate the acoustic sum rules near Γ, and that the single archived strained calculation (about +4.7 %) gives only a finite-strain secant value. Uses a consistent cell volume, including vacuum, for the stress and the stiffness.
4. Attributes the contraction to transverse flexural and torsional modes (the tension or membrane effect). States that the archive's single strained calculation cannot support a quantitative quasi-harmonic α_L, so the QHA-versus-SSCHA difference is not established by these data.
5. Distinguishes the auxiliary SSCHA frequencies (positive by construction; they define the entropy) from the free-energy Hessian frequencies (which judge stability). Identifies the small imaginary Hessian value at Γ as the axial rotation (1D rotational invariance), not an instability.
6. Derives C_V = Σ c_harm(ω)(1 − (T/ω) ∂ω/∂T) from the SSCHA entropy. States that hardening auxiliary frequencies lower C_V below the harmonic value (CuC2) and softening raises it (TaSe3), with the size of each deviation at 300 K.
7. Treats sampling and convergence explicitly. Gives stress error bars. Identifies runs stopped at the Kong–Liu limit with residual centroid forces (for example CuC2 at 400 K) and avoids or justifies them. Notes the 10k-versus-100k disagreement near the TaSe3 instability. Does not claim a CuC2 instability from small, degeneracy-breaking imaginary flexural values.
