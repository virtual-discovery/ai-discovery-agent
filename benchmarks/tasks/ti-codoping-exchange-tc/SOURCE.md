# Source and transformations

- **Source:** R. Mozumder, J. Wasmer, D. Antognini Silva, S. Blügel, P. Rüßmann, *High-throughput magnetic co-doping and design of exchange interactions in topological insulators*, Materials Cloud Archive 2024, record 6txnj-70460, DOI 10.24435/materialscloud:b7-6k; preprint arXiv:2407.04413.
- **License:** CC BY 4.0 (whole record).
- **Files used** (per-file API links; the AiiDA exports, 0.4–46 GB, are not used):

| file | bytes | md5 | sha256 |
|---|---|---|---|
| `Jij_table.csv` | 422,637 | e73c49c9da7b0653f2de3e0e4466eaa4 | 99b1c88612dca79aaaab2b6d042c0eff200de4a191f109f9cf69dc71b9a14238 |
| `Bi2Te3.cif` | 1,539 | f714130613532a4063571c361bc6dc46 | 596ede946ab3770c6633f30bd1283b52bd74af83cb9750164e4828022c036045 |
| `Tcs_Bi2Te3.txt` | 10,047 | 1c8a734b01a5b6c2aa18e95a0e2ddc49 | 7d18bfea4c661d58be6e83aa26d0f7eb6874ef1a0cc24abd7867fd8c0f920b5e |
| `data_analysis_Bi2Te3.ipynb` | 1,505,163 | 690b41fab7da399fd6828b77b36e327a | 82019779c501cd00cf5e85c917fcfed4a2ff985ea756dfb4f928bec6b774a5ea |
| `README.md` | 6,985 | 7b0d99083467b6cc0a2215d8b07873d7 | 28d3f31d0394bb0f41f77bdb1f6f297510dc69515294b5f695588812ddb24c63 |

## What is staged

- `catalog.json`: the CIF lines verbatim, the hexagonal lattice constants derived from it, the dopant list, column descriptions of `Jij_table.csv` (taken from the authors' notebook, which builds the table), the paper's Hamiltonian (Eq. 2) and moment criterion, the purchase types and the compositions.
- Lab `data.json`: all 2077 rows of `Jij_table.csv`, unchanged in value, grouped by dopant pair. A `pair` purchase returns the rows of one pair and one shell window verbatim. Missing calculations stay missing.
- The `monte_carlo` service (`dilute.py`, `mc_dilute.c`) simulates whatever couplings the agent supplies. It contains no archived values.

`Tcs_Bi2Te3.txt` (the authors' mean-field estimates) and the notebook are not staged. They were used only to verify the convention: the archived mean-field matrix equals (1/3k_B) Σ_shells N_shell J_shell, with the authors' (T_AA + T_BB + T_AB)/3 rule for co-doping, to 4·10⁻⁵ K per unit concentration. This confirms the pair-once Hamiltonian H = −½ Σ_{i≠j} J_ij ê_i·ê_j of the paper.

## Archive properties found during construction

- **`m2` column.** The notebook reads the dimer moments as `magmoms[[0, 20]]`. Index 20 is impurity 2 only in the shell-0 cluster. For shells 1–16, `m2` is about 0.01 μB for every pair, including Mn (4.4 μB), while the same species' `m1` in its own pair rows keeps the full moment. `m01` and `m02` are the single-impurity moments.
- **Missing dimers.** For example, the archive lacks Co–Co in shells 0–6 and 16, Mn–Fe in shells 6 and 16, Ti–Ti in shell 2, Fe–Fe in shell 8, Tc–Tc in shell 9 and several Ru pairs. The Co–Co shell-9 run converged to a non-magnetic state (m1 = 0.005 μB).
- Far shells (7–16) were computed only for the 11 magnetic dopants.

## What is hidden

The Curie temperatures of the five compositions (three ordering temperatures, one no-order case, one undetermined case) and the best co-dopant for Cr. These are not archived values. They are consequences of the archived exchange table, computed by the construction scripts (`benchmarks/construction/ti-codoping-exchange-tc/`) with the lab's own Monte Carlo code path, on larger supercells and with more disorder samples than the service allows. The paper reports only mean-field estimates and explicitly leaves Monte Carlo to future work.
