# Construction of zno-hydrogen-shallow-donor

Released-start task built on Aicher, Berger, Aschauer & Diwald (2025), Materials Cloud record c91vp-p9h03 (CC BY 4.0).

## Rebuild

1. `python3 -B fetch.py` streams `ZnO_ubiquitous.tar.gz` (727.5 MB, md5 d21e5cfb405c798dcd23ad41c9678614) once, verifies the md5, and extracts only members under 12 MiB (no vasprun.xml) into `~/.cache/ai-theorist/zno-hydrogen-shallow-donor/raw`, about 0.73 GB.
   - It retries on HTTP 5xx: the archive's storage returned 503 for about 70 minutes on 2026-10-26.
2. `python3 -B extract.py` writes `pool.json` (3.1 MB): 50 runs, with energies, settings, eigenvalues, core potentials and structures.
3. `python3 -B make_task.py` writes:
   - `environment/agent/catalog.json`;
   - `environment/lab/data.json`;
   - `private/hidden.json` and `private/hidden_values.json`;
   - `solution/ablation_values.json`.

   It also prints the calibration rows.
4. `.venv/bin/python -m benchmarks.tools.check_task zno-hydrogen-shallow-donor --cpuset 4,5,20,21` writes `check.json`.

`model.py` holds the thermodynamics:

- Ewald/Lany–Zunger term, verified against the simple-cubic Madelung constant to 1e-6;
- the authors' atom mapping and alignment;
- band filling;
- chemical-potential limits;
- the Fermi–Dirac charge-neutrality solver.

## Build notes

- **Reproduction.** The numpy re-implementation of the authors' notebook reproduces every formation energy in the paper's Table at E_F = 3.85 eV to about 0.01 eV. The chemical-potential vertex Δμ_O = −3.27 eV, Δμ_H = −0.08 eV comes from the Zn(OH)₂ constraint.
- **Deviation 1: band filling.** The authors' neutral OH_O/H_O energies contain about 0.65 eV of conduction-band filling, which puts OH_O(+/0) 0.53 eV above their own gap. The truth removes it, with the Lany–Zunger/Persson sum using OUTCAR k-weights and the alignment-shifted host CBM. The robustness cases keep the defensible alternatives passing:
  - unshifted CBM;
  - an "effective-mass donor at the CBM" treatment.
- **Deviation 2: Fermi level and concentrations.** These come from parabolic bands (m_e* 0.24, m_h* 0.59), stated in the instruction. The authors' py-sc-fermi supercell DOS gave E_F = 3.85 eV, above the CBM, which is unphysical; it is not used.
- **Deviation 3: dropped design elements.** The candidate design called for several things the archive cannot support, so they are not graded:
  - T,p-dependent gas chemical potentials;
  - a quench/anneal model;
  - complex-dissociation entropy;
  - BC vs AB competition;
  - the O–H zero-point energy.

  The archive holds one H configuration per defect, no phonons and only 0 K energies, and the authors also used 0 K chemical potentials. The H site is graded instead as a geometry classification (BC∥). The reservoir decision survives as the Zn(OH)₂-vs-H₂ limit.
- **Categorical risk.** "OH_O⁰ is the g = 1.96 centre" is guessable (H as a shallow donor in ZnO is famous). It is one check of 14; the numerical checks carry the grade.

## Calibration (full table in design.json → calibration)

| field | reference | robustness range | nearest wrong path | tolerance |
|---|---|---|---|---|
| gap | 3.2485 | 3.2466 (bulk 8×8×5) | 0.716 (PBE) | ±0.03 |
| Δμ_H | −0.078 | −0.078 | 0.000 (H₂ limit) | ±0.03 |
| OH_O(+/0) | 3.131 | 3.070–3.248 | 2.84 (equal k-weights); 3.78 (no band filling, outside the gap) | [3.00, 3.30] |
| H_O(+/0) | 3.213 | 3.181–3.248 | 3.85 (no band filling) | [3.00, 3.30] |
| V_O(2+/0) | 2.596 | 2.596–2.603 | 2.691 (adjacent (2+/+)) | ±0.05 |
| E_F (1000 K) | 3.182 | 3.176–3.195 | 3.224 (H₂ limit) | ±0.025 |
| log₁₀[OH_O⁺] | 18.843 | 18.843–18.905 | 19.023 (H₂ limit) | ±0.12 |
| log₁₀[OH_O⁰] | 19.098 | 18.507–19.409 | 20.568 (equal weights); 15.819 (no band filling) | ±0.7 |
| log₁₀[V_O] | 15.951 | 15.951 | 18.501 (PBE) | ±0.15 |

The tightest margins:

- E_F: the H₂-reservoir path misses by 0.042 against a tolerance of 0.025, while the robustness spread reaches 0.013.
- log₁₀[OH_O⁰]: the at-CBM robustness case sits 0.59 dex from the reference, inside the ±0.7 window.

## Oracle

- **Purchases:** 31 of 35 units, with 0 failures:
  - 25 runs: the host, six reference phases and all 18 hybrid defect cells;
  - 6 eigenvalue sets: the host, OH_O 0/+1, H_O 0, V_O +1 and Zn_i +1.
- **Runtime:** about 5 s in the gate container.
- **Results:** the gate passed with all 11 ablations failing and all 6 robustness cases passing. The Harbor oracle scored reward 1.0 (runs/build-zno-hydrogen-shallow-donor).

## Risks

- The task depends on agents recognizing band filling from the moments and occupations. If they skip the eigenvalue purchases, the OH_O/H_O ladders and [OH_O⁰] fail, which is the intended discrimination.
- The ε = 10.4 and the m* values are given; ε is the authors' choice, while the experimental static ε of ZnO is about 8.5.
- `pool.json` (3.1 MB) is kept in the construction folder so the task can be rebuilt without the 727 MB download.
