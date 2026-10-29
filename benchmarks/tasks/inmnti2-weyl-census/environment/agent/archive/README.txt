Summary of README.txt:
- HOW TO CITE
- CONTENTS
- DETAILS
- EXTERNAL REFERENCES


================================

This repository contains the input and output files required to reproduce the results of the following work:
Prediction of a novel type-I antiferromagnetic Weyl semimetal

================================
HOW TO CITE
================================


If you use these data for publication purposes, please cite:
* D. Grassano, L. Binci, N. Marzari (2022) arXiv:2208.11412


================================
CONTENTS
================================

standard/                                     - Calculations on InMnTi2 with/starting from standard DFT
    scf/                                      - DFT Self-consistent calculations
    nscf/                                     - DFT Non-self-consistent calculations on a uniform grid for the Wannierization
    nscfXYZ/                                  - DFT Non-self-consistent calculations on the q_{xyz} direction between Weyl points
    wannier/                                  - Wannierization
    wtool/                                    - Wannier post-processing and topological properties
        bulk/                                 - Bulk properties
        find/                                 - Find Weyl points
        weyl/                                 - Weyl points chirality
        weylXYZ/                              - Weyl point pair band dispersion along q_{xyz} directions
        arpes/                                - Spectral function along surface high-symmetry path
    wberri/                                   - Fermi surface

hubbard/                                      - Calculations on InMnTi2 with/starting from DFT+U
    relax/                                    - DFT+U Relax calculations
    scf/                                      - DFT+U Self-consistent calculations
    nscf/                                     - DFT+U Non-self-consistent calculations
    nscfXYZ/                                  - DFT+U Non-self-consistent calculations on the q_{xyz} direction between Weyl points
    wannier/                                  - Wannierization
        kubo/                                 - Kubo-Greenwood calculations of optical properties
    wtool/                                    - Wannier post-processing and topological properties
        bulk/                                 - Bulk properties
        find/                                 - Find Weyl points
        weyl/                                 - Weyl points chirality
        weylXYZ/                              - Weyl point pair band dispersion along q_{xyz} directions
        arpes/                                - Spectral function along surface high-symmetry path
        arpes_WW/                             - Spectral function along node pair projection direction
        farcs/                                - 2D Fermi surface at different chemical potentials
    wberri/                                   - Fermi surface
README.txt                                    - this README file


=======================================
DETAILS
=======================================

Calculations performed using the following codes:
- DFT: Quantum ESPRESSO 7.2 [1]
- Wannierization: Wannier90 [2]
- Wannier post-processing and topological properties: WannierTools [3]
- Fermi surface: WannierBerri [4]

================================
EXTERNAL REFERENCES
================================

[1] P. Giannozzi, et al. - J. Condens. Matter Phys. 21, 395502 (2009)
[2] G. Pizzi, et al. - J. Condens. Matter Phys. 32, 165902 (2020)
[3] Q. Wu, S. Zhang, H. Song, M. Troyer, A.A. Soluyanov - Comput. Phys. Commun. 224, 405 (2018)
[4] S.S. Tsirkin - Npj Comput. Mater. 7, 1 (2021)