Zenodo record 16927281, "Hinge electronic structure of strained half-Heuslers" (CC BY 4.0), file strained_hH_zenodo.zip.
Record description (verbatim): This repository contains data and code necessary to reproduce the numerical results in:
Sanjib Kumar Das, Ion Cosma Fulga, Rakshanda Dhawan, Hem C. Kandpal, Jeroen van den Brink, and Jorge I. Facio,
"Hinge electronic structure of strained half-Heuslers" arXiv:2504.20190.
The "wannier" folder contains Wannier Hamiltonians as obtained using density-functional theory. The "kwant_ham" folder
contains processed data used to construct the tight-binding models. In the "code" folder we include a python file which
reproduces some of our main results.

Staged here: code/strained_hH.py (unchanged) and the header block of each wannier/+hamdata file (everything before the
first "Tij, Hij:" hopping block).  The hopping blocks and the kwant_ham pickles (the same hoppings, processed by
process_hamdata) are evaluated by the calculation service only.
