Fonder "input_vasp+wannier90" contains VASP and Wannier90 input files in subfolders "ZrSiS", "ZrSiSe" and "ZrSiTe", respectively.

Folder "input-MR" contains WannierTools input files for magnetoresistance and Hall resistivity calculations. 

Folder input-MR contains input files for MR and Hall resistivity calculations using WannierTools. Each subfolder "ZrSiS", "ZrSiSe" and "ZrSiTe" contains input file "input.dat" that can be used to calculate the following properties by using the respective tags in the &CONTROL section:
 - BulkFS_calc: Fermi surface (Fig. 1a, Fig. 4a);
 - Boltz_OHE_calc: Hall resistivity;
 - Boltz_evolve_k: orbits in k- and r-space (Fig. 1b-f, Fig. 4b-d, Fig. 7b,d);
 - BulkFS_plane_calc: cross section of the Fermi surface (Fig. 1b-d, Fig. 4b-d, Fig. 7a,c).
In addition, orientation of the magnetic field B needs to be specified (Fig. 2, Fig. 3, Fig. 5, Fig. 6, Fig. 7a,c). Fermi surface cross sections require to specify two vectors and an initial point to identify the corresponding plane.

Files "wannier90_hr.dat_nsymm16" contain tight-binding Hamiltonians.
