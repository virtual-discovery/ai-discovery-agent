# Noncollinear Kagome Magnet Driven by Giant Antisymmetric Exchange

The SpinW script is SpinW_Co3Sn2S2_DMI.m. The results in the paper were obtained by running it with SpinW 3.1.1 in Matlab R2025a.

The diffuse magnetic scattering (DMS) data is provided in two directories narrowscans/ and widescans/.
A simple Mathematica notebook Co3Sn2S2-DMS-analysis.nb may be used to import the data and view it.
The wide angle scans were collected at T=1.5, 165, 180, and 300 K.
The narrow angle scans were collected in a narrower range of sample rotation angles close to the (003) position.

This data was produced from the raw instrumental data using the software package DMCpy.

The Mathematica notebook allows for simple plotting of the DMS data files.
The wide scans were created over several different runs and then combined.
The combined data was then exported using the python script readcombined.py.
This data is imported directly using the import function importWide in the Mathematica notebook.
The wide scans are given by file numbers:
(* T=1.5 *)wide3869 = importWide[importdir <> "3869"];
(* T=165 *) wide3871 = importWide[importdir <> "3871"]; (* combines 3870-3871 *)
(* T=180 *) wide3874 = importWide[importdir <> "3874"]; (* combines 3872-3874 *)
(* T=300 *)wide3875 = importWide[importdir <> "3875"];

The narrow scans were imported more simply from the DMCpy data, using the getDataPy function.
The narrow scans are given by file numbers:
(* T = 140 *)py3896 = Quiet@makeNorm@Normal@getDataPy[3896];
(* T = 170 *)py3885 = Quiet@makeNorm@Normal@getDataPy[3885];
(* T = 175 *)py3879 = Quiet@makeNorm@Normal@getDataPy[3879];
(* T = 178 *)py3882 = Quiet@makeNorm@Normal@getDataPy[3882];
(* T = 185 *)py3883 = Quiet@makeNorm@Normal@getDataPy[3883];
(* T = 192 *)py3887 = Quiet@makeNorm@Normal@getDataPy[3887];
(* T = 200 *)py3884 = Quiet@makeNorm@Normal@getDataPy[3884];
(* T = 220 *)py3903 = Quiet@makeNorm@Normal@getDataPy[3903];
(* T = 250 *)py3905 = Quiet@makeNorm@Normal@getDataPy[3905];

All other information about the data can be found in the paper.

