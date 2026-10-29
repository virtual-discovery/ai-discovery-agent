set style data dots
set nokey
set xrange [0: 4.22053]
set yrange [ -1.41060 : 11.34095]
set arrow from  1.01179,  -1.41060 to  1.01179,  11.34095 nohead
set arrow from  1.47998,  -1.41060 to  1.47998,  11.34095 nohead
set arrow from  2.59484,  -1.41060 to  2.59484,  11.34095 nohead
set arrow from  3.06304,  -1.41060 to  3.06304,  11.34095 nohead
set arrow from  3.75234,  -1.41060 to  3.75234,  11.34095 nohead
set xtics ("K"  0.00000,"G"  1.01179,"A"  1.47998,"K"  2.59484,"H"  3.06304,"M"  3.75234,"L"  4.22053)
 plot "wannier90_band.dat"
