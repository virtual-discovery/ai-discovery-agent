set terminal pdf enhanced color font ",20"
set output 'sigma_ahc.pdf'
set key samplen 0.8
set ylabel offset 0.0,0
set xrange [  -1.00000:   0.00000]
set xlabel "Energy (eV)"
set ylabel "\sigma (1/(Ohm*cm))"
plot 'sigma_ahc.txt' u 1:2 w l title '\sigma_{xy}' lc rgb 'red' lw 4, \
     'sigma_ahc.txt' u 1:3 w l title '\sigma_{yz}' lc rgb 'blue' lw 4, \
     'sigma_ahc.txt' u 1:4 w l title '\sigma_{zx}' lc rgb 'orange' lw 4 
