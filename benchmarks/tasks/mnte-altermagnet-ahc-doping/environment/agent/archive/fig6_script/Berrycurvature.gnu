set encoding iso_8859_1

#set terminal  pngcairo  truecolor enhanced size 1200, 1000 font "Helvetica,65"
#set output 'test.png'

set terminal postscript eps enhanced color font 'Helvetica,35'
set output 'test.eps'

#set palette maxcolors 3
#set palette rgbformulae 33,13,10 maxcolors 3
#set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#ffffff',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000') maxcolors 100

set palette defined ( 0 '#0090ff',1 '#0fffee',2 '#ffffff',3 '#ffee00',4 '#ff7000') maxcolors 100

unset key
set pm3d
set cbrange [ -10: 10 ]
set view map
set size ratio -1
set border lw 6
set xlabel 'k_x ({\305^{-1}})'
set ylabel 'k_y ({\305^{-1}})'
set colorbox
set xtics 0.1 mirror scale 1.0
set ytics 0.1 mirror scale 1.0
set xrange [-0.15: 0.15] noextend
set yrange [-0.15: 0.15] noextend

set pm3d interpolate 2,2

splot '0_0.450.dat' u 1:2:3 w pm3d

set output 'test2.eps'
splot '90_0.450.dat' u 1:2:3 w pm3d

