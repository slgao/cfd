set term postscript enhanced color
set xlabel "tend"
set ylabel "time (s)"
set output "notlog.eps"

plot "mg.dat" using 1:2 title 'SOR' with linespoints, \
"mg.dat" using 1:3 title 'Multigrid' with linespoints
set terminal x11
