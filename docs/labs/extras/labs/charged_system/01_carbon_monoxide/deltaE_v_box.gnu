#!/usr/bin/gnuplot
set key top center font ",16"
set xrange [0:0.1]
set xlabel "1/L (Ang)" font ",22"
set ylabel "E_{N-1} - E_{N} (eV)" font ",22" off -4,0
set lmargin 20
set xtics font ",22"
set ytics font ",22"
f(x) = a + b*x
fit f(x) 'FILENAME' u (1/$1):(($4-$2)*13.60569) via a,b
plot 'FILENAME' u (1/$1):(($4-$2)*13.60569) w p ps 2 pt 7 title "DATA", f(x) title "FIT"
pause mouse
