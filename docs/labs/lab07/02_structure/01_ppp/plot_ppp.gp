reset

set terminal postscript enhanced color solid "Helvetica, 16"
set output 'ppp.eps'

set encoding iso_8859_1

f(x)=a*x**2+b*x+c
a=1;b=1;c=-150
fit f(x) 'e_v_theta.txt' u 1:2 via a,b,c

set ylabel 'Energy (Ry)'
set xlabel '{/Symbol q} ({\260})'

set format y "%5.4f"

p 'e_v_theta.txt' u 1:2 ti '' w p pt 7 ps 1.5 lt -1,\
f(x) ti sprintf("y=ax^2+bx+c, a=%.3E; b=%.3E; c=%.3E", a, b, c) lt -1

