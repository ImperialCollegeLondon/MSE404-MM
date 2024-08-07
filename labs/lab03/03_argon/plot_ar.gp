reset

set terminal postscript enhanced color solid 'Helvetica, 14'
set output 'ar_dimer.eps'
set encoding iso_8859_1

set title 'Ar dimer'
set xlabel 'r ({\305})'
set ylabel 'E (Ry)'

set zeroaxis lt -1

e = system("tail -n 1 01_lda/etot_v_r.dat | awk '{print $2}'")
e_lda = e + 0
e = system("tail -n 1 02_pbe/etot_v_r.dat | awk '{print $2}'")
e_pbe = e + 0

p '01_lda/etot_v_r.dat' u 1:($2-e_lda) ti 'LDA' w lp pt 5 lw 2,\
'02_pbe/etot_v_r.dat' u 1:($2-e_pbe) ti 'PBE' w lp pt 7 lw 2
