#!/usr/bin/env gnuplot
set terminal qt enhanced size 700, 500 font "Times New Roman, 20"

#set title "Convergence Rate"
set xlabel "1/delta"
set xrange [0:30]
set logscale y
set yrange [1e-10:1]

set ylabel "Error"
set grid
set key font ",16" box at 18, 1e-7

plot "./data" u 1:($6/100) with lines lw 2 lt rgb 'black' title " |Error|", \
     "./data" u 1:($6/100) pt 7 lt rgb 'black' title "", \
     1/x           lw 2 title "  Sublinear", \
     1/(2**x)      lw 2 title "  Linear" , \
     1/(2**(2**x)) lw 2 title "  Superlinear"

pause -1