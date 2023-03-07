#!/usr/bin/env gnuplot
set terminal qt size 1000, 400 font "Times New Roman, 20"


set multiplot

set size 0.5, 1
set xlabel "X/m"
set xrange [-1:1]
set grid
set key font ",16"

set origin 0, 0
# set yrange [-1.5:1]
set ylabel "f(x)"
set key top left box at -0.75, 3
plot './SPH.dat'  u 2:3 with lines lw 2 lt rgb 'black' title " Exact", \
     './SPH.dat'  u 2:4 lw 1.5 pt 1 title " SPH" , \
     './CSPH.dat' u 2:4 lw 1.5 pt 2 title " CSPH", \
     './DSPH.dat' u 2:4 lw 1.5 pt 4 title " DSPH"

set origin 0.5, 0
# set yrange [-12:2]
set ylabel "f'(x)"
set key bottom left box at -0.75, -27.5
plot './SPH.dat'  u 2:6 with lines lw 2 lt rgb 'black' title " Exact", \
     './SPH.dat'  u 2:7 lw 1.5 pt 1 title " SPH" , \
     './CSPH.dat' u 2:7 lw 1.5 pt 2 title " CSPH", \
     './DSPH.dat' u 2:7 lw 1.5 pt 4 title " DSPH"

pause -1