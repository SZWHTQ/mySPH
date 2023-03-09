#!/usr/bin/env gnuplot
set terminal qt size 900, 600 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-1.24e-2:2.0e-2]
set yrange [-0.732e-2:2.646e-2]
set grid
set size 1, 1
unset key

do for [i=1:100:1] {
    set title sprintf('%dμs', i)
    plot sprintf('./output/Type_8_%d.dat', i) using 2:3:8 pt 6 ps 0.5 lt rgb "blue", \
         sprintf('./output/Type_-8_%d.dat', i) using 2:3:8 pt 2 ps 0.5 lt rgb "orange"
    pause 0.2
}

pause -1
