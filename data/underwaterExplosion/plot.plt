#!/usr/bin/env gnuplot
set terminal qt size 800, 800 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [0.045:0.055]
set yrange [-0.005:0.005]
set grid
set size 1, 1
unset key

do for [i=0:12:1] {
    set title sprintf('%.3fμs', i*1e-7*1e6)
    plot sprintf('./output/%d.dat', i) using 2:3 pt 6 lw 1 lt rgb "blue"
    pause 0.5
}

pause -1
