#!/usr/bin/env gnuplot
set terminal qt size 700, 700 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [0:0.001]
set yrange [0:0.001]
set grid
set size 1, 1
unset key

do for [i=100:1e4:100] {
    set title sprintf('%.3fs', i*5e-5)
    plot sprintf('./output/%d.dat', i) using 2:3 pt 6 lt rgb "blue"
    pause 0.5
}
