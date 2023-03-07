#!/usr/bin/env gnuplot
set terminal qt size 847, 500 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-10:110]
set yrange [-10:50]
set grid
set size 1, 1
unset key
# set xyplane 0.2
# set hidden3d
# set pm3d

do for [i=0:1e4:50] {
    set title sprintf('%d', i)
    plot sprintf('./output/%d.dat', i) using 2:3:8 pt 6 lt rgb "blue"
    pause 0.1
}

pause -1
