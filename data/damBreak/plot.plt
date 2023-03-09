#!/usr/bin/env gnuplot
set terminal qt size 847, 500 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-5:105]
set yrange [-5:45]
set grid
set size 1, 1
unset key
# set xyplane 0.2
# set hidden3d
# set pm3d

do for [i=1:300:1] {
    set title sprintf('%.2fs', i*0.05)
    plot sprintf('./output/Type_2_%d.dat', i) using 2:3 pt 6 ps 0.5 lw 0.2 lt rgb "blue", \
         sprintf('./output/Type_-2_%d.dat', i) using 2:3 pt 1 ps 0.5 lw 0.2 lt rgb "orange"
    pause 0.05
}

pause -1
