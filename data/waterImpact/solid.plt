#!/usr/bin/env gnuplot
set terminal qt size 1000, 700 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 1000, 700 animate delay 10 font "Times New Roman, 18"

# set output "solid.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [0.236:0.348]
set yrange [-0.012:0.085]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [0:10]
set cbrange [-1e-2:5e-2]
# unset colorbox

do for [i=20:200:1] {
    set title sprintf('%.3fs', i*5e-3)
    plot sprintf('./output/Type_2_%d.dat', i) using 4:5 pt 7 lt rgb "#23a9f2" ps 1.5 lw 1, \
         sprintf('./output/Type_103_%d.dat', i) using 4:5:16 pt 7 palette ps 1.5 lw 1, \
         sprintf('./output/Type_-2_%d.dat', i) using 4:5 pt 2 ps 1.5 lw 1 lt rgb "#ffaa00"
    pause 0.1
}

pause -1
