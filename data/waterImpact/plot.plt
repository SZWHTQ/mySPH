#!/usr/bin/env gnuplot
set terminal qt size 1000, 500 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 1000, 500 animate delay 10 font "Times New Roman, 18"

# set output "data.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-0.022:0.60]
set yrange [-0.022:0.376]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:2]
# unset colorbox

do for [i=1:95:1] {
    set title sprintf('%.3fs', 5e-3*i)
    plot sprintf('./output/Type_2_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) pt 7 ps 0.3 palette, \
         sprintf('./output/Type_-2_%d.dat', i) using 4:5 pt 1 ps 1 lt rgb "#ffaa00", \
         sprintf('./output/Type_103_%d.dat', i) using 4:5 pt 7 ps 0.3 lt rgb "#23a9f2"
    pause 0.5
}

pause -1
