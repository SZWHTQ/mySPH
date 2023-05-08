#!/usr/bin/env gnuplot
set terminal qt size 900, 850 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 900, 850 enhanced animate delay 10 font "Times New Roman, 18"

# set output "velocity.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-0.01:0.11]
set yrange [-0.058:0.062]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:250]
set cbtics 25 font ",12"

do for [i=0:100:1] {
    set title sprintf('%.1fms', i*0.1)
    plot sprintf('./output/Type_8_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) palette pt 6 ps 1
    pause 0.05
}

pause -1