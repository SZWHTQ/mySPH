#!/usr/bin/env gnuplot
set terminal qt size 900, 850 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 900, 850 enhanced animate delay 10 font "Times New Roman, 18"

# set output "velocity.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-1.24e-2:2.0e-2]
set yrange [-0.732e-2:2.646e-2]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:250]
set cbtics 25 font ",12"

do for [i=1:100:1] {
    set title sprintf('%dμs', i)
    plot sprintf('./output/Type_8_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) palette pt 6 ps 1 , \
         sprintf('./output/Type_-8_%d.dat', i) using 4:5 pt 2 ps 1 lt rgb "orange"
    pause 0.02
}

pause -1
