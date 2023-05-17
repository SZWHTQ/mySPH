#!/usr/bin/env gnuplot
set terminal qt size 900, 850 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 900, 850 enhanced animate delay 10 font "Times New Roman, 24"

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
set cbrange [-0.04:0.04]
set cbtics 0.005
set cblabel "Stress_{xx}/MPa"

do for [i=0:100:1] {
    set title sprintf('%.1fms', i*4)
    plot sprintf('./output/Type_102_%d.dat', i) using 4:5:($18*1e-6) palette pt 6 ps 1
    pause 0.1
}

pause -1
