#!/usr/bin/env gnuplot
set terminal qt size 1200, 600 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1200, 450 enhanced animate delay 10 font "Times New Roman, 24"

# set output "stress_xx.gif"

set xlabel "x/m"
set ylabel "y/m"
set xrange [-0.01:0.11]
set yrange [-0.03:0.03]
set xtics 0.02
set ytics 0.01
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [-0.05:0.05]
set cbtics 0.025
set cblabel "Stress_{xx}/MPa"

do for [i=1:50:1] {
    set title sprintf('%.1fms', i*4)
    plot sprintf('./output/Type_102_%d.dat', i) using 4:5:($18*1e-6) palette pt 7 ps 1
    pause 0.1
}

pause -1
