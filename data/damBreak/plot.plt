#!/usr/bin/env gnuplot
set terminal qt size 1000, 500 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 1000, 500 animate delay 10 font "Times New Roman, 18"

# set output "velocity.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-5:105]
set yrange [-5:45]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:20]
# unset colorbox

do for [i=1:300:1] {
    set title sprintf('%.2fs', i*0.05)
    plot sprintf('./output/Type_2_%d.dat', i) using 2:3:(($4**2+$5**2)**0.5) pt 6 palette ps 0.5 lw 0.5, \
         sprintf('./output/Type_-2_%d.dat', i) using 2:3 pt 1 ps 0.5 lw 0.2 lt rgb "orange"
    pause 0.02
}

pause -1
