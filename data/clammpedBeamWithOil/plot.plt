#!/usr/bin/env gnuplot
set terminal qt size 1000, 500 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 1000, 500 animate delay 10 font "Times New Roman, 18"

# set output "velocity.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-0.33:0.33]
set yrange [-0.03:0.37]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:20]
# unset colorbox

do for [i=0:400:1] {
    set title sprintf('%.2fs', i*0.05)
    plot sprintf('./output/Type_8_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) pt 7 palette ps 0.5 lw 0.5, \
         sprintf('./output/Type_102_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) pt 7 lt rgb "#23a9f2" ps 0.25 lw 0.25, \
         sprintf('./output/Type_-8_%d.dat', i) using 4:5 pt 1 ps 0.5 lw 0.5 lt rgb "orange"
    pause 1
}

pause -1
