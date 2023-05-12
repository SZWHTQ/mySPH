#!/usr/bin/env gnuplot
set terminal qt size 1000, 700 position 200, 50 font "Times New Roman, 18"
# set terminal gif size 1000, 500 animate delay 10 font "Times New Roman, 18"

# set output "velocity.gif"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [42:58]
set yrange [-1:11]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [0:10]
set cbrange [0:2]
# unset colorbox

do for [i=30:400:1] {
    set title sprintf('%.2fs', i*0.05)
    plot sprintf('./output/Type_2_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) pt 7 lt rgb "#23a9f2" ps 1.2 lw 1, \
         sprintf('./output/Type_8_%d.dat', i) using 4:5:16 pt 7 palette ps 1.2 lw 1, \
         sprintf('./output/Type_-2_%d.dat', i) using 4:5 pt 2 ps 1.2 lw 1 lt rgb "orange"
    # pause 0.02
}

pause -1
