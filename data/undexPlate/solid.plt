#!/usr/bin/env gnuplot
set terminal qt size 900, 900 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1000, 700 animate delay 10 font "Times New Roman, 18"

# set output "solid.gif"

array center = [-0.075,0]

startTime = 20

set xlabel "x/m"
set ylabel "y/m"
set xrange [-0.2+center[1]:0.2+center[1]]
set yrange [-0.25+center[2]:0.25+center[2]]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [0:10]
set cbrange [-7e-6:1e-6]
# unset colorbox

do for [i=20:200:1] {
    set title sprintf('%.3fμs', i*2.5)
    plot sprintf('./output/Type_6_%d.dat', i) using 4:5 pt 7 lt rgb "#7bbeff" ps 1.5 lw 1, \
         sprintf('./output/Type_101_%d.dat', i) using 4:5:16 pt 7 palette ps 1.5 lw 1, \
         sprintf('./output/Type_5_%d.dat', i) using 4:5 pt 7 ps 1.5 lw 1 lt rgb "#bf2029"
    pause 0.05
}

pause -1
