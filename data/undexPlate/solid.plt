#!/usr/bin/env gnuplot
set terminal qt size 800, 880 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1000, 700 animate delay 10 font "Times New Roman, 18"

# set output "solid.gif"

distance  = 0.3
thickness = 0.015

array center = [-distance/2-thickness/2,0]
startTime = 45
# startTime = 16

set xlabel "x/m"
set ylabel "y/m"
set xrange [-0.2+center[1]:0.2+center[1]]
set yrange [-0.25+center[2]:0.25+center[2]]
set xtics 0.1
set ytics 0.1
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [0:10]
set cbrange [-2.5e-2:0]
# unset colorbox

do for [i=startTime:200:5] {
    set title sprintf('%.1fμs', i*2)
    plot sprintf('./output/Type_104_%d.dat', i) using 4:5:16 pt 7 ps 0.75 palette, \
         sprintf('./output/Type_6_%d.dat', i)   using 4:5    pt 7 ps 0.75 lt rgb "#7bbeff", \
         sprintf('./output/Type_5_%d.dat', i)   using 4:5    pt 7 ps 0.75 lt rgb "#bf2029"
    # pause 0.02
}

pause -1
