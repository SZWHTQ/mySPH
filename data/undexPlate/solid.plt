#!/usr/bin/env gnuplot
set terminal qt size 900, 900 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1000, 700 animate delay 10 font "Times New Roman, 18"

# set output "solid.gif"

array center = [-0.15,0]
startTime = 28
# array center = [-0.05,0]
# startTime = 12
# array center = [-0.025,0]
# startTime = 7

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
set cbrange [-1e-3:1e-5]
# unset colorbox

do for [i=startTime:200:1] {
    set title sprintf('%.1fμs', (i-startTime)*2.5)
    plot sprintf('./output/Type_104_%d.dat', i) using 4:5:16 pt 7 ps 1 palette, \
         sprintf('./output/Type_6_%d.dat', i)   using 4:5    pt 7 ps 1 lt rgb "#7bbeff", \
         sprintf('./output/Type_5_%d.dat', i)   using 4:5    pt 7 ps 1 lt rgb "#bf2029"
    pause 1
}

pause -1
