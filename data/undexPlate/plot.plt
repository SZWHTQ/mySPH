#!/usr/bin/env gnuplot
set terminal qt size 1200, 700 position 200, 50 font "Times New Roman, 24"

set xlabel "x/m"
set ylabel "y/m"
set xtics
set ytics
set xrange [-0.50:0.50]
set yrange [-0.25:0.25]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [-5e2:2e3]
set cbtics 250
set cblabel "Pressure/MPa" offset 1, 0

do for [i=1:200:1] {
    set title sprintf('%.1fμs', i*2.5)
    plot sprintf('./output/Type_104_%d.dat', i) using 4:5            pt 7 ps 0.45 lt rgb "#ffaa00", \
         sprintf('./output/Type_5_%d.dat', i)   using 4:5:($10*1e-6) pt 7 ps 0.45 palette, \
         sprintf('./output/Type_6_%d.dat', i)   using 4:5:($10*1e-6) pt 7 ps 0.45 palette
    pause 0.05
}

pause -1 "press enter to exit"
