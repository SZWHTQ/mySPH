#!/usr/bin/env gnuplot
set terminal qt size 750, 700 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-0.58:0.58]
set yrange [-0.58:0.58]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:1e3]
set cbtics font ",12"
set cblabel "Pressure/MPa" offset 0,0

do for [i=0:250:1] {
    set title sprintf('%.3fμs', i*1e-5*1e6)
    plot sprintf('./output/Type_6_%d.dat', i) using 2:3:($8*1e-6) palette pt 6 lw 0.4 , \
         sprintf('./output/Type_5_%d.dat', i) using 2:3:($8*1e-6) palette pt 6 lw 0.4 , \
        #  sprintf('./output/Type_-6_%d.dat', i) using 2:3 pt 6 lw 0.4 lt rgb "green"
    pause 0.2
}

pause -1
