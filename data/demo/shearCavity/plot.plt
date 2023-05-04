#!/usr/bin/env gnuplot
set terminal qt size 700, 700 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [0:0.001001]
set yrange [0:0.001001]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [0:1e-3]
set cbtics 1e-4 font ",12"

do for [i=0:500:1] {
    set title sprintf('%.3fs', i*5e-3)
    # plot sprintf('./output/Type_3_%d.dat', i) using 4:5:(($6**2+$7**2)**0.5) pt 6 lt palette#, \
        #  sprintf('./output/Type_-3_%d.dat', i) using 4:5 pt 2 lt rgb "orange"
    plot sprintf('./output/Type_3_%d.dat', i) using 4:5:($6/5):($7/5) with vectors head lt rgb 'black'
    pause 0.05
}
