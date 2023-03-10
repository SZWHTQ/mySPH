#!/usr/bin/env gnuplot
set terminal qt size 700, 700 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-0.58:0.58]
set yrange [-0.58:0.58]
set grid
set size 1, 1
unset key

do for [i=0:250:1] {
    set title sprintf('%.3fμs', i*1e-5*1e6)
    plot sprintf('./output/Type_6_%d.dat', i) using 2:3 pt 6 lw 0.2 lt rgb "blue", \
         sprintf('./output/Type_5_%d.dat', i) using 2:3 pt 6 lw 0.2 lt rgb "red", \
        #  sprintf('./output/Type_-6_%,d.dat', i) using 2:3 pt 6 lw 0.2 lt rgb "green"
    pause 0.5
}

pause -1
