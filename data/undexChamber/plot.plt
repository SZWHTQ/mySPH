#!/usr/bin/env gnuplot
set terminal qt size 750, 700 position 200, 50 font "Times New Roman, 18"

set xlabel "x/m"
set ylabel "y/m"
set xtics  font ",12"
set ytics  font ",12"
set xrange [-0.55:0.55]
set yrange [-0.55:0.55]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [-1e2:8e2]
set cbtics 1e2
set cblabel "Pressure/MPa" offset 0,0

do for [i=1:100:1] {
    set title sprintf('%.3fμs', i*2.5e-5*1e6)
    plot sprintf('./output/Type_-6_%d.dat', i) using 4:5  pt 6 ps 1 lw 1 lt rgb "orange", \
         sprintf('./output/Type_5_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1 , \
         sprintf('./output/Type_6_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1 , \
        #  sprintf('./output/Type_101_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1 # lt rgb "orange"
    # plot sprintf('./output/Type_5_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1 , \
    #      sprintf('./output/Type_6_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1
    pause 0.5
}

pause -1 "press enter to exit"
