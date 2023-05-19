#!/usr/bin/env gnuplot
set terminal qt size 600, 900 position 200, 50 font "Times New Roman, 12"

set xlabel "x/m"
set ylabel "y/m"
# set xtics  font ",12"
# set ytics  font ",12"
set xrange [-0.06:0.06]
set yrange [-0.116:0.116]
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
set cbrange [-0.5:2.5]
set cbtics 0.5
set cblabel "Pressure/GPa" offset 0.2,0

do for [i=1:100:1] {
    set title sprintf('%.3fμs', i*5e-7*1e6)
    plot sprintf('./output/Type_6_%d.dat', i) using 4:5:($10*1e-9) palette pt 7 ps 0.35, \
         sprintf('./output/Type_9_%d.dat', i) using 4:5:($10*1e-9) palette pt 7 ps 0.35, \
         sprintf('./output/Type_104_%d.dat', i) using 4:5 pt 6 ps 0.4 lt rgb "#ffaa00" , \
        #  sprintf('./output/Type_101_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1 # lt rgb "orange"
    # plot sprintf('./output/Type_5_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1 , \
    #      sprintf('./output/Type_6_%d.dat', i) using 4:5:($10*1e-6) palette pt 6 ps 1 lw 1
    pause 0.05
}

pause -1 "press enter to exit"
