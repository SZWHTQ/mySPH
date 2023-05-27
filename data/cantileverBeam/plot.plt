#!/usr/bin/env gnuplot
set terminal qt size 1200, 450 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1200, 450 enhanced animate delay 10 font "Times New Roman, 24"

q = 43.164
l = 0.1
E = 11.9589e6
I = 5.333e-9
w(x) =  - q*x**2/(24*E*I)*(x**2-4*l*x+6*l**2) * (x >= 0 ? x <= 0.1 ? 1 : NaN : NaN )

# set output "stress_xx.gif"

set xlabel "x/m"
set ylabel "y/m"
set xrange [-0.01:0.11]
set yrange [-0.02:0.02]
set xtics 0.02
set ytics 0.01
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [-0.05:0.05]
set cbrange [0:0.12]
set cbtics 0.025
set cblabel "Stress_{xx}/MPa"

do for [i=1:100:1] {
    set title sprintf('%.1fms', i*4)
    # plot sprintf('./output/Type_102_%d.dat', i) using 4:5:($18*1e-6) palette pt 7 ps 0.75, \
    #      w(x) w l lt rgb "#000000" lw 2
    plot sprintf('./output/Type_102_%d.dat', i) using 4:5:(sqrt($18**2+$21**2-$18*$21+3*$19*$20)*1e-6) palette pt 7 ps 0.75, \
         w(x) w l lt rgb "#000000" lw 2
    pause 0.1
}

pause -1
