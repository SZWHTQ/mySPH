#!/usr/bin/env gnuplot
set terminal qt size 1200, 600 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1200, 600 enhanced animate delay 10 font "Times New Roman, 24"


E = 11.9589e6     # N/m^2
h = 0.004         # m
b = 0.002         # m
l = 0.1           # m
I = (b*h**3)/12   # m^4
rho = 1100        # kg/m^3
q = -9.81*rho*h*b # N/m
w(x, y) = q*x**2/(24*E*I)*(x**2-4*l*x+6*l**2) \
       * (x >= 0 ? x <= 0.1 ? 1 : NaN : NaN ) \
       * (y >=-0.01 ? y <= 0.01 ? 1 : NaN : NaN )

# set output "stress_xx.gif"

set xlabel "x/m"
set ylabel "y/m"
set zlabel "z/m"
set xrange [-0.01:0.11]
set yrange [-0.06:0.06]
set zrange [-0.02:0.02]
set xtics 0.02
set ytics 0.02
set ztics 0.01
set mxtics 2
set mytics 2
set mztics 2
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [-0.05:0.05]
# set cbtics 0.025
# set cblabel "Stress_{xx}/MPa"
set cbrange [0:0.16]
set cbtics 0.02
set cblabel "Mises/MPa"

do for [i=1:150:1] {
    set title sprintf('%.1fms', i*4)
    # plot sprintf('./output/Type_102_%d.dat', i) using 4:5:($18*1e-6) palette pt 7 ps 0.75, \
    #      w(x) w l lt rgb "#000000" lw 2
    splot sprintf('./output/Type_102_%d.dat', i) using 4:5:6:(sqrt(($21+$25+$29)**2-3*($25*$29+$29*$21+$21*$25-$22**2-$23**2-$28**2))*1e-4) palette pt 7 ps 0.5, \
         w(x, y) w l lt rgb "#bf2029" lw 1
    pause 0.1
}

pause -1
