#!/usr/bin/env gnuplot
set terminal qt size 1200, 600 position 200, 50 font "Times New Roman, 24"
# set terminal gif size 1200, 600 enhanced animate delay 10 font "Times New Roman, 24"

scale = 1000

E = 22e9          # N/m^2
nu = 0.3          # Poisson's ratio
h = 0.2          # m
b = 1             # m
l = 1             # m
I = (b*(h**3))/12   # m^4
G = E/(2*(1+nu))  # N/m^2
rho = 2400        # kg/m^3
# q = -9.81*rho*h*b # N/m
F = 8e3           # N
# w(x) = q*x**2/(24*E*I)*(x**2-4*l*x+6*l**2) * (x >= 0 ? x <= 0.1 ? 1 : NaN : NaN )
w(x) = -F*((x**3/(E*I)) + (3*x/(2*G*h*b))) * (x >= 0 ? x <= 1 ? scale : NaN : NaN )

# set output "stress_xx.gif"

set xlabel "x/m"
set ylabel "y/m"
set xrange [-0.1:1.1]
set yrange [-0.3:0.3]
set xtics 0.1
set ytics 0.2
set grid
set size 1, 1
unset key
set palette rgbformulae 22, 13, -31
# set cbrange [-0.05:0.05]
# set cbtics 0.025
# set cblabel "Stress_{xx}/MPa"
# set cbrange [0:0.1]
# set cbtics 0.02
# set cblabel "Mises/MPa"
set cbrange [-1.8e-4:0]
set cbtics 0.3e-4
set cblabel "Displacement/m"

# do for [i=1:5:1] {
#     set title sprintf('%.1fms', i*4)
#     # plot sprintf('./output/Type_105_%d.dat', i) using 4:5:($18*1e-6) palette pt 7 ps 1.2, \
#     # #      w(x) w l lt rgb "#000000" lw 2
#     # plot sprintf('./output/Type_105_%d.dat', i) using 4:5:(sqrt($18**2+$21**2-$18*$21+3*$19**2)*1e-6) palette pt 7 ps 1.2, \
#     #      w(x) w l lt rgb "#bf2029" lw 2
#     plot sprintf('./output/Type_105_%d.dat', i) using 4:5:($22*1e-6) palette pt 7 ps 1.2, \
#          w(x) w l lt rgb "#bf2029" lw 2
#     pause 0.1
# }

do for [i=1:30:1] {
    set title sprintf('%.1fms', i*4)
    plot sprintf('./output/Type_105_%d.dat', i) using ($4+(scale-1)*$16):($5+(scale-1)*$17):($17) palette pt 7 ps 1.2, \
         w(x) w l lt rgb "#bf2029" lw 2
    pause 0.1
}


pause -1
