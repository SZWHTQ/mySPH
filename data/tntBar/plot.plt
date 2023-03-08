#!/usr/bin/env gnuplot
set terminal qt size 1000, 800 position 200, 50 font "Times New Roman, 10"
# set terminal gif size 1000, 800 enhanced animate delay 10 font "Times New Roman, 10"
# set terminal png size 1000, 800 font "Times New Roman, 10"

# set output "data.gif"
set grid
set size 1, 1
unset key

deltaT   = 1e-9
timeStep = 100
maxStep  = 14400

do for [i=0:maxStep/100:timeStep/100] {
    # set output sprintf('./pic/%.3fμs.png', i*deltaT*1e6)
    set multiplot title sprintf('%.3fμs',  i*deltaT*1e6) font ", 15" #layout 2, 2

    set size 0.5, 0.48
    set xlabel "X/m"
    set xrange [0:0.1]

    set origin 0, 0.48
    set title 'Pressure'
    set ylabel "P/GPa"
    set yrange [0:25]
    plot sprintf('./output/Reconstruct%d.dat', i) u 2:($6*1e-9) with lines lw 2 lt rgb "black"

    set origin 0.5, 0.48
    set title 'Density'
    set ylabel "ρ/kg·m^{-3}"
    set yrange [1.2e3:2.4e3]
    plot sprintf('./output/Reconstruct%d.dat', i) u 2:5 with lines lw 2 lt rgb "black"

    set origin 0, 0
    set title 'Velocity'
    set ylabel "v/m·s^{-1}"
    set yrange [0:2e3]
    plot sprintf('./output/Reconstruct%d.dat', i) u 2:3 with lines  lw 2 lt rgb "black"#, \
       # sprintf('./output/Reconstruct%d.dat', i) u 2:3 with points pt 1 ps 0.5 lt rgb "black"

    set origin 0.5, 0
    set title 'Internal Energy'
    set ylabel "E/MJ"
    set yrange [3:7]
    plot sprintf('./output/Reconstruct%d.dat', i) u 2:($7*1e-6) with lines lw 2 lt rgb "black"

    pause 0.2
    
    unset multiplot
}

pause -1