#!/usr/bin/env gnuplot
set terminal qt enhanced size 500, 400 position 200, 50 font "Times New Roman, 12" dashlength 0.3
# set terminal png enhanced size 500, 400 font "Times New Roman, 12"

# set output "4000_NNPS2_SKF1_SLE0.png"
set grid
set size 1, 1
unset key

deltaT   = 1e-9
timeStep = 100
maxStep  = 14400

set title 'Pressure'
set xlabel "X/m"
set ylabel 'P/GPa'
set xrange [0:0.1]
set yrange [0:23]
set xtics 0.01
set ytics 2
set arrow from 0.085, 22.5 to 0.08, 21 lt rgb "black"
set label "P_{CJ}" at 0.085, 22.5
plot './output/1000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/2000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/3000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/4000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/5000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/6000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/7000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/8000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/9000.dat'  u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/10000.dat' u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/11000.dat' u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/12000.dat' u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/13000.dat' u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     './output/14000.dat' u 2:($6*1e-9) with lines lw 2 lt rgb "black", \
     21 with lines dashtype ' -' lw 2 lt rgb "black"

pause -1