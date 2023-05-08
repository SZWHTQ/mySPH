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
plot './output/Type_5_10.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_20.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_30.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_40.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_50.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_60.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_70.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_80.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_90.dat'  u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_100.dat' u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_110.dat' u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_120.dat' u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_130.dat' u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     './output/Type_5_140.dat' u 4:($8*1e-9) with lines lw 2 lt rgb "black", \
     21 with lines dashtype ' -' lw 2 lt rgb "black"

pause -1