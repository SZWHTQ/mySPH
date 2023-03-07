#!/usr/bin/env gnuplot
# set terminal pngcairo size 800, 600 font "Times New Roman, 15"
set terminal gif size 800, 600 animate delay 10 font "Times New Roman, 15"

set output "density.gif"

set xlabel "X"
set ylabel "Density"
set xrange [-0.4:0.4]
set yrange [0:1.1]
set grid
set key box left bottom

do for [i=10:1000:10] {
    # set output sprintf('pic/density%.3fs.png', i*1e-3)
    set title sprintf('%.3fs', i*1e-3)
    plot sprintf('~/OneDrive/Essay/Data/shockTube/SPH/1600/output/%d.dat', 2.5*i) u 2:5 with lines lw 2 lt rgb "green" title "SPH 1600", \
         sprintf('~/OneDrive/Essay/Data/shockTube/SPH/800/output/%d.dat',    5*i) u 2:5 with lines lw 2 dashtype "-" lt rgb "cyan" title "SPH 800", \
         sprintf('~/OneDrive/Essay/Data/shockTube/SPH/400/output/%d.dat',    5*i) u 2:5 with lines lw 2 dashtype "." lt rgb "blue" title "SPH 400", \
         sprintf('~/OneDrive/Essay/Data/shockTube/DSPH/output/%d.dat', i) u 2:5 with lines lw 2 lt rgb "grey" title "DSPH", \
         sprintf('./output/%d.dat', 5*i) u 2:5 with lines lw 2 lt rgb "red" title "Riemann"
}