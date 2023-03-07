#!/usr/bin/env gnuplot
set terminal qt size 1200, 800 position 200, 50 font "Times New Roman, 16"

set xrange [-0.4:0.4]
set grid

set multiplot title sprintf('%.3fs', 200*1e-3) font ", 20" #layout 2, 2

set origin 0, 0.49
set size 0.5, 0.49
set ylabel "Pressure"
set xrange [-0.4:0.4]
set grid
plot './exact/Riemann_0.2s.dat' u 1:4 with lines lw 4 lt rgb "black" title "Exact", \
     '~/OneDrive/Essay/Data/shockTube/SPH/1600/output/500.dat' u 2:6 with lines lw 2 lt rgb "green" title "SPH 1600", \
     '~/OneDrive/Essay/Data/shockTube/SPH/800/output/1000.dat' u 2:6 with lines lw 2 dashtype "-" lt rgb "cyan" title "SPH 800", \
     '~/OneDrive/Essay/Data/shockTube/SPH/400/output/1000.dat' u 2:6 with lines lw 2 dashtype "." lt rgb "blue" title "SPH 400", \
     './output/1000.dat' u 2:6 with lines lw 2 lt rgb "red" title "Riemann"

unset key
set origin 0.5, 0.49
set size 0.5, 0.49
set ylabel "Density"
set xrange [-0.4:0.4]
unset yrange
set grid
plot './exact/Riemann_0.2s.dat' u 1:2 with lines lw 4 lt rgb "black" title "Exact", \
     '~/OneDrive/Essay/Data/shockTube/SPH/1600/output/500.dat' u 2:5 with lines lw 2 lt rgb "green" title "SPH 1600", \
     '~/OneDrive/Essay/Data/shockTube/SPH/800/output/1000.dat' u 2:5 with lines lw 2 dashtype "-" lt rgb "cyan" title "SPH 800", \
     '~/OneDrive/Essay/Data/shockTube/SPH/400/output/1000.dat' u 2:5 with lines lw 2 dashtype "." lt rgb "blue" title "SPH 400", \
     './output/1000.dat' u 2:5 with lines lw 2 lt rgb "red" title "Riemann"

set origin 0, 0
set size 0.5, 0.49
set ylabel "Velocity"
set xrange [-0.4:0.4]
set grid
plot './exact/Riemann_0.2s.dat' u 1:3 with lines lw 4 lt rgb "black" title "Exact", \
     '~/OneDrive/Essay/Data/shockTube/SPH/1600/output/500.dat' u 2:3 with lines lw 2 lt rgb "green" title "SPH 1600", \
     '~/OneDrive/Essay/Data/shockTube/SPH/800/output/1000.dat' u 2:3 with lines lw 2 dashtype "-" lt rgb "cyan" title "SPH 800", \
     '~/OneDrive/Essay/Data/shockTube/SPH/400/output/1000.dat' u 2:3 with lines lw 2 dashtype "." lt rgb "blue" title "SPH 400", \
     './output/1000.dat' u 2:3 with lines lw 2 lt rgb "red" title "Riemann"

set origin 0.5, 0
set size 0.5, 0.49
set ylabel "Internal Energy"
set xrange [-0.4:0.4]
unset yrange
set grid
plot '~/OneDrive/Essay/Data/shockTube/SPH/1600/output/500.dat' u 2:7 with lines lw 2 lt rgb "green" title "SPH 1600", \
     '~/OneDrive/Essay/Data/shockTube/SPH/800/output/1000.dat' u 2:7 with lines lw 2 dashtype "-" lt rgb "cyan" title "SPH 800", \
     '~/OneDrive/Essay/Data/shockTube/SPH/400/output/1000.dat' u 2:7 with lines lw 2 dashtype "." lt rgb "blue" title "SPH 400", \
     './output/1000.dat' u 2:7 with lines lw 2 lt rgb "red" title "Riemann"

pause -1
