#!/usr/bin/env gnuplot
set terminal qt size 960, 720 position 200, 50 font "Times New Roman, 24" dashlength 1


set xlabel "y/L"
set ylabel "V/V_{drive}"
set xrange [0:1]
set yrange [-0.2:0.25]
set grid
set size 1, 1
set key box left top at 0.04, 0.923 width 3 height 0.5 opaque
set palette rgbformulae 22, 13, -31
set cbrange [-0.02:0.02]

do for [i=0:500:1] {
    set title sprintf('%d', i)
    plot sprintf('./output/Type_3_%d.dat', i) u ($4/1e-3):(abs(($5/1e-3)-0.5) > 0.02 ? 1/0 : ($7/1e-3)):($5/1e-3 - 0.5) pt 6 lw 2 palette title "SPH"
    pause 0.02
}
# plot './output/Type_3_50.dat' u ($4/1e-3):(abs(($5/1e-3)-0.5) > 0.02 ? 1/0 : ($7/1e-3)):($5/1e-3 - 0.5) pt 6 lw 2 lt palette title "SPH"
pause -1