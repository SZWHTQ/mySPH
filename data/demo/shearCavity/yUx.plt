#!/usr/bin/env gnuplot
set terminal qt size 960, 720 position 200, 50 font "Times New Roman, 24" dashlength 1


set xlabel "y/L"
set ylabel "V/V_{drive}"
set xrange [0:1]
set yrange [-0.35:1]
set grid
set size 1, 1
set key box left top at 0.04, 0.923 width 3 height 0.5 opaque

do for [i=0:500:1] {
    set title sprintf('%d', i)
    plot sprintf('./output/Type_3_%d.dat', i) u ($5/1e-3):(abs(($4/1e-3)-0.5) > 0.025 ? 1/0 : ($6/1e-3)) smooth bezier lw 2 lt rgb "black" title "SPH"
    pause 0.02
}
pause -1