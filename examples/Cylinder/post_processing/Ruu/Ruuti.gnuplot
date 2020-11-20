set term postscript eps enhanced color dashed
 set key right top
 set xrange [0:3.0]
 set xlabel "l/D"
 set ylabel "R_{u'u'}(l)"
 set xtics 0.5
 set output "Ruuti.eps"
 set yrange [-0.6:1.05]
 set ytics 0.2
 plot 'Ruut93' using 1:2 with lines lt 1 lw 2 lc 0 ti 't=93', \
 'Ruut97' using 1:2 with lines lt 1 lw 2 lc rgb "blue" ti 't=97'
