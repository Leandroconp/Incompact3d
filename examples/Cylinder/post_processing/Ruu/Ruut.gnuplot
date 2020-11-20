set term postscript eps enhanced color dashed
 #set logscale y
 set key right top
 set xrange [0:4.5]
 set xlabel "{/Symbol l}_{x_3}/D"
 set ylabel "Number of scans (%)"
 set xtics 0.5
 set output "Histogram.eps"
 set yrange [0.0:8.0]
 set ytics 0.5
 #set grid x y
 #set sample 5000
 plot 'hist' using 1:2 with boxes lt 1 lw 3 lc 0 ti 'fixed cylinder, Re=1250', \
 'Mansy.csv' using 1:2 with lines lt 1 lw 3 lc rgb "blue" ti 'Mansy et al (1994), Re=600'
  #plot 'hist' every 3::0 using 1:2 with lines lt 1 lw 3 lc 0 noti smooth cspline, \ 
