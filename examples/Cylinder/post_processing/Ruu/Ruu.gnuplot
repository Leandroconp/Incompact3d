 set term postscript eps enhanced color dashed 
 set style line 1 lt 1 lw 2 pt 7 ps 1.5 lc 0  #
 set style line 2 lt 1 lw 2 pt 5 ps 1.5 lc 1 #
 set style line 3 lt 1 lw 2 pt 6 ps 1.5 lc 2  #
 set style line 4 lt 1 lw 2 pt 4 ps 1.5 lc 3  #
 set key right top
 set key font ",20" 
 set key spacing 1.5
 set xrange [0:3.0]
 set xlabel "l/D" font ", 18"
 set ylabel "R_{(u_1,u_1)},R_{(u_2,u_2)}" font ", 18"
 set xtics 0.2 font ", 18"
 set output "R.eps"
 set yrange [-0.3:1.0]
 set ytics 0.1 font ", 18"
 set grid x y
 #set style func linespoints
 plot 'Ruu' using 1:2 with lines ls 1 ti 'R_{u_1,u_1},2DOF-III' , \
      'Ruu' using 1:3 with lines ls 2 ti 'R_{u_2,u_2},2DOF-III'
