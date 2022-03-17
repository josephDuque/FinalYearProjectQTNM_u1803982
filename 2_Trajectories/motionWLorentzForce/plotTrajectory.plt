#system "sh checkCourant.sh"
set for [i=1:8] linetype i dashtype i
set grid
set xtics 0.0003
set ytics 0.0003
set ztics 0.0005
set xlabel "Position X"
set ylabel "Position Y"
set zlabel "Position Z"
splot   "file2Plot.dat" u 6:7:8 lt 8 lw 2.2 pi 123 w lp t 'Trajectory'

