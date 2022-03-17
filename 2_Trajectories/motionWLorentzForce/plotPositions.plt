#system "sh checkCourant.sh"
set for [i=1:8] linetype i dashtype i
set grid
#set logscale x 10
plot   "file2Plot.dat" u 6:7 lt 8 lw 2.2 pi 123 w lp t 'Trajectory'
#replot "file2Plot.dat" u 1:3 lt 6 lw 2.2 pi 173 w lp t 'U(y)'
#replot "file2Plot.dat" u 1:4 lt 4 lw 2.2 pi 223 w lp t 'U(z)'
#plot   "logTime" u 1:2 lt 8 lw 2.2 pi 123 w lp t 'U(x)'
#replot "logTime" u 1:3 lt 6 lw 2.2 pi 173 w lp t 'U(y)'
#replot "logTime" u 1:4 lt 4 lw 2.2 pi 223 w lp t 'U(z)'
replot

