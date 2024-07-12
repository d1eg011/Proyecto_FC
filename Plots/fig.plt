set terminal epslatex color size 6.0in,5.0in standalone font "" 14
set output 'Escalability.tex'

set border linewidth 6 

set ylabel '\textbf{Speedup} $S$'
set xlabel '\textbf{Número de hilos} $P$'

set ytics scale 2
set xtics scale 2

f(x) = x

plot 'speedup.dat' u 1:2 w p pt 7 ps 2.5 lc rgb "red" t 'Escalabilidad \texttt{Ising.cpp}',\
     f(x) w l lw 3 lc rgb "red" t 'Escalabilidad Lineal'

set output
system('latex Escalability.tex')
system('dvips Escalability.dvi')
system('ps2pdf Escalability.ps')
system('rm Escalability.tex Escalability.log Escalability.aux Escalability-inc.eps Escalability.dvi Escalability.ps')
