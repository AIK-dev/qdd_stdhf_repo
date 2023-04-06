set terminal postscript eps enhanced color solid 'bold' 30
set border 15 lw 3
set output 'C60-spectr.eps'
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
set size 1.03,1.83
set ls 1 lt 1 lw 7
set ls 2 lt 2 lc rgb '#00aa00' lw 7
set multiplot
#
set size 0.8,0.8
set origin 0.2,0.2
set xlabel "energy [Ry]"
set ylabel "dipole sttrength"
set yrange [0:]
p "pspectr.C60" u 2:($4*1e3) not w l lw 5 
#
set size 0.8,0.8
set origin 0.2,1.0
set xlabel ""
set form x ""
set ylabel "dipole sttrength"
set yrange [0:]
p "xx3" u 2:($4*1e3) not w l lw 5 
#
!eps2pdf -f C60-spectr.eps
!rm -f C60-spectr.eps
#!open -a Preview C60-spectr.pdf
!evince C60-spectr.pdf
