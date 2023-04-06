set term post eps enhanced color solid "bold" 30
set out "C3-testinit.eps"
set size 2.03,1.23
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set multiplot
set key top left
#------------------------------------------------------------
set size 0.8,1.0
set origin 0.2,0.2
set xrange [0:3950]
set xtics 1000
set form y "10^{%T}"
set yrange [2e-7:1.4e2]
set xlabel "iteration number"
set ylabel "variance [eV]" offset 0
set log y
unset key
set label "{/=40 C_3}" at graph 0.88,0.93
p "init_ho_deformed/infosp.C3"  u 2:($5*13.6) t "h.o.,{/Symbol= b}=0.9,{/Symbol= g}=0,12 states" w l lw 5,\
 "init_ho_occ/infosp.C3"  u 2:($5*13.6) t "h.o.,{/Symbol= b}=0.3,{/Symbol= g}=30,12 states" w l lc rgb '#00aa00' lw 5,\
 "init_ho_unocc/infosp.C3"  u 2:($5*13.6) t "h.o.,{/Symbol= b}=0.3,{/Symbol= g}=30,18 states" w l lw 5,\
 "init_lcao/infosp.C3"  u 2:($5*13.6) t "localized,12 states" w l lw 5
#------------------------------------------------------------
set auto
unset label
set size 0.8,1.0
set origin 1.2,0.2
set xrange [0:3950]
set xtics 1000
set form y "%g"
set yrange [-433:-381]
set xlabel "iteration number"
set ylabel "energy [eV]" offset 0.9
unset log y
set key top right
set arrow from 0,-31.6710020*13.6 to 3950,-31.6710020*13.6 nohead lt 0 lw 2
p "init_ho_deformed/infosp.C3"  u 2:($3*13.6) t "h.o.,{/Symbol= b}=0.9,{/Symbol= g}=0,12 states" w l lw 5,\
 "init_ho_occ/infosp.C3"  u 2:($3*13.6) t "h.o.,{/Symbol= b}=0.3,{/Symbol= g}=30,12 states" w l lc rgb '#00aa00' lw 5,\
 "init_ho_unocc/infosp.C3"  u 2:($3*13.6) t "h.o.,{/Symbol= b}=0.3,{/Symbol= g}=30,18 states" w l lw 5,\
 "init_lcao/infosp.C3"  u 2:($3*13.6) t "localized,12 states" w l lw 5
#
!eps2pdf -f C3-testinit.eps
!rm -f C3-testinit.eps
!evince C3-testinit.pdf
q
