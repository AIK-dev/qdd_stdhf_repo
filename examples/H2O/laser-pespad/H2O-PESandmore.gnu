# sources in examples/H2O/laser-pespad
set terminal postscript eps enhanced color "TimesNewRoman,12"
set colorsequence classic
set output "H2O-PESandmore.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.05,1.39
set pm3d map interpolate 3,3
set multiplot
#
set origin 0.05,0.06
set size 0.4,0.4
set xrange [0:3*13.6]
set yrange [0:180]
set xlabel "kinetic energy E_{kin} [eV]" offset 0
set ylabel "angle {/Symbol= q} [^o]" offset -0.5
unset key
set log cb
set form cb "10^{%T}"
set cblabel "PES yield [arb.u.]" offset 0.5
set cbrange [1e-9:4e-2]
set title "PES in plane of E_{kin} and {/Symbol= q}"
sp "distavphi.H2O" u ($1*13.6):2:3 
#
set origin 0.05,0.46
set size 0.4,0.4
set xrange [0:3*13.6]
set yrange [0:180]
set xlabel "" offset 0
set ylabel "angle {/Symbol= q} [^o]" offset -0.5
unset key
set log cb
set form cb "10^{%T}"
set cblabel "PES yield [arb.u.]" offset 0.5
set cbrange [5e-10:4e-3]
set title "smoothed PES in plane of E_{kin} and {/Symbol= q}"
sp "parametric.H2O" u ($1*13.6):2:3 
#
set auto
set pm3d clip4in
set origin 0.5,0.06
set size 0.52,0.4 
set xrange [-1.8:1.8]
set yrange [0:1.8]
set xlabel "momentum p_z [a@_0^{-1}]" offset 0
set ylabel "momentum p_{x,y} [a@_0^{-1}]" offset 0
set ytics 0.5
unset key
set log cb
set form cb "10^{%T}"
set cblabel "PES yield [arb.u.]" offset 0.5
set cbrange [1e-9:2e-2]
set title "PES in plane of p_z and p_{x,y}"
sp "distavphi-map.H2O" u 2:1:3 
#
set auto
set pm3d clip4in
set origin 0.5,0.46
set size 0.52,0.4 
set xrange [-1.8:1.8]
set yrange [0:1.8]
set xlabel "" offset 0
set ylabel "momentum p_{x,y} [a@_0^{-1}]" offset 0
set ytics 0.5
unset key
set log cb
set form cb "10^{%T}"
set cblabel "PES yield [arb.u.]" offset 0.5
set cbrange [5e-10:4e-3]
set title "smoothed PES in plane of p_z and p_{x,y}"
sp "velomap.H2O" u 2:1:3 
#
set auto
set size 0.5,0.4
set origin 0.3,0.92
set xlabel "kinetic energy E_{kin} [eV]" offset 0
set ylabel "PES yield [arb.u.]" offset 0.5
set xrange [0:3*13.6]
set log y
set ytics 10
set yrange [5e-9:2e-1]
set form y "10^{%T}"
unset log cb
set label '{/=18 H_2O}' at graph 0.88,0.92
set label 'I=3{/Symbol= *}10^{14}W/cm^2' at graph 0.62,0.82
set label '{/Symbol= w}_{las}=11.4 eV, T_{pulse}=36 fs' at graph 0.62,0.74
set title "angular integrated PES along kinetic energy E_{kin}"
p "iPES.H2O" u ($1*13.6):2 w l lw 2 not
#
#!eps2pdf -f H2O-PESandmore.eps
!eps2jpg -f H2O-PESandmore.eps
!rm -f H2O-PESandmore.eps
#!evince H2O-PESandmore.pdf
!eog H2O-PESandmore.jpg
q
