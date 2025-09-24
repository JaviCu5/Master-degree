file1 = "Task2.dat"
file2 = "FPHistograms.dat"


set title font "Arial, 14"

set terminal png
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"
set xzeroaxis
set yzeroaxis



############################################
#Normal distribution with parameters N(0,1)#
############################################
set xlabel "x"
set ylabel "p(x)"
f(x)=(1/(pi*(1+x^2)))

set logscale y

set output "NormalDistr.png"
set title "Normal distribution with parameters N(0,1)"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 1 u 1:2:3 w boxerrorbars ls 5 t"Experimental","" i 0 u 1:2 smooth freq with boxes notitle, f(x) ls 5 t"Exact"

unset logscale y
unset key



###########################
#Random walk of a particle#
###########################
set xlabel "x (μm)"
set ylabel "y (μm)"

set output "random_walk.png"
set title "Random walk"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30' 
plot file1 i 1 u 1:2 w lp ls 6 lt rgb 29390 notitle



################################
#Random walk for 1000 particles#
################################
set xlabel "X Coord"
set ylabel "Y Coord"
set output "random_walk_PBC.png"
set multiplot layout 1,2 title "Random walk for N=1000 particles with PBC"


set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30' 

set title "Initial configuration"
plot file1 i 2 u 1:2 pt 7 ps 0.25 lt rgb 29390 notitle

unset ylabel
set title "Final configuration (τ=100δt)"
plot file1 i 2 u 3:4 pt 7 ps 0.25 lt rgb 29390 notitle

unset multiplot



#################
#Dffusivity vs Γ#
#################
set xlabel "Γ (μN^{2}/s)"
set ylabel "D (μm^{2}/s)"
set key inside
set logscale x
set logscale y

set output "DiffusivityNParticles.png"
set title "Diffusivity behaviour depending on Γ"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30' 
plot file1 i 3 u 2:1 title 'D vs Γ' w lp ls 6 lt rgb 'blue30'

unset logscale x
unset logscale y



####################
#Dffusion equiation#
####################
set xlabel "Δx (μm)"
set ylabel "ρ(Δx)"
set xrange [-100:100]
set output "Solution_diff_eq.png"
set title "Horizontal displacement distribution for different t"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file2 i 0 u 1:2 w l ls 1 lw 2 t"t = 1", file2 i 1 u 1:2 w l ls 2 lw 2 t"t = 10", file2 i 2 u 1:2 w l ls 3 lw 2 t"t = 30", file2 i 3 u 1:2 w l ls 4 lw 2 t"t = 100", file2 i 4 u 1:2 w l ls 5 lw 2 t"t = 1000"











