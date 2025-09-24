file1 = "practica2_E5.dat"
file2 = "practica2_E5Snap.dat"
file3 = "rhos.dat"

set title font "Arial, 14"

set terminal png
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"
set xzeroaxis
set yzeroaxis
set key left
set grid

set key right

unset logscale x
unset logscale y

set xlabel "Y (µm)"
set ylabel "ρ(Y)"

set yrange[-0.02:0.14]
set xtics 10
set output "densitatY.png"

plot file3 i 0 u 1:2 w l lc 1 lw 2 t'g = 0',  file3 i 1 u 1:2 w l lc 2 lw 2 t'g = 0.01', file3 i 2 u 1:2 w l lc 3 lw 2 t'g = 0.1', file3 i 3 u 1:2 w l lc 4 lw 2 t'g = 1', file3 i 4 u 1:2 w l lc 5 lw 2 t'g = 10'

set key left

###########
#EXECICE 5#
#################################
#MSQ versus time for different g#
#################################
set xlabel "Time s"
set ylabel "Δ^2(t) μm^2"

set logscale x
set logscale y
set xrange[1:10000]
set yrange[0.01:700]

set output "MSDGrav.png"
plot file1 i 0 u 1:2 w l ls 1 lw 2 t "g = 0", file1 i 1 u 1:2 w l ls 2 lw 2 t "g = 0.01", file1 i 2 u 1:2 w l ls 3 lw 2 t "g = 0.1", file1 i 3 u 1:2 w l ls 4 lw 2 t "g = 1", file1 i 4 u 1:2 w l ls 3 lw 2 t "g = 10"

unset key
unset title

set xrange[0:7.3]
set yrange[0:72.5]

set xtics 3

###########
#Snapshots#
###########
unset xlabel
unset ylabel

unset logscale x
unset logscale y

set terminal png size 250,1000
##############  0  ###############
set output "Grav1_0.png"
plot file2 i 0 u 1:2 w p pt 6 ps 2.5


##############  001  ################
set output "Grav2_001.png"
plot file2 i 1 u 1:2 w p pt 6 ps 2.5


##############  01  ################
set output "Grav3_01.png"
plot file2 i 2 u 1:2 w p pt 6 ps 2.5


##############  1  ################
set output "Grav4_1.png"
plot file2 i 3 u 1:2 w p pt 6 ps 2.5

##############  10  ################
set output "Grav5_10.png"
plot file2 i 4 u 1:2 w p pt 6 ps 2.5





