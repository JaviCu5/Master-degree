file1 = "practica2_data.dat"
file2 = "practica2_data2.dat"
file3 = "practica2_E2.dat"
file4 = "practica2_E2Snapshots.dat"


set title font "Arial, 14"

set terminal png
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"
set xzeroaxis
set yzeroaxis
set key left
set grid

###########
#EXECICE 1#
###################################
#MSQ versus time for different δ/σ#
###################################
set xlabel "time s"
set ylabel "Δ^2(t) μm^2"

set logscale x
set logscale y

set output "MSD_time.png"
set style fill solid 0.2 border rgb 'grey30'
#plot file1 i 0 u 1:2 w l ls 1 lw 2 t "δ/σ = 0.001", file1 i 1 u 1:2 w l ls 2 lw 2 t "δ/σ = 0.003", file1 i 2 u 1:2 w l ls 3 lw 2 t "δ/σ = 0.01", file1 i 3 u 1:2 w l ls 4 lw 2 t "δ/σ = 0.03", file1 i 4 u 1:2 w l ls 5 lw 2 t "δ/σ = 0.1", file1 i 5 u 1:2 w l ls 6 lw 2 t "δ/σ = 0.3"


unset key

unset logscale x
unset logscale y

####################
#Diffusivity vs δ/σ#
####################
set xlabel "δ/σ"
set ylabel "D μm^2/s"

set output "D_delta.png"
#plot file2 i 0 u 1:2 w l ls 1


set key left

###########
#EXECICE 2#
###################################
#MSQ versus time for different δ/σ#
###################################
set xlabel "Time s"
set ylabel "Δ^2(t) μm^2"

set logscale x
set logscale y
set xrange[1:100000]
set yrange[0.000001:200]
set output "MSDTriangular.png"
plot file3 i 0 u 1:2 w l ls 1 lw 2 t "φ = 0.05", file3 i 1 u 1:2 w l ls 2 lw 2 t "φ = 0.2", file3 i 2 u 1:2 w l ls 3 lw 2 t "φ = 0.5"


unset key
unset title

###########
#Snapshots#
###########
unset xlabel
unset ylabel

unset logscale x
unset logscale y
set terminal png size 600,600
##############  005  ###############
set xrange[0:56.1]
set yrange[0:56.1]
set output "phi005t0.png"
plot file4 i 0 u 1:2 w p pt 6 ps 1.9

set output "phi005t2.png"
plot file4 i 2 u 1:2 w p pt 6 ps 1.9

set output "phi005t6.png"
plot file4 i 5 u 1:2 w p pt 6 ps 1.9

##############  02  ################
set xrange[0:28.1]
set yrange[0:28.1]
set output "phi02t0.png"
plot file4 i 6 u 1:2 w p pt 6 ps 3.3

set output "phi02t2.png"
plot file4 i 7 u 1:2 w p pt 6 ps 3.3

set output "phi02t5.png"
plot file4 i 10 u 1:2 w p pt 6 ps 3.3


##############  05  ################
set xrange[0:17.8]
set yrange[0:17.8]
set output "phi05t0.png"
plot file4 i 12 u 1:2 w p pt 6 ps 5.2

set output "phi05t1.png"
plot file4 i 13 u 1:2 w p pt 6 ps 5.2

set output "phi05t2.png"
plot file4 i 14 u 1:2 w p pt 6 ps 5.2

#w p pt 6 ps 2





