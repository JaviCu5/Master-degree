file1="Act_Kinases_concentrationsVSTime.dat"
file2="Fig1_stationaryConc_vs_EC50.dat"
file3="Stationary_concentrationVSE1_log.dat"

set title font "Arial, 14"

set terminal png
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"

unset key
set grid

############################################
set xlabel "Time (s)"
set ylabel "Concentration (nM)"

set output "KKKP_2.png"
set title "MAPKKK-P concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 0 u 1:2 w l ls 5

############################################

set output "KKPP_2.png"
set title "MAPKK-PP concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 0 u 1:3 w l ls 5

############################################

set output "KPP_2.png"
set title "MAPK-PP concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 0 u 1:4 w l ls 5

############################################
############################################

set key inside center right
set xlabel "E1_{tot} in multiples of the EC50"
set ylabel "Predicted steady-state kinease activity"
set output "AllMapkvsE1.png"
set title "Stationary concentration vs E1 initial concentration"
set xrange [0:5]
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file2 i 0 u 1:2 w l ls 2 t'MAPKKK-P', file2 i 0 u 3:4 w l ls 3 t'MAPKK-PP',file2 i 0 u 5:6 w l ls 5 t'MAPK-PP'

############################################
############################################

set xlabel "E1 concentration a t=0 (μM)"
set xrange [0.000001:0.1]
set logscale x
set mxtics 10
set format x "%1.0tx10^{%S}" 

set output "AllMapkvsE1_2_logs.png"
set title "Stationary concentration vs E1 initial concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file3 i 0 u 1:2 w l ls 2 t'MAPKKK-P', file3 i 0 u 3:4 w l ls 3 t'MAPKK-PP',file3 i 0 u 5:6 w l ls 5 t'MAPK-PP', 0.1 lc "black" notitle, 0.9 notitle lc "black"