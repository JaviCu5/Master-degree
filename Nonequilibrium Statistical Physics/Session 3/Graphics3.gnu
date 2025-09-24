file1="IniMatrix.dat"
file2="FSMCMagne.dat"
file3="CTMCMagne.dat"

set terminal png
unset title
set title font "Arial, 16"
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"

###########
#HEAT MAPS#
###########
set title "Stocastic initial configuration for a spin matrix"
set xlabel "j"
set ylabel "i"
unset key
set border linewidth 0
set tics

Lx = 50
Ly = 50
set palette defined (-1 "pink", 1 "red")
set xrange[0.5:Lx+0.5]
set yrange[0.5:Ly+0.5]
set output "IniMatrix.png"
plot file1 u 1:2:3 w image

##############################
#MAGNETIZATION MetropolisFSMC#
##############################

unset title
unset xlabel
unset ylabel
set logscale x
set logscale y

set xzeroaxis
unset xrange
set yrange[0.01:1.2]
set xrange[1:1000]
set key top left
set output "magne1.png"
plot file2 i 0 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 0 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne2.png"
plot file2 i 1 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 1 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne3.png"
plot file2 i 2 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 2 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne4.png"
plot file2 i 3 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 3 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne5.png"
plot file2 i 4 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 4 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne6.png"
plot file2 i 5 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 5 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne7.png"
plot file2 i 6 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 6 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne8.png"
plot file2 i 7 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 7 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne9.png"
plot file2 i 8 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 8 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "magne10.png"
plot file2 i 9 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 9 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'

set output "MeanMagnet.png"
plot file2 i 10 u 1:2 w l lw 1 lc rgb 0255127 t'FSMC', file3 i 10 u 1:2 w l lw 1 lc rgb 16737792 t'CTMC'


