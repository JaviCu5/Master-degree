file1 = "MSD.dat"
file2 = "Autocorr.dat"
file3 = "XDisp.dat"
file4 = "muDiff.dat"
file5 = "Trajectories.dat"

set terminal png

set title font "Arial, 16"
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"

unset key
unset title

set xlabel "X (µm)"
set ylabel "Y (µm)"

set output "traj1.png"
plot file5 i 0 u 1:2 w l lw 1 lc'blue'

set output "traj2.png"
plot file5 i 1 u 1:2 w l lw 1 lc'blue'

set output "traj3.png"
plot file5 i 2 u 1:2 w l lw 1 lc'blue'

set xrange[0.1:100]

###############
#Autocorr EX 1#
###############

set xlabel "time (s)"
set ylabel "Autocorrelation (µm²)"

#set terminal png size 1000,500
set output "Autocorr.png"
plot file2 i 0 u 1:3 w l lw 1 lc'light-blue', file2 i 0 u 1:2 w l lw 1 lc'dark-blue', \
file2 i 1 u 1:3 w l lw 1 lc'light-green', file2 i 1 u 1:2 w l lw 1 lc'dark-green',  \
file2 i 2 u 1:3 w l lw 1 lc'light-red', file2 i 2 u 1:2 w l lw 1 lc'dark-red'

##########
#MSD EX 1#
##########
unset key
set xlabel "time (s)"
set ylabel "MSD (µm²)"
set logscale x
set logscale y

set output "MSD.png"
plot file1 i 0 u 1:3 w l lw 1 lc'light-blue', file1 i 0 u 1:2 w l lw 1 lc'dark-blue',  \
file1 i 1 u 1:3 w l lw 1 lc'light-green', file1 i 1 u 1:2 w l lw 1 lc'dark-green', \
file1 i 2 u 1:3 w l lw 1 lc'light-red', file1 i 2 u 1:2 w l lw 1 lc'dark-red'

###########
#Disp EX 1#
###########
set key inside top left
set xlabel "time (s)"
set ylabel "X displacement (µm)"

set output "disp.png"
plot file3 i 0 u 1:3 w l lt 1 lw 1 t'f=0', \
file3 i 1 u 1:3 w l lw 1 lt 2 lc 4 t'f=0.001', \
file3 i 2 u 1:3 w l lc 1 lw 1 dt 1 t'f=0.01', \
file3 i 3 u 1:3 w l lc 2 lw 1 dt 1 t'f=0.1', \
file3 i 4 u 1:3 w l lc 3 lw 1 dt 1 t'f=1'




set key inside top left
set xlabel "time (s)"
set ylabel "X MSD (µm²)"

set output "XMSD.png"
plot file3 i 0 u 1:2 w l lt 1 lw 1 t'f=0', \
file3 i 1 u 1:2 w l lw 1 lt 2 lc 4 t'f=0.001', \
file3 i 2 u 1:2 w l lc 1 lw 1 dt 1 t'f=0.01', \
file3 i 3 u 1:2 w l lc 2 lw 1 dt 1 t'f=0.1', \
file3 i 4 u 1:2 w l lc 3 lw 1 dt 1 t'f=1'



unset key
unset title
set xlabel "τ (s)"
unset ylabel
unset logscale x
unset logscale y

set ylabel "Mu (µm/N·s)"
set output "mu.png"
plot file4 i 0 u 1:2 w l lt 1 lw 1

set ylabel "D (µm²/s)"
set output "Diff.png"
plot file4 i 0 u 1:3 w l lw 1 lt 2 lc 4

