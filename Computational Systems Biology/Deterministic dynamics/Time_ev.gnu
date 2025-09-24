file1="Deterministic.dat"
file2="Gillespie.dat"

set title font "Arial, 14"

set terminal png
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"

##########################################
#Time evolution of the mRNA concentration#
##########################################
set xlabel "time (min)"
set ylabel "concentration (nM)"

set key bottom right
set grid

set output "evol_mRNA_00.png"
set title "Time evolution of the mRNA concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'

plot file1 index 0 u 1:2 w l t"Deterministic evolution", file2 index 0 u 1:2 w l ls 2 t"Stochastic evolution 1", file2 index 1 u 1:2 w l ls 3 t"Stochastic evolution 2"

##########################################
#Time evolution of the Prot concentration#
##########################################

set output "evol_prot_00.png"
set title "Time evolution of the protein concentration"
set style fill solid 0.2 border rgb 'grey30'
plot file1 index 0 u 1:3 w l t"Deterministic evolution", file2 index 0 u 1:3 w l ls 2 t"Stochastic evolution 1", file2 index 1 u 1:3 w l ls 3 t"Stochastic evolution 2"

##########################################
#Time evolution of the mRNA concentration#
##########################################
set output "evol_mRNA_rates.png"
set title "Time evolution of the mRNA concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'

plot file1 index 1 u 1:2 w l t"Deterministic evolution", file2 index 2 u 1:2 w l ls 2 t"Stochastic evolution 1", file2 index 3 u 1:2 w l ls 3 t"Stochastic evolution 2"

##########################################
#Time evolution of the Prot concentration#
##########################################

set output "evol_prot_rates.png"
set title "Time evolution of the protein concentration"
set style fill solid 0.2 border rgb 'grey30'
plot file1 index 1 u 1:3 w l t"Deterministic evolution", file2 index 2 u 1:3 w l ls 2 t"Stochastic evolution 1", file2 index 3 u 1:3 w l ls 3 t"Stochastic evolution 2"

##########################################
#Time evolution of the mRNA concentration#
##########################################
set key top right
set output "evol_mRNA_negFeedback.png"
set title "Time evolution of the mRNA concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'

plot file1 index 2 u 1:2 w l t"Deterministic evolution", file2 index 4 u 1:2 w l ls 2 t"Stochastic evolution 1", file2 index 5 u 1:2 w l ls 3 t"Stochastic evolution 2"

##########################################
#Time evolution of the Prot concentration#
##########################################
 set key bottom right
set output "evol_prot_negFeedback.png"
set title "Time evolution of the protein concentration"
set style fill solid 0.2 border rgb 'grey30'
plot file1 index 2 u 1:3 w l t"Deterministic evolution", file2 index 4 u 1:3 w l ls 2 t"Stochastic evolution 1", file2 index 5 u 1:3 w l ls 3 t"Stochastic evolution 2"