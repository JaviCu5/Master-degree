file1="Histograms.dat"

set title font "Arial, 14"

set terminal png
set ylabel font "Arial, 12"
set xlabel font "Arial, 12"
set xzeroaxis
set yzeroaxis
unset key
#################################
#Histograms for the default case#
#################################
f(x)=(1/(sqrt(2*pi)*12.87))*exp(-(x-100)**2/(2*12.87**2))
g(x)=(1/(sqrt(2*pi)*128.86))*exp(-(x-10000)**2/(2*128.86**2))

set xlabel "concentration (nM)"
set ylabel "probability density"
set output "mRNA_prob_dens.png"
set title "Probability density of mRNA concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 0 u 1:2:3 w boxerrorbars ls 5, f(x) ls 5

set output "Prot_prob_dens.png"
set title "Probability density of protein concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 0 u 4:5:6 w boxerrorbars ls 5, g(x) ls 5

##################################################################
#Histograms for the effect of transcription and translation rates#
##################################################################
f(x)=(1/(sqrt(2*pi)*40.75))*exp(-(x-1000)**2/(2*40.75**2))
g(x)=(1/(sqrt(2*pi)*128.86))*exp(-(x-10000)**2/(2*128.86**2))

set output "mRNA_prob_dens_rates.png"
set title "Probability density of mRNA concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 1 u 1:2:3 w boxerrorbars ls 5, f(x) ls 5

set output "Prot_prob_dens_rates.png"
set title "Probability density of protein concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 1 u 4:5:6 w boxerrorbars ls 5, g(x) ls 5

#############################################
#Histograms for the negative feedback effect#
#############################################
f(x)=(1/(sqrt(2*pi)*12.87))*exp(-(x-100)**2/(2*12.87**2))
g(x)=(1/(sqrt(2*pi)*128.86))*exp(-(x-10000)**2/(2*128.86**2))
set output "mRNA_prob_dens_negFeedback.png"
set title "Probability density of mRNA concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 2 u 1:2:3 w boxerrorbars ls 5, f(x) ls 5

set output "Prot_prob_dens_negFeedback.png"
set title "Probability density of protein concentration"
set style line 5 lc rgb 29390 pt 5 ps 1 lt 5 lw 1
set style fill solid 0.2 border rgb 'grey30'
plot file1 i 2 u 4:5:6 w boxerrorbars ls 5, g(x) ls 5