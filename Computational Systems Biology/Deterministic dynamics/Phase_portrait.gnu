#set term png

file1="Deterministic.dat"
file2="Gillespie.dat"

#set output 'gene_reg.png'
set title "mRNA and protein concentrations evolution"
#unset multiplot
set ylabel "mRNA (nM)"
set xlabel "protein (nM)"
set yrange [0:150]
set xrange [0:12000]
# f: mrna dynamics factor
fm=1.0
fp=1.0
alpham=100.0*fm
deltam=1.0*fm
alphap=10.0*fp
deltap=0.1*fp
# based on axes definition above, x is the protein
mnull(x)=alpham/deltam # mrna-nullcline
pnull(x)=deltap*x/alphap # protein-nullcline
set samples 20 #x-axis
set isosamples 20 #y-axis
set key bmargin center horizontal
plot file1 index 0 u 3:2 w l t"Deterministic", file2 index 0 u 3:2 w l t"Stochastic",mnull(x),pnull(x),'++' u ($1*3000):($2*50):((alphap*($2*50)-deltap*($1*3000))/(2*fm/fp)):((alpham-deltam*($2*50))/(2*fm/fp)) with vectors head size 0.08,20,60 filled lc rgb "gray" t"Vector field"

#unset multiplot

pause -1