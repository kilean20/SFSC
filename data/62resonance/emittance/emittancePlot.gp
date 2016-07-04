set xlabel "cell"
set ylabel "emittance"
set grid

# This plots the big plot
plot "emittance_Ksc0.data" using 1 w l title 'Ksc=0.0',\
