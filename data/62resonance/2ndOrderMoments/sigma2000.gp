set xlabel "cell"
set ylabel "<x^2>"
set grid

# This plots the big plot
plot "Ksc0_62resonance.data" using 1 w l title 'Ksc=0.0',\
     "Ksc1e8_62resonance.data" using 1 w l title 'Ksc=1e-8',\
     "Ksc1e9_62resonance.data" using 1 w l title 'Ksc=1e-9'
