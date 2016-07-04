set xlabel "cell"
set ylabel "<p_x^2>"
set grid

# This plots the big plot
plot "Ksc0_62resonance.data" using 2 w l title 'Ksc=0.0',\
     "Ksc1e8_62resonance.data" using 2 w l title 'Ksc=1e-8',\
     "Ksc1e9_62resonance.data" using 2 w l title 'Ksc=1e-9'
