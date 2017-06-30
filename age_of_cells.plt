 # gnuplot cambium graph
 reset
 set terminal png size 1024,768 
 set output "age_of_cells.dat.png"
 #color definitions
 set border linewidth 1.5
 set style line 1 lc rgb "#0060ad" lt 1 lw 2 # --- blue
 unset key
 # Axes
 set style line 11 lc rgb "#808080" lt 1
 set border 3 back ls 11
 set tics nomirror out scale 0.75
 # Grid
 set style line 12 lc rgb "#808080" lt 0 lw 1
 set grid back ls 12
 set title "Graph the age of cells (win 3)"
 set xlabel "Cell"
 set ylabel "Age"
 plot  "age_of_cells.dat" w l ls 1
