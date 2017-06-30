rem http://www.math.uni-leipzig.de/~hellmund/Vorlesung/gfortran.html

del a.exe
del *.dat
del *.plt
del *.plt.png
del *.mod0
del *.mod


rem gfortran -Wall -Wextra -pedantic -fbacktrace -fbounds-check -fcheck=all -g -O0 -static -static-libgfortran growth_camb7.1.f95

gfortran -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -fbounds-check -fcheck=all -Wuninitialized -g -O0 growth_camb7.1.f95

a.exe
gnuplot rate_from_position.sh
gnuplot dynamic_cambium_cells.plt
gnuplot age_of_cells.plt
gnuplot position_of_cells.plt
gnuplot inhibitor_in_cells.plt
gnuplot rate_on_concentration.plt
