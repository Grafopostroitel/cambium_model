#
rm a.out
rm *.dat
rm *.plt

# gfortran -Wall -Wextra -Wconversion -pedantic -fbacktrace -fbounds-check -fdefault-real-8 -g -O0 growth_camb7.1.f95
# http://faculty.washington.edu/rjl/uwamath583s11/sphinx/notes/html/gfortran_flags.html
# http://linux.die.net/man/1/gfortran
# https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html
gfortran -Wall -Wextra -pedantic -fbacktrace -fbounds-check  -g -O0 growth_camb7.1.f95
./a.out
gnuplot rate_from_position.sh
gnuplot dynamic_cambium_cells.plt
gnuplot age_of_cells.plt
gnuplot position_of_cells.plt
gnuplot inhibitor_in_cells.plt
gnuplot rate_on_concentration.plt
