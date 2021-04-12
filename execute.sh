gfortran -Wall -c Parameters.f90
gfortran -Wall -c Distributions.f90
gfortran -Wall -c Functions.f90
gfortran -Wall -c radiative.f90
gfortran -Wall -c main.f90
gfortran -Wall -o main Parameters.o Functions.o Distributions.o radiative.o main.o
./main
gnuplot plot_sed.gp

