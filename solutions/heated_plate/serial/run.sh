rm heat

gfortran -fopenmp -mcmodel=medium -O3 -Wall -ffixed-line-length-140 ../module.f90 heat.f90 -o heat

./heat
