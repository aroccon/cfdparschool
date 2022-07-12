rm check

gfortran -mcmodel=medium -O3 -Wall -ffixed-line-length-140 ../module.f90 check.f90 -o check

./check
