rm drop output/* vtk/*

gfortran -fopenmp -Wall -ffixed-line-length-140 module.f90 drop.f90 -o drop

./drop

gnuplot gif.plt
