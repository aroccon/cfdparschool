rm heat


gfortran -fopenmp -mcmodel=medium -O3 -Wall -ffixed-line-length-140 ../module.f90 heat.f90 -o heat
# for Mac, increase stack size
#gfortran -fopenmp -mcmodel=medium -O3 -Wall -ffixed-line-length-140 -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 ../module.f90 heat.f90 -o heat


#set OMP_NUM_THREADS=4
export OMP_NUM_THREADS=4

# for Linux machines
#ulimit -s unlimited

./heat
