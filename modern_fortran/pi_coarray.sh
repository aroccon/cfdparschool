 module load intel-oneapi-compilers/2023.2.1
 module load intel-oneapi-mpi/2021.10.0     

 export FOR_COARRAY_NUM_IMAGES=2
 ifort -coarray pi_coarray.f90 -o pi
 ./pi