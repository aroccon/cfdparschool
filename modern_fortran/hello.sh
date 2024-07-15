 module load intel-oneapi-compilers/2023.2.1
 module load intel-oneapi-mpi/2021.10.0     

 export FOR_COARRAY_NUM_IMAGES=2
 ifort -coarray hello.f90 -o hello
 ./hello