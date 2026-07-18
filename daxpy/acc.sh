module load nvhcp/24.5
nvfortran -fast -acc -Minfo=accel -gpu=managed daxpy_acc.f90  -o daxpy
