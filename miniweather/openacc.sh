module load nvhpc/24.3
nvfortran  -fast -acc -gpu=managed  -Minfo=accel miniweather_acc.f90 -o miniweather_acc
