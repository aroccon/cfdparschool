module load nvhpc/25.11
nvfortran  -fast -acc -gpu=managed  -Minfo=accel miniweather_acc.f90 -o miniweather_acc
