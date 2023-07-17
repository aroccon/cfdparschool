#!/bin/bash
#SBATCH --account="tra23_CFD"
#SBATCH --job-name="cfdparschool"
#SBATCH --time=00:01:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4  ##adjust
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=g100_usr_prod

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load intel-oneapi-compilers/2021.4.0 
module load intel-oneapi-mpi/2021.4.0 

mpiexec -np XXXXX ./a.out
