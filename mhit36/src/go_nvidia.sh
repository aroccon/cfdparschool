#!/bin/bash
#SBATCH --account="s_tra_cfd_*"
#SBATCH --job-name="toy36"
#SBATCH --time=00:05:00
#SBATCH --nodes=1     
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1   
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH -p boost_usr_prod

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module load nvhpc/23.11
module load cuda/12.1
module list

#if using HPC-SDK (OPENPMPI) use (CUDA-aware already enabled):
./mhit36
