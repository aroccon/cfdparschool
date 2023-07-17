#!/bin/bash
#SBATCH --account="tra23_CFD"
#SBATCH --job-name="cfdparschool"
#SBATCH --time=00:10:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=1  ##adjust
#SBATCH --gres=gpu:2   ###2 GPUs per node on G100
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=g100_usr_prod

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load nvhpc/22.3 

./application_name
