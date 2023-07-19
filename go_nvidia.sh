#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p g100_usr_interactive
#SBATCH -A tra23_cfd
#SBATCH -t 00:10:00 

module load nvhpc/22.3 

echo $HOSTNAME > hostname.dat
./application_name  > out.4096.C.$SLURM_JOBID.dat