#!/bin/bash

#SBATCH --job-name=genrandcircs
#SBATCH --account=quantumsimulations
#SBATCH --partition=bdwall
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=36
#SBATCH --output=genrandcircs.out
#SBATCH --error=genrandcircs.error
##SBATCH --mail-user=<your email address> # Optional if you require email
##SBATCH --mail-type=ALL                  # Optional if you require email
#SBATCH --time=01:00:00

# Run My Program
srun python gen_rand_circs_hpc.py