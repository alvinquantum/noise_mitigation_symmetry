#!/bin/bash

#SBATCH --job-name=test_circuit
#SBATCH --account=quantumsimulations
#SBATCH --partition=bdwall
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=36
#SBATCH --output=test_circuit.out
#SBATCH --error=genrandcircs.error
##SBATCH --mail-user=<your email address> # Optional if you require email
##SBATCH --mail-type=ALL                  # Optional if you require email
#SBATCH --time=01:00:00

# Run My Program
srun python test_circuit.sh