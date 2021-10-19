#!/bin/bash

#SBATCH --job-name=test_circuit
#SBATCH --account=quantumsimulations
#SBATCH --partition=knlall
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=36
#SBATCH --output=test_circuit.out
#SBATCH --error=test_circuit.error
##SBATCH --mail-user=<your email address> # Optional if you require email
##SBATCH --mail-type=ALL                  # Optional if you require email
#SBATCH --time=72:00:00

# Run My Program
# pass the number qubits and depth
srun python test_circuit.py 5 25