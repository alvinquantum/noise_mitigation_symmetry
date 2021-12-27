#!/bin/bash

#SBATCH --job-name=test_circuit
#SBATCH --account=quantumsimulations
#SBATCH --partition=bdwall
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=21
#SBATCH --output=test_circuit.%j.out
#SBATCH --error=test_circuit.%j.error
##SBATCH --mail-user=<your email address> # Optional if you require email
##SBATCH --mail-type=ALL                  # Optional if you require email
#SBATCH --time=72:00:00

# Run My Program
# pass the number qubits, cnot, and range
# srun python -m cProfile -o test_circuits.profile test_circuits.py $1 $2 $3 $4
# srun kernprof -l -v test_circuits.py $1 $2 $3 $4
srun python test_circuits.py $1 $2 $3 $4
