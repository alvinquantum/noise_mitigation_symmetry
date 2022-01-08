#!/bin/bash
# For submitting multiple jobs of the gen_rand_circs.sh
for idx in 80 #1 5 10 15 20 25 30 35 40 
do
    # number of qubits, depth, number of circuits
    sbatch gen_rand_circs.sh 10 $idx 50
done