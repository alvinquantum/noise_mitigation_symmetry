#!/bin/bash
# For submitting multiple jobs of the test_circuits.sh
for idx in 1 5 15 20
do
    sbatch test_circuits.sh 10 $idx 0 10
    # sbatch test_circuits.sh 10 $idx 11 20
    # sbatch test_circuits.sh 10 $idx 21 30
    # sbatch test_circuits.sh 10 $idx 31 40
    # sbatch test_circuits.sh 10 $idx 41 49
done