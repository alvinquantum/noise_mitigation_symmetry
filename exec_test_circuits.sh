#!/bin/bash
# For submitting multiple jobs of the test_circuits.sh
for idx in 1 5 10 15 20 25 30 35 40
do
    sbatch test_circuits.sh 5 $idx 0 4
    sbatch test_circuits.sh 5 $idx 5 9
    sbatch test_circuits.sh 5 $idx 10 14
    sbatch test_circuits.sh 5 $idx 15 19
done