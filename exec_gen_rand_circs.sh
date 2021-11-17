#!/bin/bash
# For submitting multiple jobs of the gen_rand_circs.sh
for idx in 35 40 
do
    sbatch gen_rand_circs.sh 5 $idx 20
done