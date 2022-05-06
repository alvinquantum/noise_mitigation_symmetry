#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from tkinter import N
# The number of threads is capped on Bebop.
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit import QuantumCircuit, QuantumRegister
import circgenerator
import checksfinder_new as checksfinder
import circtester
import sys
import time
import numpy as np
# from qiskit.circuit.random import random_circuit

if __name__ == "__main__":
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    NUMBER_OF_CIRCUITS=int(sys.argv[3])
    RZ_COUNT=int(sys.argv[4])

    time0=time.time()
    #Disable qiskit parallel.
    # os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=f"qubits_{NUMBER_OF_QUBITS}"
    CHECKS_SUBDIR="checks"
    RAWS_SUBDIR="raws"
    CHECKS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, RAWS_SUBDIR)

    # Generate random circuits
    file_number=0
    for _ in range(NUMBER_OF_CIRCUITS):
        time1=time.time()
        # circ=generate_a_random_circuit(number_of_qubits, cnot_count) 
        if RZ_COUNT>0:
            circ=circgenerator.generate_a_random_circuit(NUMBER_OF_QUBITS, CNOT_COUNT, RZ_COUNT)             
        else:
            circ=circgenerator.generate_a_random_clifford_circuit(NUMBER_OF_QUBITS, CNOT_COUNT)  
        #Save raw circuit
        circgenerator.write_circs_no_checks_to_qasm_file(RAWS_PATH, [circ], NUMBER_OF_QUBITS, CNOT_COUNT)

        # print(circ)
        circ_operations=circ.count_ops()
        print(circ_operations)
        
        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")