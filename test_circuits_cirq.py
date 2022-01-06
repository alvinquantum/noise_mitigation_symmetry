#!/usr/bin/env python3

import os, qiskit, cirq
from qiskit.ignis.verification import tomography
from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, Aer, transpile
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector, state_fidelity, DensityMatrix, Operator, Pauli
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit.visualization import circuit_drawer
from qiskit.circuit import QuantumRegister
import pickle, mymodule, json, psutil, time
import qiskit.providers.aer.noise as noise
from qiskit.opflow import X,Y,Z
from qiskit.ignis.verification.tomography import StateTomographyFitter
from math import e, pi
from multiprocessing import Pool, Value, Array, Manager
from os import listdir
from os.path import isfile
# import cProfile, io
# from pstats import SortKey
from cirq.contrib.qasm_import import circuit_from_qasm
from cirq.qis import fidelity

if __name__ == "__main__":
    print("running...")
    #Program parameters. These are globals in the parallel program passed from shell.
    NUMBER_OF_QUBITS=5#int(sys.argv[1])
    CNOT_COUNT=10#int(sys.argv[2])
    START_CIRC_NUMBER=0#int(sys.argv[3])
    END_CIRC_NUMBER=10#int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=False
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
#     SINGLE_QUBIT_ERROR_SPACE=[1.00000000e-05, 1.41253754e-05, 1.99526231e-05, 2.81838293e-05,
#  3.98107171e-05, 5.62341325e-05, 7.94328235e-05, 1.12201845e-04,
#  1.58489319e-04, 2.23872114e-04, 3.16227766e-04, 4.46683592e-04,
#  6.30957344e-04, 8.91250938e-04, 1.25892541e-03, 1.77827941e-03,
#  2.51188643e-03, 3.54813389e-03, 5.01187234e-03, 7.07945784e-03,
#  1.00000000e-02, 0.015, 0.02, 0.025, 0.03, 0.035, 0.06]

    # File stuff
    SUBDIR="/data/"
    # Gets the file path of the script
    # code_dir=sys.path[0]
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=CODE_DIR+SUBDIR
    # Gets the files that match the string. Files include the path
    rand_circ_files, circ_properties_files=mymodule.get_files(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT, START_CIRC_NUMBER, END_CIRC_NUMBER)

    #Get the qasm and pickle info
    for idx, file_name in enumerate(circ_properties_files):
        #Stopwatch
        time0=time.time()
        # os.environ['QISKIT_IN_PARALLEL'] = 'FALSE'

        #get the file information of the circuit. Len 3 list: qubits#, cnot#, circuit_number#.
        file_info=[int(elem) for elem in file_name.split("_") if elem.isdigit()]
        assert len(file_info) == 3, "This file doesn't have three numbers in the name."
        file_number=file_info[2]
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        if mymodule.result_exists(BASE_PATH, file_name, 0):
            continue
        
        checks_properties=mymodule.get_checks_properties(BASE_PATH, file_name)
        circ_pieces=mymodule.split_circuit_by_barrier(os.path.join(BASE_PATH, rand_circ_files[idx]))
        if len(circ_pieces)!=3:
            continue
        # If there were no Pauli checks found just skip.
        if checks_properties.p2_weights == 0: 
            continue

        created_circs=mymodule.create_circs(NUMBER_OF_QUBITS, circ_pieces)
        # print(created_circs)
        cirq_circ_with_checks =created_circs.circ_with_checks
        cirq_circ_no_checks=created_circs.circ_no_checks
        cirq_circ_initial_state=created_circs.circ_initial_state
        # print(cirq_circ_initial_state)

        keep_qubits=list(range(NUMBER_OF_QUBITS))

        rho_correct=mymodule.get_result_rho(cirq_circ_no_checks, NUMBER_OF_QUBITS, keep_qubits)
        rho_checks=mymodule.get_result_rho(cirq_circ_with_checks, NUMBER_OF_QUBITS+1, keep_qubits)
        sanity_check_fidelity=fidelity(rho_checks, rho_correct, qid_shape=(2**(NUMBER_OF_QUBITS),), validate=False)
        print("Sanity check fidelity: ", sanity_check_fidelity)
        assert sanity_check_fidelity>0.98, "Sanity check fidelity failed for circuit "+file_name

        #Noisy stuff.
        cirq_circ_compute=circ_pieces[1]
        results=[]
        test_circs=mymodule.TestCircuits(cirq_circ_with_checks, cirq_circ_no_checks, NUMBER_OF_QUBITS, rho_correct, sanity_check_fidelity, keep_qubits)
        for error_idx, single_qubit_error in enumerate(SINGLE_QUBIT_ERROR_SPACE):
            results.append(test_circs.run_test((error_idx, single_qubit_error)))

        # with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
        #     circ_info=pickle.load(circ_file)
        # mymodule.store_results(circ_compute, created_circs.circ_with_checks, BASE_PATH, file_name, NUMBER_OF_QUBITS, results, CNOT_COUNT, circ_info["rz"], circ_info["circuit_num"], 
        # circ_info["found_matches"], circ_info["max_pauli_weight"], circ_info["max_pauli_str_p1"], circ_info["max_pauli_str_p2"])

        print("execution time", time.time()-time0)
    print("Finished.")