#!/usr/bin/env python3

import os
import qiskit
os.environ['OPENBLAS_NUM_THREADS'] = 'unlimited'
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
import cirq
from cirq.contrib.qasm_import import circuit_from_qasm
from cirq.qis import fidelity

# def initialize(NUMBER_OF_QUBITS_ARGS, subdir_arg, code_dir_arg, rand_circ_files_arg, one_qubit_error_space_arg, file_info_arg, file_number_arg, circ_arg,
#     max_pauli_weight_arg, max_pauli_str_p1_arg, max_pauli_str_p2_arg, quantum_register_arg, ancilla_register_arg, all_qubits_arg,
#     circ_full_arg, circ_no_checks_arg, final_state_no_errors_with_checks_arg, final_state_no_errors_no_checks_arg, reduced_correct_state_arg, file_name_arg,
#     circ_no_checks_tomo_arg, circ_full_no_ancilla_tomo_arg, circ_full_with_ancilla_tomo_arg, cpus_pool_arg, do_tomography_arg):
#     '''Initialize globals for pool'''
#     global NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files, one_qubit_error_space, file_info, file_number, circ
#     global max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, quantum_register, ancilla_register, all_qubits
#     global circ_full, circ_no_checks, final_state_no_errors_with_checks, final_state_no_errors_no_checks, reduced_correct_state, file_name
#     global circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo, cpus_pool#, simulator
#     NUMBER_OF_QUBITS=NUMBER_OF_QUBITS_ARGS
#     subdir= subdir_arg
#     code_dir=code_dir_arg
#     rand_circ_files= rand_circ_files_arg
#     one_qubit_error_space=one_qubit_error_space_arg
#     file_info=file_info_arg 
#     file_number=file_number_arg
#     max_pauli_weight=max_pauli_weight_arg 
#     max_pauli_str_p1=max_pauli_str_p1_arg 
#     max_pauli_str_p2=max_pauli_str_p2_arg
#     quantum_register=quantum_register_arg 
#     ancilla_register=ancilla_register_arg
#     all_qubits=all_qubits_arg
#     circ_full=circ_full_arg
#     circ_no_checks=circ_no_checks_arg
#     final_state_no_errors_with_checks=final_state_no_errors_with_checks_arg
#     final_state_no_errors_no_checks=final_state_no_errors_no_checks_arg
#     reduced_correct_state=reduced_correct_state_arg
#     circ=circ_arg
#     file_name=file_name_arg
#     circ_no_checks_tomo=circ_no_checks_tomo_arg
#     circ_full_no_ancilla_tomo=circ_full_no_ancilla_tomo_arg
#     circ_full_with_ancilla_tomo=circ_full_with_ancilla_tomo_arg
#     cpus_pool=cpus_pool_arg
#     do_tomography=do_tomography_arg
#     #simulator=simulator_arg

# def build_circ_no_intial_state(cirq_circ_p1, cirq_circ_compute, cirq_circ_p2, ancilla_name):
#     '''Builds circuit without initial state'''
#     ancilla_qubit=cirq.NamedQubit(ancilla_name+"_0")
#     print(cirq_circ_p1)
#     print(cirq_circ_p2)
#     circ=cirq_circ_compute+cirq_circ_p2
#     # circ.append([cirq.H.on(ancilla_qubit)])
#     # print("test\n", circ)

if __name__ == "__main__":
    print("running...")
    #Program parameters. These are globals in the parallel program passed from shell.
    # NUMBER_OF_QUBITS=int(sys.argv[1])
    NUMBER_OF_QUBITS=5#int(sys.argv[1])
    CNOT_COUNT=1#int(sys.argv[2])
    start_circ_number=0#int(sys.argv[3])
    end_circ_number=0#int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=False
    #Error space
    one_qubit_error_space=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    # one_qubit_error_space=np.logspace(-3, 0, num=21)
#     one_qubit_error_space=[1.00000000e-05, 1.41253754e-05, 1.99526231e-05, 2.81838293e-05,
#  3.98107171e-05, 5.62341325e-05, 7.94328235e-05, 1.12201845e-04,
#  1.58489319e-04, 2.23872114e-04, 3.16227766e-04, 4.46683592e-04,
#  6.30957344e-04, 8.91250938e-04, 1.25892541e-03, 1.77827941e-03,
#  2.51188643e-03, 3.54813389e-03, 5.01187234e-03, 7.07945784e-03,
#  1.00000000e-02, 0.015, 0.02, 0.025, 0.03, 0.035, 0.06]

    # File stuff
    subdir="/data/"
    # Gets the file path of the script
    # code_dir=sys.path[0]
    code_dir=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=code_dir+subdir
    # Gets the files that match the string. Files include the path
    files_found=mymodule.get_files(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT, start_circ_number, end_circ_number)
    rand_circ_files=files_found[0]
    circ_properties_files=files_found[1]

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

        ancilla_name="q1"
        created_circs=mymodule.create_circs(NUMBER_OF_QUBITS, circ_pieces, ancilla_name)
        # print(created_circs)
        cirq_circ_with_checks =circuit_from_qasm(created_circs.circ_with_checks.qasm())
        cirq_circ_no_checks=circuit_from_qasm(created_circs.circ_no_checks.qasm())
        cirq_circ_initial_state=circuit_from_qasm(created_circs.circ_initial_state.qasm())
        print(cirq_circ_initial_state)
        # print(created_circs.circ_initial_state)
        # print(cirq_circ_initial_state)
        # #Test
        # print("circ")
        # print(circ_from_qasm)
        # print("circ full")
        # print(circ_full)
        # print("circ no checks")
        # print(circ_no_checks)
        # print(cirq_circ_full)
        # print(cirq_circ_no_checks)

        keep_qubits=list(range(NUMBER_OF_QUBITS))
        # Cirq_circ_no_checks has NUMBER_OF_QUBITS+1 qubits including the ancilla. However

        rho_correct=mymodule.get_result_rho(cirq_circ_no_checks, NUMBER_OF_QUBITS, keep_qubits)
        rho_checks=mymodule.get_result_rho(cirq_circ_with_checks, NUMBER_OF_QUBITS+1, keep_qubits)
        sanity_check_fidelity=fidelity(rho_checks, rho_correct, qid_shape=(2**(NUMBER_OF_QUBITS),), validate=False)
        # print(np.trace(rho_checks))
        # print(np.trace(rho_correct))
        # cirq.validate_density_matrix(rho_correct, qid_shape=(2**(NUMBER_OF_QUBITS),))
        # cirq.validate_density_matrix(rho_checks, qid_shape=(2**(NUMBER_OF_QUBITS),))
        # sanity_check_fidelity=qiskit.quantum_info.state_fidelity(DensityMatrix(rho_checks), DensityMatrix(rho_correct))
        print("Sanity check fidelity: ", sanity_check_fidelity)
        assert sanity_check_fidelity>0.98, "Sanity check fidelity failed for circuit "+file_name

        #Noisy stuff.
        ancilla_qubit=cirq.NamedQubit(ancilla_name+"_0")
        # Creates a channel that applies the zero projector. We use this to get the measurement zero outcome of the
        # density matrix. Since the resulting trial density matrix is unormalized we can get the percentages of outcomes that
        # we discard.
        projector0_channel=cirq.KrausChannel(
            kraus_ops=(np.array([[1,0],[0,0]]),),
            validate=False
        )
        # cirq_circ_initial_state=circuit_from_qasm(created_circs.circ_initial_state.qasm())
        # circ_p1=circ_pieces[0]
        circ_compute=circ_pieces[1]
        # circ_p2=circ_pieces[2]
        # cirq_circ_with_checks_no_input=build_circ_no_intial_state(circ_p1, circ_compute, circ_p2, ancilla_name)
        #Error space
        single_qubit_error_space=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
        cirq_circ_with_checks.append([projector0_channel.on(ancilla_qubit)])
        results=[]
        test_circs=mymodule.TestCircuits(cirq_circ_with_checks, cirq_circ_no_checks, NUMBER_OF_QUBITS, rho_correct, sanity_check_fidelity, keep_qubits)
        for error_idx, single_qubit_error in enumerate(single_qubit_error_space):
            results.append(test_circs.run_test((error_idx, single_qubit_error)))

        # with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
        #     circ_info=pickle.load(circ_file)
        # mymodule.store_results(circ_compute, created_circs.circ_with_checks, BASE_PATH, file_name, NUMBER_OF_QUBITS, results, CNOT_COUNT, circ_info["rz"], circ_info["circuit_num"], 
        # circ_info["found_matches"], circ_info["max_pauli_weight"], circ_info["max_pauli_str_p1"], circ_info["max_pauli_str_p2"])

        print("execution time", time.time()-time0)
    print("Finished.")