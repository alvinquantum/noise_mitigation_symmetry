#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

from qiskit.converters.circuit_to_dag import circuit_to_dag
from qiskit.converters.dag_to_circuit import dag_to_circuit
from scipy.sparse.extract import find
# The number of threads is capped on Bebop.
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from re import M
import numpy as np
import math, sys, pickle, ctypes, psutil, json
from qiskit import QuantumCircuit, transpile
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector, Pauli, Operator
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
import mymodule
from qiskit.visualization import circuit_drawer
from multiprocessing import Pool, Value, Array, Manager
from copy import deepcopy
from pprint import pprint
import time
# from qiskit.circuit import ControlledGate
# from qiskit.extensions import UnitaryGate
import cProfile

def initialize(ABS_TOL_ARG, circ_properties_arg,
                    pauli_labels_arg, pauli_group_positive_arg, table_length_arg):
    '''Initialize globals for pool'''
    global ABS_TOL, circ_properties, pauli_labels, pauli_group_positive, table_length
    #Passed from the main
    #Program parameters passed from main
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    circ_properties=circ_properties_arg
    ABS_TOL=ABS_TOL_ARG

    #Tables passed from main
    pauli_labels=pauli_labels_arg
    pauli_group_positive=pauli_group_positive_arg
    table_length=table_length_arg


def get_checks_parallel(pauli_group_elem):
    '''Wrapper function. Get the checks in parallel.'''
    global ABS_TOL, circ_properties, pauli_labels, pauli_group_positive, table_length
    return mymodule.find_checks_sym(pauli_group_elem, circ_properties.circ)

if __name__ == "__main__":
    # #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    # func_params=FunctionParams()
    NUMBER_OF_QUBITS=10#int(sys.argv[1])
    CNOT_COUNT=10#int(sys.argv[2])
    NUMBER_OF_CIRCUITS=10#int(sys.argv[3])
    if NUMBER_OF_QUBITS<=5:
        PARALLEL=False
    else:
        PARALLEL=True
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    ABS_TOL=.2*2**(NUMBER_OF_QUBITS-1)
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    # code_dir1=os.path.dirname(os.path.realpath("__file__"))
    CODE_DIR=sys.path[0]
    # print(code_dir1)
    # print(code_dir2)
    SUBDIR="/data/"
    # BASE_FILE_PATH=code_dir+SUBDIR+"depth"+ str(DEPTH) +"_"
    BASE_FILE_PATH=CODE_DIR+SUBDIR+"qubits_"+str(NUMBER_OF_QUBITS)+"_CNOTS_"+ str(CNOT_COUNT) +"_"

    #Create +1 phase pauli group
    #Note that we can restrict the search for p1 to the +/-1 phases. This is due to passing the phase
    #from p2 to p1 and then realizing that the eigenvalues of p2 is restricted to +/-1. The unitary
    #conjugating p1 only rotates the eigenvectors. Therefore, the eigenvalues of p1 must also be +/-1
    #which eliminates the need to check for +/-j phases. This is explained in the paper.
    pauli_table=pauli_basis(NUMBER_OF_QUBITS)
    # pauli_group_positive=deepcopy(pauli_table)
    #Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()
    pauli_group=pauli_table
    pauli_group_positive=None
    
    # #Convert to matrices
    # pauli_group_positive=list(map(lambda x:x.to_matrix(),pauli_table))
    # #Flatten the list
    # pauli_group_positive=[val for sublist in pauli_group_positive for val in sublist]
    # # Make the other phases and merge
    # pauli_group_negative=[element * -1 for element in pauli_group_positive]
    # pauli_group=pauli_group_positive+pauli_group_negative

    table_length=len(pauli_table)
    #main loop
    file_number=0
    for _ in range(NUMBER_OF_CIRCUITS):
        time0=time.time()
        output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
        while os.path.isfile(output_file_path):
            file_number+=1
            output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
        file_info_path=BASE_FILE_PATH+ "circuit_" + str(file_number)+"_.obj"
        file_qasm_path=BASE_FILE_PATH+ "circuit_" + str(file_number)+"_.qasm"
        #Random circuit.
        circ=mymodule.random_circuit_cnot(NUMBER_OF_QUBITS, CNOT_COUNT)
        # circ=mymodule.random_circuit_depth(NUMBER_OF_QUBITS, CNOT_COUNT)
        # # # Load circuit.
        # circ_file=open(BASE_FILE_PATH+ "circuit_0_.obj", "rb")
        # circ=pickle.load(circ_file)["circ"]
        # circ.x(0)
        # circ_file.close()
        # Don't delete! Qiskit random circuit
        # circ=random_circuit(5,10)

        # #Test
        # circ=QuantumCircuit(5)
        # for idx in range(5):
        #     circ.s(idx)

        # We use this to find P1 and P2.
        unitary = Operator(circ)
        #Draw to file
        circuit_drawer(circ, filename=output_file_path)
        print(circ)
        
        output_file=open(output_file_path, "a")
        output_file.write("\n")
        # Count the number of CNOT gates
        circ_operations=circ.count_ops()
        output_file.write(json.dumps(circ_operations))
        print(circ_operations)
        # max_pauli_str_p2=Array(ctypes.c_char, "+1"+"I"*NUMBER_OF_QUBITS)

        #Results
        count=0
        p2_weights=[0]
        pauli_str_p1s=[""]
        pauli_str_p2s=[""]

        circ_properties=mymodule.CircuitProperties(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS, circ, circ_operations)
        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        if PARALLEL:
            with Pool(psutil.cpu_count(logical=False), initialize, initargs=(ABS_TOL, circ_properties,
                    pauli_labels, pauli_group_positive, table_length)) as pool:
                #In some cases pool.imap_unordered needs to be wrapped in list in order to return properly. 
                #see: https://stackoverflow.com/questions/5481104/multiprocessing-pool-imap-broken
                for result in pool.imap_unordered(get_checks_parallel, pauli_labels, chunksize=1000):
                # pool.map(find_p1s_p2s, enumerate(pauli_group))                    
                    #Store the results.
                    if result:
                        count+=1
                        if result[0]>p2_weights[0]:
                            p2_weights[0]=result[0]
                            pauli_str_p1s[0]=result[1]
                            pauli_str_p2s[0]=result[2]
                        p2_weights.append(result[0])
                        pauli_str_p1s.append(result[1])
                        pauli_str_p2s.append(result[2])
                        # Terminate the pool when found is true, i.e., found=1.
                        if result[0]==NUMBER_OF_QUBITS:
                            print("terminating...")
                            pool.close() #Stops passing jobs to processes.
                            pool.join() #Waits for processes to finish.
                            break
        else:
            for elem in pauli_labels:
                # temp_p2_circ=find_p2s(elem, circ)
                result=mymodule.find_checks_sym(elem, circ)
                if result:
                    count+=1
                    if result[0]>p2_weights[0]:
                        p2_weights[0]=result[0]
                        pauli_str_p1s[0]=result[1]
                        pauli_str_p2s[0]=result[2]
                    p2_weights.append(result[0])
                    pauli_str_p1s.append(result[1])
                    pauli_str_p2s.append(result[2])
                    # Terminate the pool when found is true, i.e., found=1.
                    if result[0]==NUMBER_OF_QUBITS:
                        break
                # print(pauli_to_circuit(elem))
                #Test
                # cProfile.run("find_p1s_p2s(elem)", filename="rand_circ_stats.txt")
        checks_properties=mymodule.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
        # Combine the circuit and checks.
        #todo.
        mymodule.write_outputs(circ_properties, checks_properties, file_number, file_info_path, file_qasm_path, output_file)

        print("execution time", time.time()-time0)
    print("done")