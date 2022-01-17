#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
# The number of threads is capped on Bebop.
# os.environ['OPENBLAS_NUM_THREADS'] = '1'

import sys
import psutil
import json
from qiskit.quantum_info import Operator
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
import circgenerator
import utilities
from qiskit.visualization import circuit_drawer
from multiprocessing import Pool
import time
import checksfinder
import cProfile

def initialize(ABS_TOL_ARG, circ_properties_arg,
                    pauli_labels_arg, table_length_arg):
    '''Initialize globals for pool'''
    global ABS_TOL, circ_properties, pauli_labels, table_length
    #Passed from the main
    #Program parameters passed from main
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    circ_properties=circ_properties_arg
    ABS_TOL=ABS_TOL_ARG

    #Tables passed from main
    pauli_labels=pauli_labels_arg
    table_length=table_length_arg

def get_checks_parallel(pauli_group_elem):
    '''Wrapper function. Get the checks in parallel.'''
    global ABS_TOL, circ_properties, pauli_labels, table_length
    return circgenerator.find_checks_sym(pauli_group_elem, circ_properties.circ)

if __name__ == "__main__":
    time0=time.time()
    # #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    # func_params=FunctionParams()
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    NUMBER_OF_CIRCUITS=int(sys.argv[3])
    if NUMBER_OF_QUBITS<=5:
        PARALLEL=False
    else:
        PARALLEL=True
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    ABS_TOL=.2*2**(NUMBER_OF_QUBITS-1)
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    # code_dir1=os.path.dirname(os.path.realpath("__file__"))
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    # print(code_dir1)
    # print(code_dir2)
    SUBDIR="data"
    # BASE_FILE_PATH=code_dir+SUBDIR+"depth"+ str(DEPTH) +"_"
    # BASE_FILE_PATH=CODE_DIR+SUBDIR+"qubits_"+str(NUMBER_OF_QUBITS)+"_CNOTS_"+ str(CNOT_COUNT) +"_"
    BASE_PATH=os.path.join(CODE_DIR, SUBDIR)
    
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

    table_length=len(pauli_table)
    #main loop
    file_number=0
    for _ in range(NUMBER_OF_CIRCUITS):
        time1=time.time()
        # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
        output_file_name=f"qubits_{NUMBER_OF_QUBITS}_CNOTS_{CNOT_COUNT}_circuit_{file_number}_.txt"
        while os.path.isfile(os.path.join(BASE_PATH, output_file_name)):
            file_number+=1
            # output_file_name=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
            output_file_name=f"qubits_{NUMBER_OF_QUBITS}_CNOTS_{CNOT_COUNT}_circuit_{file_number}_.txt"
        info_file_name=f"qubits_{NUMBER_OF_QUBITS}_CNOTS_{CNOT_COUNT}_circuit_{file_number}_.obj"
        qasm_file_name=f"qubits_{NUMBER_OF_QUBITS}_CNOTS_{CNOT_COUNT}_circuit_{file_number}_.qasm"
        #Random circuit.
        circ=circgenerator.random_circuit_cnot(NUMBER_OF_QUBITS, CNOT_COUNT)

        # We use this to find P1 and P2.
        unitary = Operator(circ)
        #Draw to file
        circuit_drawer(circ, filename=os.path.join(BASE_PATH, output_file_name))
        print(circ)
        
        output_file=open(os.path.join(BASE_PATH, output_file_name), "a")
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

        circ_properties=circgenerator.CircuitProperties(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS, circ, circ_operations)
        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        if PARALLEL:
            with Pool(psutil.cpu_count(logical=False), initialize, initargs=(ABS_TOL, circ_properties,
                    pauli_labels, table_length)) as pool:
                #In some cases pool.imap_unordered needs to be wrapped in list in order to return properly. 
                #see: https://stackoverflow.com/questions/5481104/multiprocessing-pool-imap-broken
                for result in pool.imap_unordered(get_checks_parallel, reversed(pauli_labels), chunksize=100):
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
                result=checksfinder.find_checks_sym(elem, circ)
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
        checks_properties=utilities.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
        # Combine the circuit and checks.
        if p2_weights[0]>0:
            circ_properties.circ=checksfinder.append_checks_to_circ(circ_properties, checks_properties)
        checksfinder.write_outputs(circ_properties, checks_properties, file_number, os.path.join(BASE_PATH, info_file_name),os.path.join(BASE_PATH, qasm_file_name), output_file)

        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    print("done")