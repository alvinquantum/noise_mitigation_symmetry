#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from tkinter import N
# The number of threads is capped on Bebop.
# os.environ['OPENBLAS_NUM_THREADS'] = '1'

from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit import QuantumCircuit
import circgenerator
import checksfinder
import sys
import time
from sys import maxsize
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import QuantumRegister
from qiskit.converters.dag_to_circuit import dag_to_circuit
from qiskit.quantum_info import random_clifford, decompose_clifford
from qiskit.converters import circuit_to_dag
from multiprocessing import Pool
import psutil
import time
import os
import utilities
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit.quantum_info import Operator
from qiskit.visualization import circuit_drawer
from cirq.contrib.qasm_import import circuit_from_qasm
import json
# import pickle


if __name__ == "__main__":
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    MAX_LAYERS=int(sys.argv[3])
    NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))

    time0=time.time()
    #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    # func_params=FunctionParams()
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=f"qubits_{NUMBER_OF_QUBITS}"
    CHECKS_SUBDIR=os.path.join("min_weight_checks","checks")
    RAWS_SUBDIR="raws"
    CHECKS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, RAWS_SUBDIR)
    
    # #Create +1 phase pauli group
    # #Note that we can restrict the search for p2 to the +1 phase. The other results
    # #can all be recovered by multiplying by the corresponding phase, e.g., -i. 
    # pauli_table=pauli_basis(number_of_qubits)
    # # pauli_group_positive=deepcopy(pauli_table)
    # #Labels will be used to print pauli strings in the loop
    # pauli_labels=pauli_table.to_labels()

    # Generate random circuits
    # file_number=0
    # for _ in range(number_of_circuits):
    #     circ=generate_a_random_circuit(number_of_qubits, cnot_count)            
    #     #Save raw circuit
    #     write_circs_no_checks_to_qasm_file([circ], number_of_qubits, cnot_count, subdir=os.path.join(MAIN_SUBDIR, RAWS_SUBDIR))

    #     print(circ)
    #     circ_operations=circ.count_ops()
    #     print(circ_operations)

    # if number_of_qubits<=5:
    #     PARALLEL=False
    # else:
    #     PARALLEL=True
    PARALLEL=True #Can't use true because of the feed in parallel protocol.
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))
    #Load qasm files and find checks.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=circgenerator.filter_by_cnot_qubit(input_qasm_file_names, CNOT_COUNT, NUMBER_OF_QUBITS)
    for input_qasm_file_name in input_qasm_file_names:
        # print(input_qasm_file_name)

        # Skip the raw qasm if we already have an output.
        # if checksfinder.checks_file_exists(CHECKS_PATH, input_qasm_file_name):
        #     continue

        circ=QuantumCircuit.from_qasm_file(os.path.join(RAWS_PATH, input_qasm_file_name))
        circ_operations=circ.count_ops()
        circ_properties=circgenerator.CircuitProperties(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS, circ, circ_operations)
        checks_finder=checksfinder.ChecksFinder(NUMBER_OF_QUBITS, circ_properties.circ)
        # print(circ)
        # # Iterate over the qubits and try to get two p2 checks for each qubit
        layers_count=0
        check_layers=[]
        for qubit_idx in range(NUMBER_OF_QUBITS):
            time1=time.time()
            pauli_labels=checksfinder.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, NUMBER_OF_QUBITS)
            # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
            # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
            # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
            print(pauli_labels)
            # if PARALLEL:
            if False:
                count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_parallel_clifford(pool, pauli_labels)
                # with Pool(psutil.cpu_count(logical=False)) as pool:
                print(f"count: {count}")
            else:
                # count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels)
                count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels)


            if count>0:
                # count=0 #Reset the count. Since we only max of 2 checks.
                for idx, elem in enumerate(pauli_str_p1s):
                    # if idx<=1:
                    check_layers.append((elem, pauli_str_p2s[idx]))
                        # count+=1
                layers_count+=count
        if layers_count>0:
            # checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            # loop through the layers.
            increment_size=1
            start_size=0
            # for final_layers_idx in list(range(layers_count))[4:MAX_LAYERS:5]:
            #     temp_layers_count=final_layers_idx+1
            # circ_properties.circ=checksfinder.append_multilayer_checks_to_circ(circ, temp_layers_count, check_layers[0:temp_layers_count])
            print(layers_count)
            sets_of_check_layers=[check_layers[0:final_layers_idx+1] for final_layers_idx in list(range(layers_count))[start_size:MAX_LAYERS:increment_size]]
            circs=[circ]*len(sets_of_check_layers)
            temp_layers_counts=[len(elem) for elem in sets_of_check_layers]
            result_circs =pool.starmap(checksfinder.append_multilayer_checks_to_circ_clifford, zip(circs, temp_layers_counts, sets_of_check_layers))
            # circ_properties.circ=checksfinder.append_multilayer_checks_to_circ_clifford(circ, temp_layers_count, check_layers[0:temp_layers_count])

            # print(circ_properties.circ)
            for idx, circ in enumerate(result_circs):
                circ_properties.circ=circ
                temp_layers_count=temp_layers_counts[idx]
                #Draw to file
                output_qasm_file_name=checksfinder.generate_checks_file_name_from_raw_qasm(input_qasm_file_name, temp_layers_count)
                print(output_qasm_file_name)
                info_file_name=checksfinder.generate_obj_file_name_from_checks_file_name(output_qasm_file_name)
                txt_file_name=checksfinder.generate_txt_file_name_from_checks_file_name(output_qasm_file_name)
                circuit_drawer(circ_properties.circ, filename=os.path.join(CHECKS_PATH, txt_file_name))        
                output_txt_file=open(os.path.join(CHECKS_PATH, txt_file_name), "a")
                circuit_number=checksfinder.get_circuit_number(output_qasm_file_name)
                checksfinder.write_outputs_multilayer(circ_properties, os.path.join(CHECKS_PATH, output_qasm_file_name))
                
                output_txt_file.close()
        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    if PARALLEL:
        pool.close() #Stops passing jobs to processes.
        pool.join() #Waits for processes to finish.
    print("done")