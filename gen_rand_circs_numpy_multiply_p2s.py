#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
# The number of threads is capped on Bebop. This also gives a performance boost in general.
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from re import M
import numpy as np
import math, sys, pickle, ctypes, psutil, json
from qiskit import QuantumCircuit
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

def create_controlU(npmat, number_of_qubits):
    '''Returns a controlled operation. Returns: np.array'''
    return np.kron(np.array([[1,0],[0,0]]),np.eye(2**number_of_qubits))+np.kron(np.array([[0,0],[0,1]]),npmat)

def get_weight(pauli_string):
    '''Gets the weight of a Pauli string. Returns: int'''
    count=0
    for character in pauli_string:
        if character!="I":
            count+=1
    return count    

def check_p2(control_p1, control_p2, unitary, number_of_qubits):
    '''Sanity check for p2. I\otimes U- ControlP2(I\otimes U)ControlP1==0 Returns: bool'''
    return np.allclose(np.kron(np.eye(2),unitary.data)-control_p2.dot(np.kron(np.eye(2),unitary.data)).dot(control_p1), np.zeros(2**(number_of_qubits+1)))

def find_p1s_p2s(pauli_group_tuple):
    # print("here")
    # print(pauli_group_tuple)
    '''For multiprocessing.'''
    global unitary, NUMBER_OF_QUBITS, CNOT_COUNT, ABS_TOL 
    global pauli_labels, pauli_group_positive, table_length
    global count, p2_weights, pauli_str_p1s, pauli_str_p2s
    global unitary
    idx1=pauli_group_tuple[0]
    # to_matrix() returns a list.
    p1=pauli_group_tuple[1]#.to_matrix()[0]
    # if idx1 >= table_length:
    #     p1 = p1*-1
    #U.p1=p2.U ---->U.p1.U^\dagger=p2. Operator class so we need .data to access numpy array.
    p2=unitary.dot(p1).dot(unitary.adjoint()).data
    #Sanity check. Can comment out.
    control_p1=create_controlU(p1, NUMBER_OF_QUBITS)
    control_p2=create_controlU(p2, NUMBER_OF_QUBITS)
    if not check_p2(control_p1, control_p2, unitary, NUMBER_OF_QUBITS):
        return
    # Check if p2 is traceless. All elements of the pauli group are traceless except identity.
    if not math.isclose(0.0,np.trace(p2), abs_tol=ABS_TOL):
        return
    #Only need p2 with +1 phase since the global phase can be absorbed into p1. Faster this way.
    for idx2, element in enumerate(pauli_group_positive):
        # to_matrix() returns a list.
        # element=label_element.to_matrix()[0]
        # allclose checks if the values are within tolerance of atol=10^-8.
        if np.allclose(p2, element):
            #Have to check which part of the table p1 belongs to so we can print the correct phase.
            if idx1-table_length<0:
                p1_str="+1"+pauli_labels[idx1 % table_length]
            else:
                p1_str="-1"+pauli_labels[idx1 % table_length]
            # elif idx1-2*table_length<0:
            #     p1_str="-1"+pauli_labels[idx1 % table_length]
            # elif idx1-3*table_length<0:
            #     p1_str="+j"+pauli_labels[idx1 % table_length]
            # else:
            #     p1_str="-j"+pauli_labels[idx1 % table_length]
            #P2 is always +1 phase.
            p2_str="+1"+pauli_labels[idx2]
            print("p1: ", p1_str)
            print("p2: ", p2_str)

            # Print the weight. We care about P2 weight since we commute p1 through U.
            p2_weight=get_weight(pauli_labels[idx2])
            print("Pauli weight P2: ", p2_weight)
            print()
            #Need to lock the value so it doesn't change while checking. Since we write only
            #after locking the weight we don't need to lock the other values.
            with count.get_lock():
                count.value+=1
                if p2_weight>p2_weights[0]:
                    # We store the max values in the beginning.
                    p2_weights[0]=p2_weight
                    pauli_str_p1s[0]=p1_str
                    pauli_str_p2s[0]=p2_str
                pauli_str_p1s.append(p1_str)
                pauli_str_p2s.append(p2_str)
                p2_weights.append(p2_weight)

            # #Need to lock the value so it doesn't change. Nonatomic operation.
            # with count.get_lock():
            #     count.value+=1
            #We found p2 so just return.
            return

def initialize(unitary_arg, NUMBER_OF_QUBITS_ARG, CNOT_COUNT_ARG, ABS_TOL_ARG,
    pauli_labels_arg, pauli_group_positive_arg, table_length_arg, count_arg,
    p2_weights_arg, pauli_str_p1s_arg, pauli_str_p2s_arg):
    '''Initialize globals for pool'''
    global unitary, NUMBER_OF_QUBITS, CNOT_COUNT, ABS_TOL 
    global pauli_labels, pauli_group_positive, table_length
    global count, p2_weights, pauli_str_p1s, pauli_str_p2s
    global unitary
    #Passed from the main
    unitary=unitary_arg
    #Program parameters passed from main
    NUMBER_OF_QUBITS=NUMBER_OF_QUBITS_ARG
    CNOT_COUNT=CNOT_COUNT_ARG

    #Tables passed from main
    pauli_labels=pauli_labels_arg
    pauli_group_positive=pauli_group_positive_arg
    table_length=table_length_arg

    #Variables that can change
    p2_weights=p2_weights_arg
    pauli_str_p1s=pauli_str_p1s_arg
    pauli_str_p2s=pauli_str_p2s_arg
    # Number of found p2 matrices.
    count=count_arg
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    ABS_TOL=ABS_TOL_ARG


if __name__ == "__main__":
    # #Disable qiskit parallel.
    # os.environ['QISKIT_IN_PARALLEL'] = 'FALSE'
    print("running...")
    #Program parameters
    PARALLEL=True
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    NUMBER_OF_CIRCUITS=int(sys.argv[3])
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
    # pauli_group=pauli_table._add(pauli_table)
    
    #Convert to matrices
    pauli_group_positive=list(map(lambda x:x.to_matrix(),pauli_table))
    #Flatten the list
    pauli_group_positive=[val for sublist in pauli_group_positive for val in sublist]
    # Make the other phases and merge
    pauli_group_negative=[element * -1 for element in pauli_group_positive]
    pauli_group=pauli_group_positive+pauli_group_negative

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

        #Parallel stuff
        #We need to use Value type for locking and Manager for persistence.
        manager0=Manager()
        manager1=Manager()
        manager2=Manager()
        # max_pauli_weight=Value("i",0)
        p2_weights=manager0.list([0])
        count=Value("i",0)
        pauli_str_p1s=manager1.list([""])
        pauli_str_p2s=manager2.list([""])

        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        if PARALLEL: #we were using int(psutil.cpu_count(logical=False)*85/100) for the number of processes before (works on windows). 
                     #os.sched_getaffinity(0) works on linux
            with Pool(psutil.cpu_count(logical=False), initialize, initargs=(unitary, NUMBER_OF_QUBITS, CNOT_COUNT, ABS_TOL,
                    pauli_labels, pauli_group_positive, table_length, count,
                    p2_weights, pauli_str_p1s, pauli_str_p2s)) as pool:
                #needs to be wrapped in list in order to execute. see: https://stackoverflow.com/questions/5481104/multiprocessing-pool-imap-broken
                list(pool.imap_unordered(find_p1s_p2s, enumerate(pauli_group)))
                # pool.map(find_p1s_p2s, enumerate(pauli_group))
        else:
            for elem in enumerate(pauli_group):
                #Test
                find_p1s_p2s(elem)
                # cProfile.run("find_p1s_p2s(elem)", filename="rand_circ_stats.txt")
        assert len(pauli_str_p1s)==len(pauli_str_p2s), "number of p1's and p2's don't match."
        assert len(pauli_str_p2s)==len(p2_weights), "number of p2's and weights don't match."
        # Count doesn't count the extra storage of max vales in the beginning of the lists of solutions (e.g. pauli_str_p1s)
        assert len(pauli_str_p2s)==count.value+1, "number of p2's and counts of solutions don't match."

        #Outputs
        if count.value==0:
            output_file.write("\n")
            output_file.write("nothing found: trivial solution\n")
            print("nothing found: trivial solution")
        else:
            for index1, strp1 in enumerate(pauli_str_p1s):
                if index1!=0: #The max stuff are stored in the beginning, which we print at the end.
                    output_file.write("\n")
                    output_file.write("p1: "+ str(strp1)+ "\n")
                    output_file.write("p2: "+ str(pauli_str_p2s[index1])+ "\n")
                    output_file.write("Pauli weight P2: "+ str(p2_weights[index1])+ "\n")
            output_file.write("\n")
            output_file.write("Found Matches: "+ str(count.value)+ "\n")
            output_file.write("Max Weight: "+ str(p2_weights[0]) + "\n")
            output_file.write("P1 that creates max P2: "+ str(pauli_str_p1s[0])+ "\n")
            output_file.write("Max P2: " + str(pauli_str_p2s[0])+ "\n")
            print("Found Matches: ", count.value)   
            print("Max P2 Weight: ", p2_weights[0])
            print("P1 that creates max P2: ", pauli_str_p1s[0])
            print("Max P2: ", pauli_str_p2s[0])

        cnot_count=0
        if "cx" in circ_operations:
            cnot_count=circ_operations["cx"]

        rz_count=0
        if "rz" in circ_operations:
            rz_count=circ_operations["rz"]
        output_file.write("Circuit no: " + str(file_number)+ "\n")
        output_file.write("Qubits: " + str(NUMBER_OF_QUBITS)+ "\n")
        output_file.write("CNOT count: "+ str(cnot_count) + "\n")
        output_file.write("RZ count: "+ str(rz_count) + "\n")

        print("Qubits: ", str(NUMBER_OF_QUBITS))
        print("Circuit no: ", str(file_number))
        print("CNOT count: ", cnot_count)
        print("RZ count: ", rz_count)

        # # Dump all the info into a pickle
        circ_file=open(file_info_path, "wb")
        pickle.dump({"cx": cnot_count, "rz": rz_count, "qubits": NUMBER_OF_QUBITS, "circuit_num" : file_number, "found_matches": count.value,
             "max_pauli_weight": p2_weights[0], "max_pauli_str_p1": pauli_str_p1s[0], "max_pauli_str_p2": pauli_str_p2s[0]}, circ_file)
        circ.qasm(filename=file_qasm_path)
        # # Close the files. 
        circ_file.close()
        output_file.close()
        print("execution time", time.time()-time0)
    print("done")