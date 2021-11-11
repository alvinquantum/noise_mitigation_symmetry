#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from re import M
import numpy as np
import math, sys, os, pickle, ctypes, psutil, json
from qiskit import QuantumCircuit
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector, Pauli, Operator
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from mymodule import my_random_circuit
from qiskit.visualization import circuit_drawer
from multiprocessing import Pool, Value, Array, Manager
from copy import deepcopy
from pprint import pprint
# from qiskit.circuit import ControlledGate
# from qiskit.extensions import UnitaryGate

def create_controlU(npmat, number_of_qubits):
    '''Returns a controlled operation. Type is np.array'''
    return np.kron(np.array([[1,0],[0,0]]),np.eye(2**number_of_qubits))+np.kron(np.array([[0,0],[0,1]]),npmat)

def get_weight(pauli_string):
    '''Gets the weight of a Pauli string.'''
    count=0
    for character in pauli_string:
        if character!="I":
            count+=1
    return count    

def check_p2(control_p1, control_p2, unitary, number_of_qubits):
    '''Sanity check for p2. I\otimes U- ControlP2^\dagger(I\otimes U)ControlP1==0'''
    return np.allclose(np.kron(np.eye(2),unitary.data)-control_p2.dot(np.kron(np.eye(2),unitary.data)).dot(control_p1), np.zeros(2**(number_of_qubits+1)))

def find_p1s_p2s(pauli_group_tuple):
    '''For multiprocessing.'''
    global NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS, ABS_TOL, code_dir, subdir, base_file_path
    global pauli_labels, pauli_group_positive, pauli_group, table_length
    global count, max_pauli_str_p1, max_pauli_str_p2, max_pauli_weight
    global unitary
    idx1=pauli_group_tuple[0]
    p1=pauli_group_tuple[1].to_matrix()[0]
    # print(type(p1))
    # print(p1)
    # for idx1, p1 in enumerate(pauli_group):
    #U.p1=p2.U ---->U.p1.U^\dagger=p2. Operator class so we need .data to access numpy array.
    p2=unitary.dot(p1).dot(unitary.adjoint()).data
    # print(idx1)
    #Sanity check. Can comment out.
    control_p1=create_controlU(p1, NUMBER_OF_QUBITS)
    control_p2=create_controlU(p2, NUMBER_OF_QUBITS)
    if not check_p2(control_p1, control_p2, unitary, NUMBER_OF_QUBITS):
        return
    #Check if p2 is traceless. All elements of the pauli group are traceless except identity.
    if not math.isclose(np.trace(p2),0.0, abs_tol=ABS_TOL):
        return
    #Only need p2 with +1 phase since the global phase can be absorbed into p1. Faster this way.
    for idx2, label_element in enumerate(pauli_group_positive):
        # allclose checks if the values are within tolerance of 10^-8.
        element=label_element.to_matrix()[0]
        if idx2 >= table_length:
            element = element*-1
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
             
            #The phase for p2 can be passed to p1.
            p2_str="+1"+pauli_labels[idx2]
            # print("p1: ", p1_str, file=output_file)
            print("p1: ", p1_str)
            # print(p1, file=file)
            # print("p2: ", p2_str, file=output_file)
            print("p2: ", p2_str)
            # print(p2, file=file)

            # Print the weight. We care about P2 weight since we commute p1 through U.
            p2_weight=get_weight(pauli_labels[idx2])
            # print("Pauli weight P2: ", p2_weight, file=output_file)
            print("Pauli weight P2: ", p2_weight)
            #Need to lock the value so it doesn't change while checking. Since we write only
            #after locking the weight we don't need to lock the other values.
            with count.get_lock():
                count.value+=1
                if p2_weight>max_pauli_weight[0]:
                    max_pauli_weight[0]=p2_weight
                    max_pauli_str_p1[0]=p1_str
                    max_pauli_str_p2[0]=p2_str
                max_pauli_str_p1.append(p1_str)
                max_pauli_str_p2.append(p2_str)
                max_pauli_weight.append(p2_weight)
            # print(file=output_file)
            print()

            # #Need to lock the value so it doesn't change. Nonatomic operation.
            # with count.get_lock():
            #     count.value+=1
            #We found p2 so just return.
            return

def initialize(unitary_arg, NUMBER_OF_QUBITS_ARG, CNOT_COUNT_ARG, NUMBER_OF_CIRCUITS_ARG, ABS_TOL_ARG,
    pauli_labels_arg, pauli_group_positive_arg, pauli_group_arg, table_length_arg, count_arg,
    max_pauli_weight_arg, max_pauli_str_p1_arg, max_pauli_str_p2_arg):
    '''Initialize globals for pool'''
    global NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS, ABS_TOL, code_dir, subdir, base_file_path
    global pauli_labels, pauli_group_positive, pauli_group, table_length
    global count, max_pauli_str_p1, max_pauli_str_p2, max_pauli_weight
    global unitary
    #Passed from the main
    unitary=unitary_arg
    #Program parameters passed from main
    NUMBER_OF_QUBITS=NUMBER_OF_QUBITS_ARG
    CNOT_COUNT=CNOT_COUNT_ARG
    NUMBER_OF_CIRCUITS=NUMBER_OF_CIRCUITS_ARG
    # POOL=Pool(psutil.cpu_count(logical=False))

    #Tables passed from main
    pauli_labels=pauli_labels_arg
    pauli_group_positive=pauli_group_positive_arg
    pauli_group=pauli_group_arg
    table_length=table_length_arg

    #Variables that can change
    max_pauli_weight=max_pauli_weight_arg
    max_pauli_str_p1=max_pauli_str_p1_arg
    max_pauli_str_p2=max_pauli_str_p2_arg
    # Number of found p2 matrices.
    count=count_arg
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    ABS_TOL=ABS_TOL_ARG


if __name__ == "__main__":
    print("running...")
    #Program parameters
    NUMBER_OF_QUBITS=5#int(sys.argv[1])
    CNOT_COUNT=50#int(sys.argv[2])
    NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    # Absolute tolerance for checking if the trace of p2 is close to zero with the isclose function.
    ABS_TOL=.2*2**(NUMBER_OF_QUBITS-1)
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    # code_dir1=os.path.dirname(os.path.realpath("__file__"))
    code_dir=sys.path[0]
    # print(code_dir1)
    # print(code_dir2)
    subdir="/data/"
    # base_file_path=code_dir+subdir+"depth"+ str(DEPTH) +"_"
    base_file_path=code_dir+subdir+"qubits_"+str(NUMBER_OF_QUBITS)+"_CNOTS_"+ str(CNOT_COUNT) +"_"

    #Create +1 phase pauli group
    pauli_table=pauli_basis(NUMBER_OF_QUBITS)
    pauli_group_positive=list(deepcopy(pauli_table))
    #Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()
    pauli_group=list(pauli_table._add(pauli_table))
    # print(type(pauli_group))
    # print(pauli_group_positive)
    # print(pauli_group)
    # print(pauli_table)
    # print(len(pauli_table))
    # print(len(pauli_table._add(pauli_table)))
    # print(pauli_table._add(pauli_table))
    #Convert to matrices
    # pauli_group_positive=list(map(lambda x:x.to_matrix(),pauli_table))
    #Flatten the list
    # pauli_group_positive=[val for sublist in pauli_group_positive for val in sublist]
    # Make the other phases and merge
    # pauli_group_negative=[element * -1 for element in pauli_group_positive]
    # pauli_group_positiveI=[element * 1j for element in pauli_group_positive]
    # pauli_group_negativeI=[element * -1j for element in pauli_group_positive]
    # Note that we can restrict the search for p1 to the +/-1 phases. This is due to passing the phase
    # from p2 to p1 and then realizing that the eigenvalues of p2 is restricted to +/-1. The unitary
    # conjugating p1 only rotates the eigenvectors. Therefore, the eigenvalues of p1 must also be +/-1
    # which eliminates the need to check for +/-j phases. This is explained in the paper.
    # pauli_group=pauli_group_positive+pauli_group_negative
    # pauli_group=pauli_group_positive+pauli_group_negative+pauli_group_positiveI+pauli_group_negativeI


    table_length=len(pauli_table)
    #main loop
    file_number=0
    for idx in range(NUMBER_OF_CIRCUITS):
        output_file_path=base_file_path+ "circuit_" + str(file_number) +"_.txt"
        while os.path.isfile(output_file_path):
            file_number+=1
            output_file_path=base_file_path+ "circuit_" + str(file_number) +"_.txt"
        circ_file_path=base_file_path+ "circuit_" + str(file_number)+"_.obj"
        #Modified Random circuit using H, CNOT, S, T.
        circ=my_random_circuit(NUMBER_OF_QUBITS, CNOT_COUNT)
        # # Load circuit.
        # circ_file=open(base_file_path+"3.obj", "rb")
        # circ=pickle.load(circ_file)["circ"]
        # circ_file.close()
        # Don't delete! Qiskit random circuit
        # circ=random_circuit(5,10)

        # We use this to find P1 and P2.
        unitary = Operator(circ)
        #Draw
        # circ.draw(filename=output_file_path)
        circuit_drawer(circ, filename=output_file_path)
        # print("here1")
        # pprint(vars(circ.data))
        print(circ)
        
        output_file=open(output_file_path, "a")
        # with open(output_file_path, "a") as file:
        # print(file=output_file)
        output_file.write("\n")
        # Count the number of CNOT gates
        circ_operations=circ.count_ops()
        # print(circ_operations, file=output_file)
        output_file.write(json.dumps(circ_operations))
        print(circ_operations)
        #We need to use these types for locking and persistence.
        managerp0=Manager()
        # max_pauli_weight=Value("i",0)
        max_pauli_weight=managerp0.list([0])
        count=Value("i",0)
        # max_pauli_str_p2=Array(ctypes.c_char, "+1"+"I"*NUMBER_OF_QUBITS)

        # max_pauli_str_p1[0]='+1jxxxxx'.encode()
        # print(max_pauli_str_p1[0])
        # print(max_pauli_str_p1[0].decode())

        # #Parallel stuff
        managerp1=Manager()
        managerp2=Manager()
        max_pauli_str_p1=managerp1.list([""])
        max_pauli_str_p2=managerp2.list([""])

        # unitary_arg, NUMBER_OF_QUBITS_ARG, DEPTH_ARG, NUMBER_OF_CIRCUITS_ARG, ABS_TOL_ARG,
        # pauli_labels_arg, pauli_group_positive_arg, pauli_group_arg, table_length_arg, count_arg,
        # max_pauli_weight_arg, max_pauli_str_p1_arg, max_pauli_str_p2_arg
        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        with Pool(psutil.cpu_count(logical=False), initialize, initargs=(unitary, NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS, ABS_TOL,
            pauli_labels, pauli_group_positive, pauli_group, table_length, count, max_pauli_weight, max_pauli_str_p1, 
            max_pauli_str_p2)) as pool:
            pool.map(find_p1s_p2s, enumerate(pauli_group))
        # print("Max P2: ", max_pauli_str_p2[0].decode())
        #Outputs
        if count.value==0:
            # print("nothing found: trivial solution", file=output_file)
            output_file.write("\n")
            output_file.write("nothing found: trivial solution\n")
            print("nothing found: trivial solution")
        else:
            for index1, strp1 in enumerate(max_pauli_str_p1):
                if index1!=0:
                    output_file.write("\n")
                    output_file.write("p1: "+ str(strp1)+ "\n")
                    output_file.write("p2: "+ str(max_pauli_str_p2[index1])+ "\n")
                    output_file.write("Pauli weight P2: "+ str(max_pauli_weight[index1])+ "\n")
            # print("Found Matches: ", count.value, file=output_file)
            output_file.write("\n")
            output_file.write("Found Matches: "+ str(count.value)+ "\n")
            print("Found Matches: ", count.value)
            # print("Max Weight: ", max_pauli_weight.value, file=output_file)
            output_file.write("Max Weight: "+ str(max_pauli_weight[0]) + "\n")
            print("Max Weight: ", max_pauli_weight[0])
            # print("P1 that creates max P2: ", max_pauli_str_p1[0], file=output_file)
            output_file.write("P1 that creates max P2: "+ str(max_pauli_str_p1[0])+ "\n")
            print("P1 that creates max P2: ", max_pauli_str_p1[0])
            # print("Max P2: ", max_pauli_str_p2[0], file=output_file)
            output_file.write("Max P2: " + str(max_pauli_str_p2[0])+ "\n")
            print("Max P2: ", max_pauli_str_p2[0])
        cnot_count=0
        if "cx" in circ_operations:
            cnot_count=circ_operations["cx"]
        # print("CNOT count: ", cnot_count, file=output_file)
        output_file.write("CNOT count: "+ str(cnot_count) + "\n")
        print("CNOT count: ", cnot_count)

        # Dump all the info into a pickle
        circ_file=open(circ_file_path, "wb")
        pickle.dump({"circ": circ, "max_pauli_weight": max_pauli_weight[0], "max_pauli_str_p1": max_pauli_str_p1[0], "max_pauli_str_p2": max_pauli_str_p2[0]}, circ_file)
        # Close the files. 
        circ_file.close()
        output_file.close()

    print("done")