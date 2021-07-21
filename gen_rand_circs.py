#!/usr/bin/env python3

import numpy as np
from qiskit import QuantumCircuit
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector
from qiskit.quantum_info import Operator
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit.quantum_info import Pauli
from mymodule import my_random_circuit
from qiskit.visualization import circuit_drawer
import os
import pickle
# from qiskit.circuit import ControlledGate
# from qiskit.extensions import UnitaryGate

#Program parameters
NUMBER_OF_QUBITS=2
DEPTH=6
NUMBER_OF_CIRCUITS=1
#Paths for outputs and pickle file of circuit
code_dir=os.path.dirname(os.path.realpath('__file__'))
subdir="/data/"
base_file_path=code_dir+subdir+"depth"+ str(DEPTH) +"_"

def create_controlU(npmat, number_of_qubits):
    '''Returns a controlled operation. Type is np.array'''
    return np.kron(np.eye(2**number_of_qubits),np.array([[1/2,1/2],[1/2,1/2]]))+np.kron(npmat, np.array([[1/2,-1/2],[-1/2,1/2]]))

def get_weight(pauli_string):
    '''Gets the weight of a Pauli string.'''
    count=0
    for character in pauli_string:
        if character!="I":
            count+=1
    return count    

def check_p2(control_p1, control_p2, unitary, number_of_qubits):
    '''Sanity check for p2. U\otimes I- ControlP2^\dagger(U\otimes I)ControlP1==0'''
    assert np.allclose(np.kron(unitary.data,np.eye(2))-control_p2.T.conjugate().dot(np.kron(unitary.data,np.eye(2))).dot(control_p1), np.zeros(2**(number_of_qubits+1))), "wrong p2"

if __name__ == "__main__":
    print("running...")
    #Create +1 phase pauli group
    pauli_table=pauli_basis(NUMBER_OF_QUBITS)
    #Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()
    #Convert to matrices
    pauli_group_positive=list(map(lambda x:x.to_matrix(),pauli_table))
    #Flatten the list
    pauli_group_positive=[val for sublist in pauli_group_positive for val in sublist]
    # Make the other phases and merge
    pauli_group_negative=[element * -1 for element in pauli_group_positive]
    pauli_group_positiveI=[element * 1j for element in pauli_group_positive]
    pauli_group_negativeI=[element * -1j for element in pauli_group_positive]
    pauli_group=pauli_group_positive+pauli_group_negative+pauli_group_positiveI+pauli_group_negativeI
    table_length=len(pauli_table)
    #main loop
    file_number=0
    for idx in range(NUMBER_OF_CIRCUITS):
        output_file_path=base_file_path+ str(file_number) +".txt"
        while os.path.isfile(output_file_path):
            file_number+=1
            output_file_path=base_file_path+ str(file_number) +".txt"
        circ_file_path=base_file_path+ str(file_number)+".obj"

        #Modified Random circuit using H, CNOT, S, T.
        circ=my_random_circuit(NUMBER_OF_QUBITS,DEPTH)
        # Don't delete! Qiskit random circuit
        # circ=random_circuit(5,10)

        unitary = Operator(circ)
        max_pauli_weight=0
        max_pauli_str_p1=""
        max_pauli_str_p2=""
        # Number of found p2 matrices.
        count=0

        #Draw
        circ.draw('mpl')
        circuit_drawer(circ, filename=output_file_path)
        print(circ)
        
        with open(output_file_path, "a") as file:
            print(file=file)
            # Count the number of CNOT gates
            circ_operations=circ.count_ops()
            print(circ_operations, file=file)
            print(circ_operations)
            #We need to iterate through the entire pauli group for p1.
            for idx1, p1 in enumerate(pauli_group):
                #U.p1=p2.U ---->U.p1.U^\dagger=p2. Operator class so we need .data to access numpy array.
                p2=unitary.dot(p1).dot(unitary.adjoint()).data
                #Only need p2 with +1 phase since the global phase can be absorbed into p1. Faster this way.
                for idx2, element in enumerate(pauli_group_positive):
                    # allclose checks if the values are within tolerance of 10^-8.
                    if np.allclose(p2, element):
                        #Have to check which part of the table p1 belongs to so we can print the correct phase.
                        if idx1-table_length<=0:
                            p1_str="+1"+pauli_labels[idx1 % table_length]
                        elif idx1-2*table_length<=0:
                            p1_str="-1"+pauli_labels[idx1 % table_length]
                        elif idx1-3*table_length<=0:
                            p1_str="+j"+pauli_labels[idx1 % table_length]
                        else:
                            p1_str="-j"+pauli_labels[idx1 % table_length]
                        
                        #The phase for p2 can be passed to p1.
                        p2_str="+1"+pauli_labels[idx2]
                        print("p1: ", p1_str, file=file)
                        print("p1: ", p1_str)
                        # print(p1, file=file)
                        print("p2: ", p2_str, file=file)
                        print("p2: ", p2_str)
                        # print(p2, file=file)

                        # Print the weight. We care about P2 weight since we commute p1 through U.
                        p2_weight=get_weight(pauli_labels[idx2])
                        print("Pauli weight P2: ", p2_weight, file=file)
                        print("Pauli weight P2: ", p2_weight)
                        if p2_weight>max_pauli_weight:
                            max_pauli_weight=p2_weight
                            max_pauli_str_p1=p1_str
                            max_pauli_str_p2=p2_str
                        print(file=file)
                        print()
                        #Sanity check. Can comment out.
                        # control_p1=create_controlU(p1, NUMBER_OF_QUBITS)
                        # control_p2=create_controlU(p2, NUMBER_OF_QUBITS)
                        # check_p2(control_p1, control_p2, unitary, NUMBER_OF_QUBITS)

                        count+=1
                        break

            #Outputs
            if count==1:
                print("nothing found: trivial solution", file=file)
                print("nothing found: trivial solution")
            else:
                print("Found Matches: ", count, file=file)
                print("Found Matches: ", count)
                cnot_count=0
                if "cx" in circ_operations:
                    cnot_count=circ_operations["cx"]
                print("CNOT count: ", cnot_count, file=file)
                print("CNOT count: ", cnot_count)
                print("Max Weight: ", max_pauli_weight, file=file)
                print("Max Weight: ", max_pauli_weight)
                print("P1 that creates max P2: ", max_pauli_str_p1, file=file)
                print("P1 that creates max P2: ", max_pauli_str_p1)
                print("Max P2: ", max_pauli_str_p2, file=file)
                print("Max P2: ", max_pauli_str_p2)
        # Dump all the info into a pickle
        circ_file=open(circ_file_path, "wb")
        pickle.dump({"circ": circ, "max_pauli_weight": max_pauli_weight, "max_pauli_str_p1": max_pauli_str_p1, "max_pauli_str_p2": max_pauli_str_p2}, circ_file)
    print("done")