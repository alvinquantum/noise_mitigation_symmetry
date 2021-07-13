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


if __name__ == "__main__":
    #Program parameters
    number_of_qubits=5
    depth=10
    number_of_circuits=1

    #Create +1 phase pauli group
    pauli_table=pauli_basis(number_of_qubits)
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
    for idx in range(number_of_circuits):
        #Modified Random circuit using H, CNOT, S, T.
        circ=my_random_circuit(number_of_qubits,depth)
        # Don't delete! Qiskit random circuit
        # circ=random_circuit(5,10,2)
        #Draw
        circ.draw('mpl')
        print(circ)
        unitary = Operator(circ)
        found=False
        # Number of found p2 matrices.
        count=0

        for idx1, p1 in enumerate(pauli_group):
            #U.p1=p2.U ---->U.p1.U^\dagger=p2. Operator class so we need .data to access numpy array.
            p2=unitary.dot(p1).dot(unitary.adjoint()).data
            #Only need p2 with +1 phase since the global phase can be absorbed into p1. Faster this way.
            for idx2, element in enumerate(pauli_group_positive):
                # allclose checks if the values are within tolerance of 10^-8.
                if np.allclose(p2, element):
                    #Have to check which part of the table p1 belongs to so we can output the correct phase.
                    if idx1-table_length<=0:
                        print("p1: +1", pauli_table[idx1 % table_length])
                    elif idx1-2*table_length<=0:
                        print("p1: -1",pauli_table[idx1 % table_length])
                    elif idx1-3*table_length<=0:
                        print("p1: +j", pauli_table[idx1 % table_length])
                    else:
                        print("p1: -j", pauli_table[idx1 % table_length])

                    print("p2: +1",pauli_table[idx2])
                    print()
                    found=True
                    count+=1
        if not found:
            print("nothing found")
        print("Found: ", count)
        plt.show()