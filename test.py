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
    #Creat pauli group
    pauli_table=pauli_basis(5)
    # pauli_labels=list(map(lambda x:x.to_labels(),pauli_table))
    # pauli_labels=[val for sublist in pauli_labels for val in sublist] #Flatten the list
    # pauli_group=list(map(Pauli,pauli_labels))
    pauli_group=list(map(lambda x:x.to_matrix(),pauli_table)) #Convert to matrices
    pauli_group=[val for sublist in pauli_group for val in sublist] #Flatten the list
    for idx in range(2):
        #Random circuit using H, CNOT, S, T
        circ=my_random_circuit(5,3) #First number is num of qubits, second is depth
        # circ.draw('mpl')
        unitary = Operator(circ)
        print(circ)
        found=False
        for p1 in pauli_group:
            p2=unitary.dot(p1).dot(unitary.conjugate())
            if any(np.array_equal(p2.data,element) for element in pauli_group):
                print("p1:",p1)
                print("p2:",p2)
                found=True
        if not found:
            print("nothing found")
        # print(pauli_group[0])
        # print()
        # print(type(unitary.dot(pauli_group[0]).data))
        # plt.show()
        #Qiskit random circuit
        # circ2=random_circuit(5,10,2)
        # circ2.draw('mpl')
        # plt.show()