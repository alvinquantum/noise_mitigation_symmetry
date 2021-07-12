#!/usr/bin/env python3

import numpy as np
from qiskit import QuantumCircuit
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector
from qiskit.quantum_info import Operator
from qiskit.quantum_info import Pauli
from qiskit.circuit.random import random_circuit
from mymodule import my_random_circuit

# Create random cirquit

if __name__ == "__main__":
    #Random circuit using H, CNOT, S, T
    circ=my_random_circuit(5,30,2)
    print(circ)
    circ.draw('mpl')
    plt.show()
    #Qiskit random circuit
    # circ2=random_circuit(5,10,2)
    # circ2.draw('mpl')
    # plt.show()