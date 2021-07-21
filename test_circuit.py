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
from qiskit.circuit import QuantumRegister
import os
import pickle
# from qiskit.circuit import ControlledGate
# from qiskit.extensions import UnitaryGate

from qiskit.circuit.library.standard_gates import (IGate, U1Gate, U2Gate, U3Gate, XGate,
                                                   YGate, ZGate, HGate, SGate, SdgGate, TGate,
                                                   TdgGate, RXGate, RYGate, RZGate, CXGate,
                                                   CYGate, CZGate, CHGate, CRZGate, CU1Gate,
                                                   CU3Gate, SwapGate, RZZGate,
                                                   CCXGate, CSwapGate)

#Program parameters
NUMBER_OF_QUBITS=5
DEPTH=40
NUMBER_OF_CIRCUITS=3
#Paths for outputs and pickle file of circuit
code_dir=os.path.dirname(os.path.realpath('__file__'))
subdir="/data/"
base_file_path=code_dir+subdir+"depth"+ str(DEPTH) +"_"

if __name__ == "__main__":
    # print(base_file_path)
    circ_file=open(base_file_path+"1.obj", "rb")
    circ_info=pickle.load(circ_file)
    circ=circ_info["circ"]
    ancilla=QuantumRegister(size=1, name="q_a")
    circ.add_register(ancilla)
    print(circ)
    print(circ_info["max_pauli_weight"])
    print(circ_info["max_pauli_str_p1"])
    print(circ_info["max_pauli_str_p2"])

    zero_state=Statevector.from_int(0, 2**(NUMBER_OF_QUBITS+1))
    # rand_state=zero_state.evolve(RXGate)