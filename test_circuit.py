#!/usr/bin/env python3

from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, execute, Aer, transpile
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector, state_fidelity
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit.visualization import circuit_drawer
from qiskit.circuit import QuantumRegister
import os
import pickle
import qiskit.providers.aer.noise as noise
from qiskit.opflow import X,Y,Z
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
import mymodule
# from qiskit.circuit import ControlledGate
# from qiskit.extensions import UnitaryGate

#Error Stuff

from qiskit.circuit.library.standard_gates import (IGate, U1Gate, U2Gate, U3Gate, XGate,
                                                   YGate, ZGate, HGate, SGate, SdgGate, TGate,
                                                   TdgGate, RXGate, RYGate, RZGate, CXGate,
                                                   CYGate, CZGate, CHGate, CRZGate, CU1Gate,
                                                   CU3Gate, SwapGate, RZZGate,
                                                   CCXGate, CSwapGate)

def create_controlU(npmat, number_of_qubits):
    '''Returns a controlled operation. Type is np.array'''
    return np.kron(np.array([[1,0],[0,0]]),np.eye(2**number_of_qubits))+np.kron(np.array([[0,0],[0,1]]),npmat)

def check_p2(control_p1, control_p2, unitary, number_of_qubits):
    '''Sanity check for p2. U\otimes I- ControlP2^\dagger(U\otimes I)ControlP1==0'''
    assert np.allclose(np.kron(np.eye(2),unitary.data)-control_p2.dot(np.kron(np.eye(2),unitary.data)).dot(control_p1), np.zeros(2**(number_of_qubits+1))), "wrong p2"

def add_controlU(circ, pauli_str, number_of_qubits, quantum_register, ancilla_register):
    '''Adds a controlled Pauli to circ.'''
    print(pauli_str)
    if len(pauli_str)>number_of_qubits:
        phase=pauli_str[:2]
        pauli_str=pauli_str[2:]
        qubit_pos=0
        for pos, char in enumerate(pauli_str):
            if char=="X":
                circ.cx(ancilla_register, quantum_register[qubit_pos])
            elif char=="Y":
                circ.sdg(quantum_register[qubit_pos])
                circ.cx(ancilla_register, quantum_register[qubit_pos])
                circ.s(quantum_register[qubit_pos])
            elif char=="Z":
                circ.h(quantum_register[qubit_pos])
                circ.cx(ancilla_register, quantum_register[qubit_pos])
                circ.h(quantum_register[qubit_pos])
            qubit_pos+=1
        # print(pauli_str)
    return circ

#Program parameters
NUMBER_OF_QUBITS=5
DEPTH=10
NUMBER_OF_CIRCUITS=3
#Paths for outputs and pickle file of circuit
code_dir=os.path.dirname(os.path.realpath('__file__'))
subdir="/data/"
base_file_path=code_dir+subdir+"depth"+ str(DEPTH) +"_"

#Error probabilities
one_qubit_err=0.0025
two_qubit_err=0.025



if __name__ == "__main__":
    #Get the pickle info
    circ_file=open(base_file_path+"8.obj", "rb")
    circ_info=pickle.load(circ_file)
    circ=circ_info["circ"]
    max_pauli_weight=circ_info["max_pauli_weight"]
    max_pauli_str_p1=circ_info["max_pauli_str_p1"]
    max_pauli_str_p2=circ_info["max_pauli_str_p2"]

    #Random data qubits
    # rand_data_state=random_statevector(2**(NUMBER_OF_QUBITS))
    # zero_ancilla_qubit=Statevector.from_int(0, 2)
    # full_intial_state=rand_data_state.expand(zero_ancilla_qubit)
    # print(rand_data_state)
    # print(full_intial_state)

    #Initial zero state. The randomness is introduced with the circuit.
    zero_state=Statevector.from_int(0, 2**(NUMBER_OF_QUBITS+1))

    #Create the complete circuit without a measurement
    quantum_register=QuantumRegister(5)
    ancilla_register=QuantumRegister(1)
    classical_register=ClassicalRegister(1)
    all_qubits=[]
    for reg in quantum_register:
        all_qubits.append(reg)
    all_qubits.append(ancilla_register[0])
    circ_full=QuantumCircuit(quantum_register, ancilla_register)
    #Create the circuit with no checks
    circ_no_checks=QuantumCircuit(quantum_register, ancilla_register)
    #Insert random state generator
    random_params = np.random.uniform(size=(NUMBER_OF_QUBITS, 3))
    for i,qreg in enumerate(quantum_register):
        for j, pauli in enumerate([X,Y,Z]):
            rand_tuple=(random_params[i][j])
            circ_full = circ_full.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg])
            circ_no_checks=circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg])

    # Finish the full circuit
    circ_full.h(ancilla_register[0])
    # p1=Pauli(max_pauli_str_p1)
    # p2=Pauli(max_pauli_str_p2)
    circ_full=add_controlU(circ_full, max_pauli_str_p1, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
    circ_full.barrier(all_qubits)
    circ_full.append(circ, qargs=quantum_register)
    circ_full.barrier(all_qubits)
    circ_full=add_controlU(circ_full, max_pauli_str_p2, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
    circ_full.h(ancilla_register[0])
    print(circ_full.decompose())

    #Finish the circuit with no checks
    circ_no_checks.append(circ, qargs=[0,1,2,3,4])
    print(circ_no_checks.decompose())

    # zero_state=Statevector.from_int(0, 2**(NUMBER_OF_QUBITS+1))

    #Alvin: don't delete...tests
    # testqc=QuantumCircuit(2)
    # testqc.x([1])
    # test_st=Statevector(testqc)
    # print(testqc)
    # print(test_st)
    # print("measure 1", test_st.measure([1]))
    # print("measure 0", test_st.measure([0]))
    # print(test_st)

    # Print the final state. Note that evolve doesn't change the object.
    final_state_no_errors_checks=zero_state.evolve(circ_full)
    final_state_no_errors_no_checks=zero_state.evolve(circ_no_checks)
    # print("No errors full state no checks: ", final_state_no_errors_no_checks)
    # print("No errors full state with checks: ", final_state_no_errors_checks)
    print("Fidelity no errors with checks: ", state_fidelity(final_state_no_errors_checks,final_state_no_errors_no_checks))

    #Add measurements to the cricuits.
    # circ_full.measure(ancilla_register[0], classical_register)

    # Depolarizing quantum errors
    error_1 = noise.depolarizing_error(one_qubit_err, 1)
    error_2 = noise.depolarizing_error(two_qubit_err, 2)

    # Add errors to noise model
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_1, ['h', 's', 't', 'u3', 'id'])
    noise_model.add_all_qubit_quantum_error(error_2, ['cx'])
    basis_gates=noise_model.basis_gates
    simulator=Aer.get_backend('qasm_simulator')
    # No checks tomography
    circ_no_checks=state_tomography_circuits(circ_no_checks, all_qubits)
    # With Checks tomography with classical register
    circ_full_no_ancilla=state_tomography_circuits(circ_full, quantum_register)
    circ_full_with_ancilla=deepcopy(circ_full_no_ancilla)
    for circ_temp in circ_full_with_ancilla:
        circ_temp.add_register(classical_register)
        circ_temp.measure(ancilla_register[0], classical_register)
    print(circ_full_with_ancilla[0])

    #No checks
    # print("running job...")
    # result_no_checks=execute(circ_no_checks, Aer.get_backend('qasm_simulator'),
    #              basis_gates=basis_gates,
    #              noise_model=noise_model, shots=5).result()
    # print("job done.")
    # tomo_fitter_no_checks = StateTomographyFitter(result_no_checks, circ_no_checks)
    # rho_fit_no_checks = tomo_fitter_no_checks.fit(method='lstsq')
    # print('State Fidelity: ', state_fidelity(rho_fit_no_checks, final_state_no_errors_no_checks))

    # With Checks
    print("running job...")
    result_with_checks=execute(circ_full_with_ancilla, Aer.get_backend('qasm_simulator'),
                 basis_gates=basis_gates,
                 noise_model=noise_model, shots=10).result()
    print("job done.")
    print(result_with_checks.results[0])
    result_with_checks=mymodule.post_select_on_ancilla(result_with_checks, "1", 5)
    tomo_fitter_checks = StateTomographyFitter(result_with_checks, circ_full_no_ancilla)
    rho_fit_checks = tomo_fitter_checks.fit(method='lstsq')
    print('State Fidelity: ', state_fidelity(rho_fit_checks, final_state_no_errors_no_checks))