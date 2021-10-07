#!/usr/bin/env python3

from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, execute, Aer, transpile
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector, state_fidelity, DensityMatrix, partial_trace
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
from math import pi
# from qiskit.circuit import ControlledGate
# from qiskit.extensions import UnitaryGate

#Error Stuff
from qiskit.circuit.library import CRXGate
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
    c_minus_x = QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(1.0592554e-13,3.288277,2.8116442) q[0];
    u3(pi,5.8991791,2.7575864) q[1];
    cx q[0],q[1];
    u3(2*pi,5.4995268,4.1085153) q[0];
    u3(pi,0.38400625,3.5255989) q[1];
    """.strip())

    c_minus_j_x = QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(6.9759301e-13,4.3371506,0.076820611) q[0];
    u3(6.9089241e-09,5.477472,3.947306) q[1];
    cx q[0],q[1];
    u3(2*pi,0.60772643,2.832284) q[0];
    u3(5.5644668e-09,5.4477338,3.9770441) q[1];
    """.strip())

    c_plus_j_x =QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(pi,4.6477413,6.258122) q[0];
    u3(1.7599424,3*pi/2,3*pi/2) q[1];
    cx q[0],q[1];
    u3(pi,3.1666559,3.2062404) q[0];
    u3(4.9015351,3*pi/2,3*pi/2) q[1];
    """.strip())

    if len(pauli_str)>number_of_qubits:
        phase=pauli_str[:2]
        pauli_str=pauli_str[2:]
        qubit_pos=number_of_qubits-1
        phase_added=False
        for pos, char in enumerate(pauli_str):
            # First case is most common
            if phase_added or phase=="+1":
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
            # -1 phase
            elif phase=="-1":
                # circ.cz(ancilla_register, quantum_register[qubit_pos])
                if char=="X":
                    circ.compose(c_minus_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                elif char=="Y":
                    circ.sdg(quantum_register[qubit_pos])
                    circ.compose(c_minus_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    # print(RXGate(pi*j).to_matrix())
                    circ.h(quantum_register[qubit_pos])
                    circ.compose(c_minus_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    # circ.crx(-j*pi,ancilla_register, quantum_register[qubit_pos])
                    circ.h(quantum_register[qubit_pos])
                # circ.cz(ancilla_register, quantum_register[qubit_pos])
                phase_added=True
            # +j phase
            elif phase=="+j":
                if char=="X":
                    circ.compose(c_plus_j_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                elif char=="Y":
                    circ.sdg(quantum_register[qubit_pos])
                    circ.compose(c_plus_j_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    circ.h(quantum_register[qubit_pos])
                    circ.compose(c_plus_j_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.h(quantum_register[qubit_pos])
                phase_added=True
            # -j phase
            elif phase=="-j":
                if char=="X":
                    circ.compose(c_minus_j_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                elif char=="Y":
                    circ.sdg(quantum_register[qubit_pos])
                    circ.compose(c_minus_j_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    circ.h(quantum_register[qubit_pos])
                    circ.compose(c_minus_j_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.h(quantum_register[qubit_pos])
                phase_added=True
            # print("pos ", pos)
            # print("qubit_pos ", qubit_pos)
            qubit_pos-=1
        # print(pauli_str)
    return circ

#Program parameters
NUMBER_OF_QUBITS=5
DEPTH=10
NUMBER_OF_CIRCUITS=3
#Paths for outputs and pickle file of circuit
code_dir=os.path.dirname(os.path.realpath('__file__'))
file_number="32"
subdir="/data/"
base_file_path=code_dir+subdir+"depth"+ str(DEPTH) +"_"

#Error probabilities
one_qubit_err=0.0025
two_qubit_err=0.025

error_space=np.logspace(0, 1, base=1)

if __name__ == "__main__":
    # print(error_space)
    #Get the pickle info
    circ_file=open(base_file_path+ file_number + ".obj", "rb")
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
    circ_full.barrier(all_qubits)
    circ_full=add_controlU(circ_full, max_pauli_str_p1, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
    circ_full.barrier(all_qubits)
    circ_full.append(circ, qargs=quantum_register)
    circ_full.barrier(all_qubits)
    circ_full=add_controlU(circ_full, max_pauli_str_p2, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
    circ_full.h(ancilla_register[0])
    # print(circ_full.decompose())

    #Finish the circuit with no checks
    circ_no_checks.append(circ, qargs=quantum_register)
    # print(circ_no_checks.decompose())

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
    #Reduced final since we will be conducting tomography excluding the ancilla.
    reduced_correct_state=partial_trace(final_state_no_errors_checks, [5])

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
    circ_no_checks=state_tomography_circuits(circ_no_checks, quantum_register)
    # With Checks tomography with classical register
    circ_full_no_ancilla=state_tomography_circuits(circ_full, quantum_register)
    circ_full_with_ancilla=deepcopy(circ_full_no_ancilla)
    for circ_temp in circ_full_with_ancilla:
        circ_temp.add_register(classical_register)
        circ_temp.measure(ancilla_register[0], classical_register)
    # print(circ_no_checks[0].decompose())
    # print(circ_full_no_ancilla[0].decompose())
    print(circ_full_with_ancilla[0].decompose())

    #No checks
    print("running job...")
    result_no_checks=execute(circ_no_checks, Aer.get_backend('qasm_simulator'),
                 basis_gates=basis_gates,
                 noise_model=noise_model, 
                 shots=10000).result()
    print("job done.")
    tomo_fitter_no_checks = StateTomographyFitter(result_no_checks, circ_no_checks)
    rho_fit_no_checks = tomo_fitter_no_checks.fit(method='lstsq')
    print('State fidelity no checks with errors: ', state_fidelity(DensityMatrix(rho_fit_no_checks), reduced_correct_state))

    # With Checks
    print("running job...")
    result_with_checks=execute(circ_full_with_ancilla, Aer.get_backend('qasm_simulator'),
                 basis_gates=basis_gates,
                 noise_model=noise_model, 
                 shots=10000).result()
    print("job done.")
    # print(result_with_checks.get_counts())
    result_with_checks=mymodule.post_select_on_ancilla(result_with_checks, "0", 5)
    # print(result_with_checks.get_counts())
    tomo_fitter_checks = StateTomographyFitter(result_with_checks, circ_full_no_ancilla)
    rho_fit_checks = tomo_fitter_checks.fit(method='lstsq')
    print('State fidelity with checks and errors', state_fidelity(DensityMatrix(rho_fit_checks), reduced_correct_state))
    # Sanity check
    print("Sanity check fidelity no errors with checks: ", state_fidelity(final_state_no_errors_checks,final_state_no_errors_no_checks))
    print(circ.count_ops())

    #Pickle the results
    temp_file_number=0
    output_file_obj=base_file_path+ file_number+"_"+ str(temp_file_number) +"_result.obj"
    while os.path.isfile(output_file_obj):
        temp_file_number+=1
        output_file_obj=base_file_path+ file_number+ "_" + str(temp_file_number) +"_result.obj"
    output_file_txt_path=base_file_path+ file_number+ "_" + str(temp_file_number) +"_result.txt"
    # Dump all the info into a pickle
    circ_file=open(output_file_obj, "wb")
    pickle.dump({"circ": circ, "circ_no_checks": circ_no_checks, "circ_full_no_ancilla": circ_full_no_ancilla, 
    "circ_full_with_ancilla": circ_full_with_ancilla, "rho_fit_checks": DensityMatrix(rho_fit_checks), 
    "rho_fit_no_checks":DensityMatrix(rho_fit_no_checks), "reduced_correct_state": reduced_correct_state,
    "max_pauli_weight": max_pauli_weight, "max_pauli_str_p1": max_pauli_str_p1, "max_pauli_str_p2": max_pauli_str_p2}, circ_file)
    circ_file.close()

    #Print text results to file
    circ.draw('mpl')
    circuit_drawer(circ_full_with_ancilla[0].decompose(), filename=output_file_txt_path)
    output_file_txt=open(output_file_txt_path, "a")
    print(file=output_file_txt)
    print('State fidelity no checks with errors: ', state_fidelity(DensityMatrix(rho_fit_no_checks), reduced_correct_state), file=output_file_txt)
    print('State fidelity with checks and errors', state_fidelity(DensityMatrix(rho_fit_checks), reduced_correct_state), file=output_file_txt)
    print("Sanity check fidelity no errors with checks: ", state_fidelity(final_state_no_errors_checks,final_state_no_errors_no_checks), file=output_file_txt)
    print(circ.count_ops(), file=output_file_txt)
    output_file_txt.close()