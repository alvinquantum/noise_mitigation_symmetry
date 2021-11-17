#!/usr/bin/env python3

# Don't change the following two lines. This gives a performance boost and lets it run on Bebop.
import os
os.environ['OPENBLAS_NUM_THREADS'] = 'unlimited'

from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, execute, Aer, transpile
from matplotlib import pyplot as plt
from qiskit.quantum_info import Statevector, state_fidelity, DensityMatrix, partial_trace, Operator, Pauli
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit.visualization import circuit_drawer
from qiskit.circuit import QuantumRegister
import pickle, mymodule, glob, sys, json, psutil, random, time
import qiskit.providers.aer.noise as noise
from qiskit.opflow import X,Y,Z
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
from math import pi
from multiprocessing import Pool, Value, Array, Manager
from os import listdir
from os.path import isfile

def create_controlU(npmat, number_of_qubits):
    '''Returns a controlled operation. Type is np.array'''
    return np.kron(np.array([[1,0],[0,0]]),np.eye(2**number_of_qubits))+np.kron(np.array([[0,0],[0,1]]),npmat)

def check_p2(control_p1, control_p2, unitary, number_of_qubits):
    '''Sanity check for p2. U\otimes I- ControlP2^\dagger(U\otimes I)ControlP1==0'''
    assert np.allclose(np.kron(np.eye(2),unitary.data)-control_p2.dot(np.kron(np.eye(2),unitary.data)).dot(control_p1), np.zeros(2**(number_of_qubits+1))), "wrong p2"

def add_controlU(circ, pauli_str, number_of_qubits, quantum_register, ancilla_register):
    '''Adds a controlled Pauli to circ.'''
    # print(pauli_str)
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

    #Check the length of the pauli string. It maybe empty when there were no matches found.
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
            else:
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
            qubit_pos-=1
    return circ

def result_exists(base_path, file_name, result_num):
    '''Checks if the result file exists.'''
    name_split=file_name.split("_")
    result_name="_".join(name_split[:-1])+"_result_"+str(result_num)+"_.txt"
    # print(file_name)
    # print(result_name)
    output_file_txt_path=os.path.join(base_path, result_name)
    # print(output_file_txt_path)
    return os.path.isfile(output_file_txt_path)

def get_results(error_tuple):
    '''Main function for pool'''
    global NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files, one_qubit_error_space, file_info, file_number, circ
    global max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, quantum_register, ancilla_register, all_qubits
    global circ_full, circ_no_checks, final_state_no_errors_with_checks, final_state_no_errors_no_checks, reduced_correct_state, file_name
    global circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo
    error_idx=error_tuple[0]
    one_qubit_err=error_tuple[1]
    two_qubit_err = one_qubit_err*10

    # Depolarizing quantum errors
    error_1 = noise.depolarizing_error(one_qubit_err, 1)
    error_2 = noise.depolarizing_error(two_qubit_err, 2)

    # Add errors to noise model
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_1, ['id','rz','sx'])
    noise_model.add_all_qubit_quantum_error(error_2, ['cx'])
    basis_gates=noise_model.basis_gates
    # simulator=Aer.get_backend('qasm_simulator')
    # simulator.set_option("method", "statevector")
    # # No checks tomography
    # circ_no_checks_tomo=state_tomography_circuits(circ_no_checks, quantum_register)
    # # With Checks tomography with classical register
    # circ_full_no_ancilla_tomo=state_tomography_circuits(circ_full, quantum_register)
    # circ_full_with_ancilla_tomo=deepcopy(circ_full_no_ancilla_tomo)
    # classical_register=ClassicalRegister(1)
    # for circ_temp in circ_full_with_ancilla_tomo:
    #     circ_temp.add_register(classical_register)
    #     circ_temp.measure(ancilla_register[0], classical_register)
    # print(circ_no_checks[0].decompose())
    # print(circ_full_no_ancilla[0].decompose())
    # print(circ_full_with_ancilla_tomo[0].decompose())

    #No checks
    result_no_checks=run_job(circ_no_checks_tomo, basis_gates, noise_model, 10000)
    tomo_fitter_no_checks = StateTomographyFitter(result_no_checks, circ_no_checks_tomo)
    rho_fit_no_checks = tomo_fitter_no_checks.fit(method='lstsq')

    # With Checks
    result_with_checks=run_job(circ_full_with_ancilla_tomo, basis_gates, noise_model, 10000)
    result_with_checks=mymodule.post_select_on_ancilla(result_with_checks, "0", NUMBER_OF_QUBITS)
    tomo_fitter_checks = StateTomographyFitter(result_with_checks, circ_full_no_ancilla_tomo)
    rho_fit_checks = tomo_fitter_checks.fit(method='lstsq')

    # outputs
    state_fidelity_no_checks_with_errors=state_fidelity(DensityMatrix(rho_fit_no_checks), reduced_correct_state)
    state_fidelity_with_checks_with_errors= state_fidelity(DensityMatrix(rho_fit_checks), reduced_correct_state)
    state_fidelity_with_checks_no_errors=state_fidelity(final_state_no_errors_with_checks, final_state_no_errors_no_checks)
    assert state_fidelity_with_checks_no_errors>0.98, "Sanity check fidelity "+ str(state_fidelity_with_checks_no_errors) +" failed for "+ file_name

    print("Error_idx", error_idx)
    print("One qubit error: "+str(one_qubit_err))
    print("One qubit error: "+str(two_qubit_err))
    print('State fidelity no checks with errors: ', state_fidelity_no_checks_with_errors)
    print('State fidelity with checks and errors: ', state_fidelity_with_checks_with_errors)
    # Sanity check
    print("Sanity check fidelity with checks no errors: ", state_fidelity_with_checks_no_errors)
    print()
    # print("One_qubit_err:", one_qubit_err)
    # print(circ.count_ops())

    return {"error_idx": error_idx, "one_qubit_err": one_qubit_err, "two_qubit_err": two_qubit_err, "state_fidelity_no_checks_with_errors": state_fidelity_no_checks_with_errors, 
    "state_fidelity_with_checks_with_errors": state_fidelity_with_checks_with_errors, "state_fidelity_with_checks_no_errors": state_fidelity_with_checks_no_errors}

def run_job(job, basis_gates, noise_model, shots_arg):
    '''Runs the given job. Returns result.'''
    print("running job...")
    result= execute(job, Aer.get_backend('qasm_simulator'),
                basis_gates=basis_gates,
                noise_model=noise_model, 
                shots=shots_arg, optimization_level=0).result()
    print("job done.")
    
    return result

def initialize(NUMBER_OF_QUBITS_ARGS, subdir_arg, code_dir_arg, rand_circ_files_arg, one_qubit_error_space_arg, file_info_arg, file_number_arg, circ_arg,
    max_pauli_weight_arg, max_pauli_str_p1_arg, max_pauli_str_p2_arg, quantum_register_arg, ancilla_register_arg, all_qubits_arg,
    circ_full_arg, circ_no_checks_arg, final_state_no_errors_with_checks_arg, final_state_no_errors_no_checks_arg, reduced_correct_state_arg, file_name_arg,
    circ_no_checks_tomo_arg, circ_full_no_ancilla_tomo_arg, circ_full_with_ancilla_tomo_arg):
    '''Initialize globals for pool'''
    global NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files, one_qubit_error_space, file_info, file_number, circ
    global max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, quantum_register, ancilla_register, all_qubits
    global circ_full, circ_no_checks, final_state_no_errors_with_checks, final_state_no_errors_no_checks, reduced_correct_state, file_name
    global circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo
    NUMBER_OF_QUBITS=NUMBER_OF_QUBITS_ARGS
    subdir= subdir_arg
    code_dir=code_dir_arg
    rand_circ_files= rand_circ_files_arg
    one_qubit_error_space=one_qubit_error_space_arg
    file_info=file_info_arg 
    file_number=file_number_arg
    max_pauli_weight=max_pauli_weight_arg 
    max_pauli_str_p1=max_pauli_str_p1_arg 
    max_pauli_str_p2=max_pauli_str_p2_arg
    quantum_register=quantum_register_arg 
    ancilla_register=ancilla_register_arg
    all_qubits=all_qubits_arg
    circ_full=circ_full_arg
    circ_no_checks=circ_no_checks_arg
    final_state_no_errors_with_checks=final_state_no_errors_with_checks_arg
    final_state_no_errors_no_checks=final_state_no_errors_no_checks_arg
    reduced_correct_state=reduced_correct_state_arg
    circ=circ_arg
    file_name=file_name_arg
    circ_no_checks_tomo=circ_no_checks_tomo_arg
    circ_full_no_ancilla_tomo=circ_full_no_ancilla_tomo_arg
    circ_full_with_ancilla_tomo=circ_full_with_ancilla_tomo_arg

def store_results(circ_full_with_ancilla, base_path, circ_file_name, number_of_qubits, results, cnot_count, rz_count, circuit_num,
    count, p2_weight, pauli_str_p1, pauli_str_p2):
    '''Stores all the results.'''
    #File naming stuff. 
    #Strip the extension.
    split_circ_file_name=circ_file_name.split("_")
    file_name_no_extension="_".join(split_circ_file_name[:-1])

    # print(output_file_name_obj)
    # print(os.path.join(base_path,output_file_name_obj))
    temp_file_number=0
    output_file_name_obj=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.obj"
    while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
        temp_file_number+=1
        output_file_name_obj=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.obj"
    output_file_name_txt=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.txt"
    output_file_name_qasm=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.qasm"
    circ_full_with_ancilla.qasm(filename=os.path.join(base_path, output_file_name_qasm))

    # Dump all the results into a pickle
    with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
        pickle.dump({"cx": cnot_count, "rz": rz_count, "qubits": number_of_qubits, "circuit_num" : circuit_num, "found_matches: ": count,
             "max_pauli_weight": p2_weight, "max_pauli_str_p1": pauli_str_p1, "max_pauli_str_p2": pauli_str_p2 ,"results": results}, circ_file)

    #Print text results to file
    output_file_txt_path=os.path.join(base_path, output_file_name_txt)
    circuit_drawer(circ_full_with_ancilla.decompose(), filename=output_file_txt_path)
    with open(output_file_txt_path, "a") as output_file_txt:
        output_file_txt.write("\n")
        output_file_txt.write(json.dumps(circ.count_ops()))
        for result in results:
            output_file_txt.write("\n")
            output_file_txt.write("Error idx: "+str(result["error_idx"])+"\n")
            output_file_txt.write("One qubit error: "+str(result["one_qubit_err"])+"\n")
            # print("One_qubit_err:", one_qubit_err)
            output_file_txt.write("Two qubit error: "+str(result["two_qubit_err"])+"\n")
            # print("Two qubit error:", two_qubit_err)
            output_file_txt.write("State fidelity no checks and with errors: "+str(result["state_fidelity_no_checks_with_errors"])+"\n")
            output_file_txt.write("State fidelity with checks and with errors: "+str(result["state_fidelity_with_checks_with_errors"])+"\n")
            output_file_txt.write("Sanity check fidelity with checks and no errors: "+str(result["state_fidelity_with_checks_no_errors"])+"\n")

if __name__ == "__main__":
    print("running...")
    #Disable parallel. Qiskit aer uses parallel inside which clashes with the outer pool when running in Linux. note: FALSE means enabled.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    #Program parameters. These are globals in the parallel program passed from shell.
    # NUMBER_OF_QUBITS=int(sys.argv[1])
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    start_circ_number=int(sys.argv[3])
    end_circ_number=int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=True
    #Error space
    one_qubit_error_space=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    # one_qubit_error_space=np.logspace(-3, 0, num=21)

    # File stuff
    subdir="/data/"
    # Gets the file path of the script
    code_dir=sys.path[0]
    BASE_PATH=code_dir+subdir
    # Gets the files that match the string. Files include the path
    # rand_circ_files=[elem for elem in glob.glob(code_dir+subdir+"qubits"+str(NUMBER_OF_QUBITS)+"_depth"+ str(DEPTH)+"*.obj") if "result" not in elem]
    all_files=[f for f in listdir(BASE_PATH) if isfile(os.path.join(BASE_PATH, f))]
    rand_circ_files=[]
    circ_properties=[]
    for file in all_files:
        name_split=file.split("_")
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if ".qasm" in name_split and "result" not in name_split and name_split_nums[1]==CNOT_COUNT and start_circ_number<=name_split_nums[2]<=end_circ_number:
            rand_circ_files.append(file)
            circ_properties.append("_".join(name_split[:-1])+"_.obj")
    # print(rand_circ_files)
    # print(circ_properties)

    #Get the qasm and pickle info
    for idx, file_name in enumerate(circ_properties):
        #Stopwatch
        time0=time.time()
        #get the file information of the circuit. Len 3 list: qubits#, cnot#, circuit_number#.
        file_info=[int(elem) for elem in file_name.split("_") if elem.isdigit()]
        assert len(file_info) == 3, "This file doesn't have three numbers in the name."
        file_number=file_info[2]
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        if result_exists(BASE_PATH, file_name, 0):
            continue

        with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
            circ_info=pickle.load(circ_file)
            # circ=circ_info["circ"]
        max_pauli_weight=circ_info["max_pauli_weight"]
        max_pauli_str_p1=circ_info["max_pauli_str_p1"]
        max_pauli_str_p2=circ_info["max_pauli_str_p2"]
        circ=QuantumCircuit.from_qasm_file(os.path.join(BASE_PATH, rand_circ_files[idx]))
        # If there were no Pauli checks found just skip.
        if max_pauli_weight == 0: 
            continue

        unitary=Operator(circ)
        control_p1=create_controlU(Pauli(max_pauli_str_p1).to_matrix(), NUMBER_OF_QUBITS)
        control_p2=create_controlU(Pauli(max_pauli_str_p2).to_matrix(), NUMBER_OF_QUBITS)
        check_p2(control_p1,control_p2, unitary, NUMBER_OF_QUBITS)
        print(np.allclose(unitary.dot(Pauli(max_pauli_str_p1).to_matrix()).data, Pauli(max_pauli_str_p2).to_matrix().dot(unitary)))
        print(unitary.dot(Pauli(max_pauli_str_p1).to_matrix()).data-Pauli(max_pauli_str_p2).to_matrix().dot(unitary))
        print(max_pauli_str_p1)
        print(max_pauli_str_p2)
    #     #Random data qubits
    #     # rand_data_state=random_statevector(2**(NUMBER_OF_QUBITS))
    #     # zero_ancilla_qubit=Statevector.from_int(0, 2)
    #     # full_intial_state=rand_data_state.expand(zero_ancilla_qubit)
    #     # print(rand_data_state)
    #     # print(full_intial_state)

    #     #Initial zero state. The randomness is introduced with the circuit.
    #     zero_state=Statevector.from_int(0, 2**(NUMBER_OF_QUBITS+1))

    #     #Create the complete circuit without a measurement
    #     quantum_register=QuantumRegister(NUMBER_OF_QUBITS)
    #     ancilla_register=QuantumRegister(1)
    #     all_qubits=[]
    #     for reg in quantum_register:
    #         all_qubits.append(reg)
    #     all_qubits.append(ancilla_register[0])
    #     circ_full=QuantumCircuit(quantum_register, ancilla_register)
    #     #Create the circuit with no checks
    #     circ_no_checks=QuantumCircuit(quantum_register, ancilla_register)
    #     #Insert random state generator
    #     random_params = np.random.uniform(size=(NUMBER_OF_QUBITS, 3))
    #     for i,qreg in enumerate(quantum_register):
    #         for j, pauli in enumerate([X,Y,Z]):
    #             rand_tuple=(random_params[i][j])
    #             circ_full = circ_full.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg])
    #             circ_no_checks=circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg])

    #     # Finish the full circuit
    #     circ_full.h(ancilla_register[0])
    #     circ_full.barrier(all_qubits)
    #     circ_full=add_controlU(circ_full, max_pauli_str_p1, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
    #     circ_full.barrier(all_qubits)
    #     circ_full.append(circ, qargs=quantum_register)
    #     circ_full.barrier(all_qubits)
    #     circ_full=add_controlU(circ_full, max_pauli_str_p2, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
    #     circ_full.h(ancilla_register[0])
    #     # print(circ_full.decompose())

    #     #Finish the circuit with no checks
    #     circ_no_checks.append(circ, qargs=quantum_register)
    #     # print(circ_no_checks.decompose())

    #     # The final state. Note that evolve doesn't change the object.
    #     final_state_no_errors_with_checks=zero_state.evolve(circ_full)
    #     final_state_no_errors_no_checks=zero_state.evolve(circ_no_checks)
    #     # print("No errors full state no checks: ", final_state_no_errors_no_checks)
    #     # print("No errors full state with checks: ", final_state_no_errors_checks)
    #     #Reduced final since we will be conducting tomography excluding the ancilla.
    #     reduced_correct_state=partial_trace(final_state_no_errors_no_checks, [NUMBER_OF_QUBITS])

    #     # Tomography stuff.
    #     # No checks tomography
    #     circ_no_checks_tomo=state_tomography_circuits(circ_no_checks, quantum_register)
    #     # With Checks tomography with classical register
    #     circ_full_no_ancilla_tomo=state_tomography_circuits(circ_full, quantum_register)
    #     circ_full_with_ancilla_tomo=deepcopy(circ_full_no_ancilla_tomo)
    #     classical_register=ClassicalRegister(1)
    #     for circ_temp in circ_full_with_ancilla_tomo:
    #         circ_temp.add_register(classical_register)
    #         circ_temp.measure(ancilla_register[0], classical_register)
    #     print(circ_full.decompose())
    #     print(circ.count_ops())
    #     # Parallel code.
    #     if PARALLEL:
    #         with Pool(psutil.cpu_count(logical=False), initialize, initargs=(NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files,
    #             one_qubit_error_space, file_info, file_number, circ, max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, 
    #             quantum_register, ancilla_register, all_qubits, circ_full, circ_no_checks, final_state_no_errors_with_checks,
    #             final_state_no_errors_no_checks, reduced_correct_state, file_name, circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo)) as pool:
    #             results=list(pool.imap_unordered(get_results, enumerate(one_qubit_error_space)))
    #     else:
    #         # Non parallel.
    #         results=[]
    #         for error_idx, one_qubit_err in enumerate(one_qubit_error_space):
    #             results.append(get_results((error_idx, one_qubit_err)))
    #     #Store results.
    #     store_results(circ_full, BASE_PATH, file_name, NUMBER_OF_QUBITS, results, CNOT_COUNT, circ_info["rz"], circ_info["circuit_num"], 
    #         circ_info["found_matches"], max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2)

    #     print("execution time", time.time()-time0)
    # print("Finished.")
    