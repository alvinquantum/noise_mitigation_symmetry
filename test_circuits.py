#!/usr/bin/env python3

# Don't change the following two lines. This gives a performance boost and lets it run on Bebop.
import os

from qiskit.ignis.verification import tomography
os.environ['OPENBLAS_NUM_THREADS'] = 'unlimited'

from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, execute, Aer, result, transpile
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
from math import e, pi
from multiprocessing import Pool, Value, Array, Manager
from os import listdir
from os.path import isfile
import math
# import cProfile, io
# from pstats import SortKey
import line_profiler, kernprof
import cirq

# @profile
def get_results(error_tuple):
    '''Main function for pool'''
    global NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files, one_qubit_error_space, file_info, file_number, circ
    global max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, quantum_register, ancilla_register, all_qubits
    global circ_full, circ_no_checks, final_state_no_errors_with_checks, final_state_no_errors_no_checks, reduced_correct_state, file_name
    global circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo, cpus_pool, do_tomography#, simulator
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

    simulator=Aer.get_backend('qasm_simulator')
    simulator.set_option("method", "density_matrix")
    # simulator.set_option("max_parallel_experiments", math.floor(psutil.cpu_count(logical=False)/cpus_pool))
    simulator.set_option("max_parallel_experiments", psutil.cpu_count(logical=False))

    #No checks
    result_no_checks=mymodule.run_job(circ_no_checks_tomo, basis_gates, noise_model, 10000, simulator)
    tomo_fitter_no_checks = StateTomographyFitter(result_no_checks, circ_no_checks_tomo)
    rho_fit_no_checks = tomo_fitter_no_checks.fit(method='lstsq')

    # With Checks
    result_with_checks=mymodule.run_job(circ_full_with_ancilla_tomo, basis_gates, noise_model, 10000, simulator)
    #Test
    num_results_before_postselect=mymodule.get_number_of_results(result_with_checks)

    result_with_checks=mymodule.post_select_on_ancilla(result_with_checks, "0", NUMBER_OF_QUBITS)
    #Test
    num_results_after_postselect=mymodule.get_number_of_results(result_with_checks)
    

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

    return {"num_results_before_postselect": num_results_before_postselect, "num_results_after_postselect": num_results_after_postselect, "error_idx": error_idx, "one_qubit_err": one_qubit_err, "two_qubit_err": two_qubit_err, "state_fidelity_no_checks_with_errors": state_fidelity_no_checks_with_errors, 
    "state_fidelity_with_checks_with_errors": state_fidelity_with_checks_with_errors, "state_fidelity_with_checks_no_errors": state_fidelity_with_checks_no_errors}

def initialize(NUMBER_OF_QUBITS_ARGS, subdir_arg, code_dir_arg, rand_circ_files_arg, one_qubit_error_space_arg, file_info_arg, file_number_arg, circ_arg,
    max_pauli_weight_arg, max_pauli_str_p1_arg, max_pauli_str_p2_arg, quantum_register_arg, ancilla_register_arg, all_qubits_arg,
    circ_full_arg, circ_no_checks_arg, final_state_no_errors_with_checks_arg, final_state_no_errors_no_checks_arg, reduced_correct_state_arg, file_name_arg,
    circ_no_checks_tomo_arg, circ_full_no_ancilla_tomo_arg, circ_full_with_ancilla_tomo_arg, cpus_pool_arg, do_tomography_arg):
    '''Initialize globals for pool'''
    global NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files, one_qubit_error_space, file_info, file_number, circ
    global max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, quantum_register, ancilla_register, all_qubits
    global circ_full, circ_no_checks, final_state_no_errors_with_checks, final_state_no_errors_no_checks, reduced_correct_state, file_name
    global circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo, cpus_pool#, simulator
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
    cpus_pool=cpus_pool_arg
    do_tomography=do_tomography_arg
    #simulator=simulator_arg

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
    circuit_drawer(circ_full_with_ancilla, filename=output_file_txt_path)
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

def get_files(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT, start_circ_number, end_circ_number):
    '''Get the desired files for testing.'''

    # Gets the files that match the string. Files include the path
    # rand_circ_files=[elem for elem in glob.glob(code_dir+subdir+"qubits"+str(NUMBER_OF_QUBITS)+"_depth"+ str(DEPTH)+"*.obj") if "result" not in elem]
    all_files=[f for f in listdir(BASE_PATH) if isfile(os.path.join(BASE_PATH, f))]
    rand_circ_files=[]
    circ_properties_files=[]
    for file in all_files:
        name_split=file.split("_")
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if (".qasm" in name_split and "result" not in name_split and name_split_nums[1]==CNOT_COUNT and name_split_nums[0]==NUMBER_OF_QUBITS 
            and start_circ_number<=name_split_nums[2]<=end_circ_number):
            rand_circ_files.append(file)
            circ_properties_files.append("_".join(name_split[:-1])+"_.obj")

    return (rand_circ_files, circ_properties_files)

def add_rand_input_state(NUMBER_OF_QUBITS, quantum_register, circ_full, circ_no_checks):
    '''Create a random state.'''
    #Insert random state generator
    random_params = np.random.uniform(size=(NUMBER_OF_QUBITS, 3))
    for i,qreg in enumerate(quantum_register):
        for j, pauli in enumerate([X,Y,Z]):
            rand_tuple=(random_params[i][j])
            circ_full.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)
            circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)

def split_circuit_by_barrier(qasm_file_path):
    '''Testing circuits: Split circuits by barrier.'''
    with open(qasm_file_path, "r") as file:
        qasm = file.read()
        prelude = []
        circuits = [[]]
        for line in qasm.splitlines():
            if any([line.startswith(t) for t in ['OPENQASM', 'include', 'qreg', 'creg']]):
                prelude.append(line)
            elif line.startswith('barrier'):
                circuits.append([])
            else:
                circuits[-1].append(line)
        circuits_with_prelude = [prelude+circuit for circuit in circuits]
        return list(map(lambda x: QuantumCircuit.from_qasm_str("\n".join(x)), circuits_with_prelude))

# def sso_results(error_tuple):
#     global circ_full, circ_no_checks, final_state_no_errors_with_checks, final_state_no_errors_no_checks, reduced_correct_state
#     error_idx=error_tuple[0]
#     one_qubit_err=error_tuple[1]
#     two_qubit_err = one_qubit_err*10

#     circ_no_checks.measure_all()
#     circ_full.measure_all()

#     # Depolarizing quantum errors
#     error_1 = noise.depolarizing_error(one_qubit_err, 1)
#     error_2 = noise.depolarizing_error(two_qubit_err, 2)

#     # Add errors to noise model
#     noise_model = noise.NoiseModel()
#     noise_model.add_all_qubit_quantum_error(error_1, ['id','rz','sx'])
#     noise_model.add_all_qubit_quantum_error(error_2, ['cx'])
#     basis_gates=noise_model.basis_gates

#     simulator=Aer.get_backend('qasm_simulator')
#     simulator.set_option("method", "density_matrix")
#     # simulator.set_option("max_parallel_experiments", math.floor(psutil.cpu_count(logical=False)/cpus_pool))
#     simulator.set_option("max_parallel_experiments", psutil.cpu_count(logical=False))
#     simulator.set_option("noise_model", noise_model)
#     result_with_checks=simulator.run(circ_full).result()
#     result_no_checks=simulator.run(circ_no_checks).result()

#     print(result_no_checks)
#     print(result_with_checks)

#     num_results_before_postselect=mymodule.get_number_of_results(result_with_checks)
#     result_with_checks=mymodule.post_select_on_ancilla(result_with_checks, "0", NUMBER_OF_QUBITS)
#     num_results_after_postselect=mymodule.get_number_of_results(result_with_checks)
#     state_fidelity_no_checks_with_errors=get_sso(result_no_checks, final_state_no_errors_no_checks)
#     state_fidelity_with_checks_with_errors=get_sso(result_with_checks, final_state_no_errors_no_checks)
#     state_fidelity_with_checks_no_errors=state_fidelity(final_state_no_errors_with_checks, final_state_no_errors_no_checks)

#     print("Error_idx", error_idx)
#     print("One qubit error: "+str(one_qubit_err))
#     print("One qubit error: "+str(two_qubit_err))
#     print('State fidelity no checks with errors: ', state_fidelity_no_checks_with_errors)
#     print('State fidelity with checks and errors: ', state_fidelity_with_checks_with_errors)
#     # Sanity check
#     print("Sanity check fidelity with checks no errors: ", state_fidelity_with_checks_no_errors)
#     print()
#     # print("One_qubit_err:", one_qubit_err)
#     # print(circ.count_ops())

#     return {"num_results_before_postselect": num_results_before_postselect, "num_results_after_postselect": num_results_after_postselect, "error_idx": error_idx, "one_qubit_err": one_qubit_err, "two_qubit_err": two_qubit_err, "state_fidelity_no_checks_with_errors": state_fidelity_no_checks_with_errors, 
#     "state_fidelity_with_checks_with_errors": state_fidelity_with_checks_with_errors, "state_fidelity_with_checks_no_errors": state_fidelity_with_checks_no_errors}

if __name__ == "__main__":
    # #Test
    # pr=cProfile.Profile()
    # pr.enable()
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
    PARALLEL=False
    do_tomography=True
    #Error space
    one_qubit_error_space=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    # one_qubit_error_space=np.logspace(-3, 0, num=21)
#     one_qubit_error_space=[1.00000000e-05, 1.41253754e-05, 1.99526231e-05, 2.81838293e-05,
#  3.98107171e-05, 5.62341325e-05, 7.94328235e-05, 1.12201845e-04,
#  1.58489319e-04, 2.23872114e-04, 3.16227766e-04, 4.46683592e-04,
#  6.30957344e-04, 8.91250938e-04, 1.25892541e-03, 1.77827941e-03,
#  2.51188643e-03, 3.54813389e-03, 5.01187234e-03, 7.07945784e-03,
#  1.00000000e-02, 0.015, 0.02, 0.025, 0.03, 0.035, 0.06]

    # File stuff
    subdir="/data/"
    # Gets the file path of the script
    # code_dir=sys.path[0]
    code_dir=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=code_dir+subdir
    # # Gets the files that match the string. Files include the path
    # # rand_circ_files=[elem for elem in glob.glob(code_dir+subdir+"qubits"+str(NUMBER_OF_QUBITS)+"_depth"+ str(DEPTH)+"*.obj") if "result" not in elem]
    # all_files=[f for f in listdir(BASE_PATH) if isfile(os.path.join(BASE_PATH, f))]
    files_found=get_files(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT, start_circ_number, end_circ_number)
    rand_circ_files=files_found[0]
    circ_properties_files=files_found[1]

    #Get the qasm and pickle info
    for idx, file_name in enumerate(circ_properties_files):
        #Stopwatch
        time0=time.time()
        # os.environ['QISKIT_IN_PARALLEL'] = 'FALSE'

        #get the file information of the circuit. Len 3 list: qubits#, cnot#, circuit_number#.
        file_info=[int(elem) for elem in file_name.split("_") if elem.isdigit()]
        assert len(file_info) == 3, "This file doesn't have three numbers in the name."
        file_number=file_info[2]
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        if mymodule.result_exists(BASE_PATH, file_name, 0):
            continue

        with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
            circ_info=pickle.load(circ_file)
            # circ=circ_info["circ"]
        max_pauli_weight=circ_info["max_pauli_weight"]
        max_pauli_str_p1=circ_info["max_pauli_str_p1"]
        max_pauli_str_p2=circ_info["max_pauli_str_p2"]
        circ_from_qasm=QuantumCircuit.from_qasm_file(os.path.join(BASE_PATH, rand_circ_files[idx]))
        circ_pieces=split_circuit_by_barrier(os.path.join(BASE_PATH, rand_circ_files[idx]))
        circ=circ_pieces[1]
        # If there were no Pauli checks found just skip.
        if max_pauli_weight == 0: 
            continue

        #Initial zero state. The randomness is introduced with the circuit.
        zero_state=Statevector.from_int(0, 2**(NUMBER_OF_QUBITS+1))

        #Create the complete circuit without a measurement
        quantum_register=QuantumRegister(NUMBER_OF_QUBITS)
        ancilla_register=QuantumRegister(1)
        all_qubits=[]
        for reg in quantum_register:
            all_qubits.append(reg)
        all_qubits.append(ancilla_register[0])
        circ_full=QuantumCircuit(quantum_register, ancilla_register)
        #Create the circuit with no checks
        circ_no_checks=QuantumCircuit(quantum_register, ancilla_register)
        add_rand_input_state(NUMBER_OF_QUBITS, quantum_register, circ_full, circ_no_checks)

        print(circ_no_checks)
        print(circ_full)

        circ_full.barrier()
        circ_full.compose(circ_from_qasm, inplace=True)
        circ_no_checks.barrier()
        circ_no_checks.compose(circ_pieces[1], inplace=True)

        #Test
        # print(list(range(NUMBER_OF_QUBITS)))
        comp_qubit_positions=list(range(NUMBER_OF_QUBITS))

        # Finish the full circuit
        # circ_full.h(ancilla_register[0])
        # circ_full.barrier(all_qubits)
        # circ_full=mymodule.add_controlU(circ_full, max_pauli_str_p1, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
        # circ_full.barrier(all_qubits)
        # circ_full.compose(circ, qubits=comp_qubit_positions, inplace=True)
        # circ_full.barrier(all_qubits)
        # circ_full=mymodule.add_controlU(circ_full, max_pauli_str_p2, NUMBER_OF_QUBITS, quantum_register, ancilla_register)
        # circ_full.barrier(all_qubits)
        # circ_full.h(ancilla_register[0])
        # print(circ_full.decompose())

        #Finish the circuit with no checks
        # circ_no_checks.compose(circ, qubits=comp_qubit_positions, inplace=True)
        # print(circ_no_checks.decompose())

        #Test
        print("circ")
        print(circ_from_qasm)
        print("circ full")
        print(circ_full)
        print("circ no checks")
        print(circ_no_checks)

        # The final state. Note that evolve doesn't change the object.
        final_state_no_errors_with_checks=zero_state.evolve(circ_full)
        final_state_no_errors_no_checks=zero_state.evolve(circ_no_checks)
        # print("No errors full state no checks: ", final_state_no_errors_no_checks)
        # print("No errors full state with checks: ", final_state_no_errors_checks)
        #Reduced final since we will be conducting tomography excluding the ancilla.
        reduced_correct_state=partial_trace(final_state_no_errors_no_checks, [NUMBER_OF_QUBITS])

        #Test.
        basis_gates=['id', 'rz', 'sx', 'cx']
        circ_full=transpile(circ_full, basis_gates=basis_gates, optimization_level=0)
        circ_no_checks=transpile(circ_no_checks, basis_gates=basis_gates, optimization_level=0)

        # Tomography stuff.
        # No checks tomography
        if do_tomography:
            circ_no_checks_tomo=state_tomography_circuits(circ_no_checks, quantum_register)
            # With Checks tomography with classical register
            circ_full_no_ancilla_tomo=state_tomography_circuits(circ_full, quantum_register)
        else:
            circ_no_checks_tomo=None
            circ_full_no_ancilla_tomo=None
        # Test
        # circ_no_checks_tomo=transpile(circ_no_checks_tomo, basis_gates=basis_gates, optimization_level=0)
        # circ_full_no_ancilla_tomo=transpile(circ_full_no_ancilla_tomo, basis_gates=basis_gates, optimization_level=0)
        # circ_full_with_ancilla_tomo=transpile(circ_full_with_ancilla_tomo, basis_gates=basis_gates, optimization_level=0)
        # print("transpiled circs")
        # print(circ_no_checks)
        # print(circ_no_checks)
        if do_tomography:
            circ_full_with_ancilla_tomo=deepcopy(circ_full_no_ancilla_tomo)
            classical_register=ClassicalRegister(1)
            for circ_temp in circ_full_with_ancilla_tomo:
                circ_temp.add_register(classical_register)
                circ_temp.measure(ancilla_register[0], classical_register)
            # print(circ_full.decompose())
        else:
            circ_full_with_ancilla_tomo=None
        print(circ.count_ops())

        # simulator=Aer.get_backend('qasm_simulator')
        # simulator.set_option("method", "density_matrix")
        # # # simulator.set_option("max_parallel_experiments", math.floor(psutil.cpu_count(logical=False)/cpus_pool))
        # simulator.set_option("max_parallel_experiments", psutil.cpu_count(logical=False))

        # Parallel code.
        if PARALLEL:
            # with Pool(psutil.cpu_count(logical=False), initialize, initargs=(NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files,
            #     one_qubit_error_space, file_info, file_number, circ, max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, 
            #     quantum_register, ancilla_register, all_qubits, circ_full, circ_no_checks, final_state_no_errors_with_checks,
            #     final_state_no_errors_no_checks, reduced_correct_state, file_name, circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo)) as pool:
            #     results=list(pool.imap_unordered(get_results, enumerate(one_qubit_error_space)))

            # For HPC we cap cpus_pool to 21 since that's the number of error indices. We can't do this for machines with less than 21 cpus.
            cpus_pool=psutil.cpu_count(logical=False)
            # with Pool(cpus_pool, initialize, initargs=(NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files,
            #     one_qubit_error_space, file_info, file_number, circ, max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, 
            #     quantum_register, ancilla_register, all_qubits, circ_full, circ_no_checks, final_state_no_errors_with_checks,
            #     final_state_no_errors_no_checks, reduced_correct_state, file_name, circ_no_checks_tomo, circ_full_no_ancilla_tomo, circ_full_with_ancilla_tomo, cpus_pool, simulator)) as pool:
            #     results=list(pool.imap_unordered(get_results, enumerate(one_qubit_error_space)))

            with Pool(cpus_pool, initialize, initargs=(NUMBER_OF_QUBITS, subdir, code_dir, rand_circ_files,
                one_qubit_error_space, file_info, file_number, circ, max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2, 
                quantum_register, ancilla_register, all_qubits, circ_full, circ_no_checks, final_state_no_errors_with_checks,
                final_state_no_errors_no_checks, reduced_correct_state, file_name, circ_no_checks_tomo, circ_full_no_ancilla_tomo, 
                circ_full_with_ancilla_tomo, cpus_pool, do_tomography)) as pool:
                if do_tomography:
                    results=list(pool.imap_unordered(get_results, enumerate(one_qubit_error_space)))
                else:
                    print("here")
                    results=list(pool.imap_unordered(get_results, enumerate(one_qubit_error_space)))
        else:
            # Non parallel.
            results=[]
            for error_idx, one_qubit_err in enumerate(one_qubit_error_space):
                if  error_idx==11:
                    results.append(get_results((error_idx, one_qubit_err)))
                #Test
                # if error_idx==1:
                    break
        #Store results.
        store_results(circ_full, BASE_PATH, file_name, NUMBER_OF_QUBITS, results, CNOT_COUNT, circ_info["rz"], circ_info["circuit_num"], 
            circ_info["found_matches"], max_pauli_weight, max_pauli_str_p1, max_pauli_str_p2)

        print("execution time", time.time()-time0)
    print("Finished.")

    # #Test
    # #Profiling stuff
    # pr.disable()
    # s = io.StringIO()
    # # sortby = SortKey.CUMULATIVE
    # # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # # ps.print_stats()
    # file_name="stats.txt"
    # with open(file_name, "w") as file:
    #     file.write(s.getvalue())    