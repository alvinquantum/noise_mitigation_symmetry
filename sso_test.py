#!/usr/bin/env python3

import os
from cirq.contrib.qasm_import.qasm import circuit_from_qasm
from copy import deepcopy
from qiskit.circuit.quantumcircuit import QuantumCircuit, QuantumRegister
from qiskit.circuit.quantumregister import Qubit
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import circtester 
import psutil
import time
from multiprocessing import Pool
import sys
from qiskit import transpile

def make_circuits_from_result(number_of_qubits, circ_file):
    qiskit_circ_with_checks_store=QuantumCircuit.from_qasm_file(circ_file)
    qiskit_temp=deepcopy(qiskit_circ_with_checks_store.remove_final_measurements(inplace=False))
    circ_pieces=circtester.CircuitMaker.split_qasm_str_by_barrier(qiskit_temp.qasm())
    qubits_label="q"
    quantum_register=QuantumRegister(number_of_qubits+1, qubits_label)    
    qiskit_circ_no_checks=QuantumCircuit(quantum_register)
    qiskit_circ_with_checks=QuantumCircuit(quantum_register)
    circ_initial=circ_pieces[0]
    circ_compute=circ_pieces[3]
    # Rebuild qiskit circ
    for piece in circ_pieces:
        qiskit_circ_with_checks.compose(piece, inplace=True)
    qiskit_circ_no_checks.compose(circ_initial, inplace=True)
    qiskit_circ_no_checks.compose(circ_compute, inplace=True)
    basis_gates=['u1', 'u2', 'u3', 'cx']
    qiskit_circ_with_checks=transpile(qiskit_circ_with_checks, basis_gates=basis_gates, optimization_level=0)
    qiskit_circ_no_checks=transpile(qiskit_circ_no_checks, basis_gates=basis_gates, optimization_level=0)
    cirq_circ_with_checks=circuit_from_qasm(qiskit_circ_with_checks.qasm())
    cirq_circ_no_checks=circuit_from_qasm(qiskit_circ_no_checks.qasm())
    cirq_circ_measurements_with_checks=deepcopy(cirq_circ_with_checks)
    cirq_circ_measurements_no_checks=deepcopy(cirq_circ_no_checks)
    circtester.CircuitMaker.add_measurements(cirq_circ_measurements_with_checks)
    circtester.CircuitMaker.add_measurements(cirq_circ_measurements_no_checks)
    return (circtester.NoiselessCircuits(qiskit_circ_with_checks_store, cirq_circ_with_checks, cirq_circ_no_checks, qubits_label),
            circtester.NoiselessMeasurementCircuits(cirq_circ_measurements_with_checks, cirq_circ_measurements_no_checks, qubits_label))

if __name__ == "__main__":
    time0=time.time()
    print("running...")
    #Program parameters.
    NUMBER_OF_QUBITS=5#int(sys.argv[1])
    CNOT_COUNT=25#int(sys.argv[2])
    START_CIRC_NUMBER=0#int(sys.argv[3])
    END_CIRC_NUMBER=0#int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=False
    # File stuff
    SUBDIR="data"
    # Gets the file path of the script
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    files_manipulator=circtester.FilesManipulator(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT)
    # Gets the files that match the string. Files include the path.
    # The returned files correspond accordingly, e.g., circ_file[0] and circ_properties_files[0] refer
    # to the same circuit.
    file_name="qubits_5_CNOTS_40_circuit_45_result_0_.qasm"
    if files_manipulator.result_sso_exists(file_name):
        print(f"total execution time {time.time()-time0}")
        print("Finished.")
        quit()
        
    circ_file, circ_property_file=files_manipulator.get_result_qasm_and_original_prop(BASE_PATH, file_name)
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))

    #Stopwatch
    # circ_pieces=circtester.CircuitMaker.split_circuit_by_barrier()
    noiseless_circs, noiseless_circs_meas=make_circuits_from_result(NUMBER_OF_QUBITS, os.path.join(BASE_PATH, circ_file))    

    circ_tester=circtester.CircuitRunner(noiseless_circs_meas, NUMBER_OF_QUBITS)
    # assert circ_tester.sanity_check_sso>0.90, f"Sanity check sso {circ_tester.sanity_check_sso} failed for circuit {file_name}"
    print(f"sanity check sso: {circ_tester.sanity_check_sso}")

    if PARALLEL:
        results=circ_tester.run_all_tests_parallel(pool, SINGLE_QUBIT_ERROR_SPACE)
    else:
        results=circ_tester.run_all_tests(SINGLE_QUBIT_ERROR_SPACE)
    files_manipulator.store_sso_results(circ_property_file, noiseless_circs, results)

    if PARALLEL:
        pool.close()
        pool.join()
    print(f"total execution time {time.time()-time0}")
    print("Finished.")