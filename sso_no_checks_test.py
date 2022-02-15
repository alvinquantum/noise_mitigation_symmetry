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
from qiskit import transpile

if __name__ == "__main__":
    time0=time.time()
    print("running...")
    #Program parameters.
    NUMBER_OF_QUBITS=5#int(sys.argv[1])
    CNOT_COUNT=40#int(sys.argv[2])
    START_CIRC_NUMBER=0#int(sys.argv[3])
    END_CIRC_NUMBER=0#int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=True
    # File stuff
    SUBDIR="sso_checks_and_no_checks"
    # Gets the file path of the script
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    files_manipulator=circtester.FilesManipulator(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT)
    # Gets the files that match the string. Files include the path.
    # The returned files correspond accordingly, e.g., circ_file[0] and circ_properties_files[0] refer
    # to the same circuit.
  
    qasm_files=circtester.FilesManipulator.get_qasm_files(BASE_PATH)
    # file_name="qubits_5_CNOTS_40_circuit_45_result_0_nochecks_.qasm"
    for file_name in qasm_files:
        circtester.FilesManipulator.remove_measurements_in_qasm_file(BASE_PATH, file_name)
    
    # Strip the measurements
    # file_name_nochecks="qubits_5_CNOTS_25_circuit_0_result_0_nochecks_.qasm"
    # file_name_withchecks="qubits_5_CNOTS_25_circuit_0_result_0_checks_.qasm"
    file_name_nochecks="qubits_5_CNOTS_40_circuit_45_result_0_nochecks_.qasm"
    file_name_withchecks="qubits_5_CNOTS_40_circuit_45_result_0_checks_.qasm"

    circtester.FilesManipulator.remove_checks_in_complete_qasm(BASE_PATH, file_name_withchecks, file_name_nochecks)

    #Transpile
    basis_gates=['u1', 'u2', 'u3', 'cx']
    qiskit_circ_with_checks=transpile(QuantumCircuit.from_qasm_file(os.path.join(BASE_PATH,file_name_withchecks)), optimization_level=0)
    qiskit_circ_no_checks=transpile(QuantumCircuit.from_qasm_file(os.path.join(BASE_PATH,file_name_nochecks)), optimization_level=0)

    cirq_circ_measurements_no_checks=circtester.CircuitMaker.qasm_string_to_cirq_circ(qiskit_circ_no_checks.qasm())
    cirq_circ_measurements_with_checks=circtester.CircuitMaker.qasm_string_to_cirq_circ(qiskit_circ_with_checks.qasm())
    circtester.CircuitMaker.add_measurements(cirq_circ_measurements_no_checks)
    circtester.CircuitMaker.add_measurements(cirq_circ_measurements_with_checks)
    
    noiseless_circs_meas=circtester.NoiselessMeasurementCircuits(cirq_circ_measurements_with_checks, cirq_circ_measurements_no_checks, "q")
    if files_manipulator.result_sso_exists(file_name):
        print(f"total execution time {time.time()-time0}")
        print("Finished.")
        quit()
        
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))    

    circ_tester=circtester.CircuitRunner(noiseless_circs_meas, NUMBER_OF_QUBITS, 4000)
    # # assert circ_tester.sanity_check_sso>0.90, f"Sanity check sso {circ_tester.sanity_check_sso} failed for circuit {file_name}"
    print(f"sanity check sso: {circ_tester.sanity_check_sso}")

    if PARALLEL:
        results=circ_tester.run_all_tests_parallel(pool, SINGLE_QUBIT_ERROR_SPACE)
    else:
        results=circ_tester.run_all_tests(SINGLE_QUBIT_ERROR_SPACE)
    files_manipulator.store_basic_sso_result(BASE_PATH, file_name_withchecks, results)

    if PARALLEL:
        pool.close()
        pool.join()

    # Add measurements to qasm files
    qasm_files=circtester.FilesManipulator.get_qasm_files(BASE_PATH)
    for file_name in qasm_files:
        circtester.FilesManipulator.add_all_measurements_in_qasm_file(BASE_PATH, file_name)
    print(f"total execution time {time.time()-time0}")
    print("Finished.")