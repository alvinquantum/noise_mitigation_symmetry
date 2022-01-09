#!/usr/bin/env python3

import os 
from copy import deepcopy
import numpy as np
import circtester 
import psutil
import pickle
import time
from multiprocessing import Pool, Value, Array, Manager
from cirq.contrib.qasm_import import circuit_from_qasm

if __name__ == "__main__":
    time0=time.time()
    print("running...")
    #Program parameters.
    NUMBER_OF_QUBITS=5#int(sys.argv[1])
    CNOT_COUNT=10#int(sys.argv[2])
    START_CIRC_NUMBER=0#int(sys.argv[3])
    END_CIRC_NUMBER=10#int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=False
    # File stuff
    SUBDIR="/data/"
    # Gets the file path of the script
    # code_dir=sys.path[0]
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=CODE_DIR+SUBDIR
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    # SINGLE_QUBIT_ERROR_SPACE=[
    #     1.00000000e-05, 1.41253754e-05, 1.99526231e-05, 2.81838293e-05,
    #     3.98107171e-05, 5.62341325e-05, 7.94328235e-05, 1.12201845e-04,
    #     1.58489319e-04, 2.23872114e-04, 3.16227766e-04, 4.46683592e-04,
    #     6.30957344e-04, 8.91250938e-04, 1.25892541e-03, 1.77827941e-03,
    #     2.51188643e-03, 3.54813389e-03, 5.01187234e-03, 7.07945784e-03,
    #     1.00000000e-02, 0.015, 0.02, 0.025, 0.03, 0.035, 0.06]
    # Gets the files that match the string. Files include the path
    files_manipulator=circtester.FilesManipulator(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT)
    circ_files, circ_properties_files=files_manipulator.get_files(START_CIRC_NUMBER, END_CIRC_NUMBER)

    #Get the qasm and pickle info
    for file_idx, file_name in enumerate(circ_files):
        #Stopwatch
        time1=time.time()
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        if files_manipulator.result_exists(file_name):
            continue

        circuit_maker=circtester.CircuitMaker(os.path.join(BASE_PATH, file_name), NUMBER_OF_QUBITS)
        noiseless_circs=circuit_maker.make_noiseless_circs()

        keep_qubits=list(range(NUMBER_OF_QUBITS))
        circ_tester=circtester.CircuitTester(noiseless_circs, NUMBER_OF_QUBITS, keep_qubits)
        assert circ_tester.sanity_check_fidelity, f"Sanity check fidelity failed for circuit {file_name}"

        results=[]
        for error_idx, single_qubit_error in enumerate(SINGLE_QUBIT_ERROR_SPACE):
            results.append(circ_tester.run_test((error_idx, single_qubit_error)))

        files_manipulator.store_results(circ_properties_files[file_idx], noiseless_circs, results)

        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    print("Finished.")