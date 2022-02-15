#!/usr/bin/env python3

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import circtester 
import psutil
import time
from multiprocessing import Pool
import sys

if __name__ == "__main__":
    time0=time.time()
    print("running...")
    #Program parameters.
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    START_CIRC_NUMBER=int(sys.argv[3])
    END_CIRC_NUMBER=int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=True
    # File stuff
    SUBDIR=os.path.join("data", "results")
    # Gets the file path of the script
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    files_manipulator=circtester.FilesManipulator(BASE_PATH, NUMBER_OF_QUBITS, CNOT_COUNT)
    # Gets the files that match the string. Files include the path.
    # The returned files correspond accordingly, e.g., circ_file[0] and circ_properties_files[0] refer
    # to the same circuit.
    circ_files, circ_properties_files=files_manipulator.get_files(START_CIRC_NUMBER, END_CIRC_NUMBER)
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))

    #Get the qasm and pickle info
    for file_idx, file_name in enumerate(circ_files):
        #Stopwatch
        time1=time.time()
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        if files_manipulator.result_exists(file_name):
            continue

        circuit_maker=circtester.CircuitMaker(os.path.join(BASE_PATH, file_name), NUMBER_OF_QUBITS)
        noiseless_circs, noiseless_circs_meas=circuit_maker.make_noiseless_circs()

        keep_qubits=list(range(NUMBER_OF_QUBITS))
        circ_tester=circtester.CircuitSimulator(noiseless_circs, NUMBER_OF_QUBITS, keep_qubits)
        assert circ_tester.sanity_check_fidelity>0.98, f"Sanity check fidelity {circ_tester.sanity_check_fidelity} failed for circuit {file_name}"
        print(f"sanity check fidelity: {circ_tester.sanity_check_fidelity}")

        if PARALLEL:
            results=circ_tester.simulate_all_tests_parallel(pool, SINGLE_QUBIT_ERROR_SPACE)
        else:
            results=circ_tester.simulate_all_tests(SINGLE_QUBIT_ERROR_SPACE)
        files_manipulator.store_fidelity_results_verbose(circ_properties_files[file_idx], noiseless_circs, results)

        print(f"file execution time {time.time()-time1}")
    if PARALLEL:
        pool.close()
        pool.join()
    print(f"total execution time {time.time()-time0}")
    print("Finished.")