import checksfinder
import time
import os
import utilities
import circgenerator
import circtester
import sys


if __name__=="__main__":
    time0=time.time()
    #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    # NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    # func_params=FunctionParams()
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=f"qubits_{NUMBER_OF_QUBITS}"
    CHECKS_SUBDIR=os.path.join(MAIN_SUBDIR, "min_weight_checks", "checks")
    RESULTS_SUBDIR=os.path.join(MAIN_SUBDIR, "min_weight_checks", "results")
    INITIAL_STATES_SUBDIR=os.path.join(MAIN_SUBDIR, "initial_states")
    RAWS_SUBDIR=os.path.join(MAIN_SUBDIR, "raws")
    CHECKS_PATH=os.path.join(CODE_DIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, RAWS_SUBDIR)
    RESULTS_PATH=os.path.join(CODE_DIR, RESULTS_SUBDIR)
    INITIAL_STATES_PATH=os.path.join(CODE_DIR, INITIAL_STATES_SUBDIR)

    # PARALLEL=False
    # if PARALLEL:
    #     pool=Pool(psutil.cpu_count(logical=False))
    #Load qasm files and find checks.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=circgenerator.filter_by_cnot_qubit(input_qasm_file_names, CNOT_COUNT, NUMBER_OF_QUBITS)
    for input_qasm_file_name in input_qasm_file_names:
        time1=time.time()

        # print(qasm_file_name)
        # print(qasm_file_name)
        end_string="_raw_.qasm"
        output_file_name=f"{input_qasm_file_name[0:-len(end_string)]}_inputstate_0_.qasm"
        print(output_file_name)
        # Skip the raw qasm if we already have an output.
        if checksfinder.init_circ_exists(CHECKS_PATH, input_qasm_file_name):
            continue
        circ=circtester.CircuitMaker.make_rand_input_state_multilayer(NUMBER_OF_QUBITS)
        # print(circ)
        # print(INITIAL_STATES_PATH)
        # print(os.path.join(INITIAL_STATES_PATH, output_file_name))
        circ.qasm(filename= os.path.join(INITIAL_STATES_PATH, output_file_name))

        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    # if PARALLEL:
    #     pool.close() #Stops passing jobs to processes.
    #     pool.join() #Waits for processes to finish.
    print("done")