#!/usr/bin/env python3
import os
import graphplotter
import time
import checksfinder
import circgenerator
import circtester

def run_sims():
    '''running everything from scratch. requires the rigid
    folder structure described in the readme.'''
    # Qubits 5 rz count 5 low weight checks of 6 layers.
    num_circs=50 #number of circs for each cnot count.
    start_circ=0
    end_circ=num_circs-1
    cnots=[1,5,10,15,20,25,30,35,40]
    qubits=5
    target_layers=6#*qubits
    max_layers=target_layers+target_layers/3 #some weight one checks will be redundant.
    rz_count=5
    for cnot_count in cnots:
        circgenerator.gen_rand_circs(qubits, cnot_count, num_circs, rz_count)
        circgenerator.gen_rand_input_state(qubits, cnot_count, rz_count)
    for cnot_count in cnots: 
        #More efficient to separate the simulations this since it uses mulitprocessing
        #accross the all cnot counts that are equal.
        checksfinder.create_low_weight_checks_circs(qubits, cnot_count, target_layers, max_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "low_weight_checks")
    rz_count=10
    for cnot_count in cnots:
        circgenerator.gen_rand_circs(qubits, cnot_count, num_circs, rz_count)
        circgenerator.gen_rand_input_state(qubits, cnot_count, rz_count)
    for cnot_count in cnots:
        checksfinder.create_low_weight_checks_circs(qubits, cnot_count, target_layers, max_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "low_weight_checks")
    rz_count=15
    for cnot_count in cnots:
        circgenerator.gen_rand_circs(qubits, cnot_count, num_circs, rz_count)
        circgenerator.gen_rand_input_state(qubits, cnot_count, rz_count)
    for cnot_count in cnots:
        checksfinder.create_low_weight_checks_circs(qubits, cnot_count, target_layers, max_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "low_weight_checks")
    #Only one layer for high weight
    target_layers=1
    max_layers=1
    rz_count=5
    for cnot_count in cnots:
        checksfinder.create_high_weight_checks_circs(qubits, cnot_count, target_layers, max_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "high_weight_checks")

    # Qubits 10
    num_circs=50 #number of circs for each cnot count.
    qubits=10
    cnots=[1,5,10,15,20,25,30,35,40, 80]
    rz_count=5
    # We only simulate one layer. It's too slow to do more.
    target_layers=1
    max_layers=1
    for cnot_count in cnots:
        circgenerator.gen_rand_circs(qubits, cnot_count, num_circs, rz_count)
        circgenerator.gen_rand_input_state(qubits, cnot_count, rz_count)
    for cnot_count in cnots: 
        #More efficient to separate the simulations this since it uses mulitprocessing
        #accross the all cnot counts that are equal.
        checksfinder.create_low_weight_checks_circs(qubits, cnot_count, target_layers, max_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "low_weight_checks")
    for cnot_count in cnots:
        checksfinder.create_high_weight_checks_circs(qubits, cnot_count, target_layers, max_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "high_weight_checks")

    #Qubits 2
    num_circs=50
    start_circ=0
    end_circ=num_circs-1
    cnots=[1,2,4,8,16,32,64,128,256,512,1024]
    qubits=2
    target_layers=4
    max_layers=4
    rz_count=0
    for cnot_count in cnots:
        circgenerator.gen_rand_circs(qubits, cnot_count, num_circs, rz_count)
        circgenerator.gen_rand_input_state(qubits, cnot_count, rz_count)
    for cnot_count in cnots: 
        #More efficient to separate the simulations this since it uses mulitprocessing
        #accross the all cnot counts that are equal.
        checksfinder.create_weight_one_check_circs(qubits, cnot_count, target_layers, rz_count)
        circtester.exec_simulations(qubits, rz_count, cnot_count, start_circ, end_circ, target_layers, "low_weight_checks")

    # Qubits 2 example unit fidelity for post selected state (noiseless checks.)
    num_circs=1
    start_circ=0
    end_circ=num_circs-1
    cnots=[30]
    qubits=2
    target_layers=4
    max_layers=4
    rz_count=0
    for cnot_count in cnots:
        circgenerator.gen_rand_circs_example(qubits, cnot_count, num_circs, rz_count)
        circgenerator.gen_rand_input_state_example(qubits, cnot_count, rz_count)
    for cnot_count in cnots: 
        # More efficient to separate the simulations this since it uses mulitprocessing
        # accross the all cnot counts that are equal.
        checksfinder.create_weight_one_check_circs_example(qubits, cnot_count, target_layers, rz_count)
        circtester.exec_simulations_noiseless_checks(qubits, rz_count, cnot_count, start_circ, end_circ, 
            target_layers, "low_weight_checks")


# For automating the process on a local machine. Don't use for HPC.
if __name__ == "__main__":
    time0=time.time()
    print("running...")
    # run_sims()
    
    # 5 compute qubits graphs.
    rz_count=5
    graphplotter.get_5qubits_fid_and_counts_vs_error_low_weightchks(rz_count, [6])
    graphplotter.get_5qubits_fid_vs_error_max_weightchks(rz_count)
    graphplotter.get_5qubits_fidelity_vs_layers(rz_count)
    graphplotter.get_5qubits_raw_fidel_and_counts_vs_error(rz_count)
    rz_count=10
    graphplotter.get_5qubits_fid_and_counts_vs_error_low_weightchks(rz_count, [6])
    rz_count=15
    graphplotter.get_5qubits_fid_and_counts_vs_error_low_weightchks(rz_count, [6])
    graphplotter.get_5qubits_fid_and_counts_vs_error_low_weightchks(rz_count, [1])

    # 2 compute qubits graphs.
    graphplotter.get_2qubits_unit_fidel_ex()
    graphplotter.get_2qubits_fidel_counts_vs_cnots()

    # # 10 compute qubits graphs.
    rz_count=5
    graphplotter.get_10qubits_fid_and_counts_vs_error_low_weightchks(rz_count, [1])
    graphplotter.get_10qubits_fid_vs_error_max_weightchks(rz_count)
    print(f"total execution time {time.time()-time0}")

    print("done.")