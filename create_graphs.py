#!/usr/bin/env python3

import os, sys, glob, pickle, operator
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile

from qiskit.quantum_info.synthesis import one_qubit_decompose

def create_fidelity_plot_cnots(number_of_qubits,cnot_count):
    '''Creates plots based on CNOT counts.'''
    # global SUBDIR, BASE_PATH, ONE_QUBIT_ERROR_SPACE, NUMBER_OF_ERROR_POINTS
    # Gets the file path of the script.
    code_dir=sys.path[0]
    all_files=[f for f in listdir(BASE_PATH) if isfile(os.path.join(BASE_PATH, f))]

    target_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if ".obj" in name_split and "result" in name_split and name_split_nums[1]==cnot_count and name_split_nums[0]==number_of_qubits: #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            target_files.append(file)
    
    if target_files:
        # Load the results
        all_circs_results_sorted=[]
        for file_name in target_files:
            #We're looking at one circuit
            # Load data from each pickle and get the results object for the circuit. 
            with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
                    pickle_info=pickle.load(circ_file)
                    circ_results=pickle_info["results"]
                    assert pickle_info["qubits"]==number_of_qubits, file_name+" does not have the required qubit count."
                    assert pickle_info["cx"]==cnot_count, file_name+" does not have the required CNOT count."
            
            # Get the results in order. 
            circ_results_sorted=sorted(circ_results, key=lambda d: d['error_idx'])
            assert len(circ_results_sorted)==len(circ_results), "error when sorting outcomes."
            all_circs_results_sorted.append(circ_results_sorted)

        assert len(all_circs_results_sorted)==len(target_files), "lost some files from sorting."
        
        average_fidelity_gains=[]
        fidelity_gains=[]
        for circ_results in all_circs_results_sorted:
            circ_outcomes_temp=[]
            # print(len(circ))
            for idx, rnd in enumerate(circ_results):
                # print(idx)
                circ_outcomes_temp.append(rnd["state_fidelity_with_checks_with_errors"]-rnd["state_fidelity_no_checks_with_errors"])
            fidelity_gains.append(circ_outcomes_temp)
        
        average_fidelity_gains=[sum(x)/len(x) for x in zip(*fidelity_gains)]
        print("cnot count: ", cnot_count)
        print(average_fidelity_gains)
        std_deviations=[np.std(x) for x in zip(*fidelity_gains)]
        print(std_deviations)

        plt.clf()
        plt.errorbar(ONE_QUBIT_ERROR_SPACE, average_fidelity_gains, yerr=std_deviations, capsize=3.0)
        plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Fidelity Gain")
        plt.xscale("log")
        plt.savefig(BASE_PATH+"Qubits_"+str(number_of_qubits)+ "_CNOTS_" + str(cnot_count)+".png")
        # plt.show()
        return (cnot_count, average_fidelity_gains, std_deviations)

def create_all_fidelity_plot_cnots(number_of_qubits, cnot_counts):
    '''Wrapper for create_fidelity_plot_cnots(). Makes individual plots with all the lines and a single plot'''
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_fidelity_plot_cnots(number_of_qubits, cnot_count)
        if result:
            all_results.append(result)


    for elem in all_results:
        plt.errorbar(ONE_QUBIT_ERROR_SPACE, elem[1], yerr=elem[2], capsize=3.0, label="CNOT Count: "+str(elem[0]))
    plt.title("Qubits: "+str(number_of_qubits)+ " All CNOTS")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Fidelity Gain")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(BASE_PATH+"Qubits_"+str(number_of_qubits)+ "_allCNOTS_.png")
    # plt.show()

def create_count_percentage_plot(number_of_qubits, cnot_count):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    code_dir=sys.path[0]
    all_files=[f for f in listdir(BASE_PATH) if isfile(os.path.join(BASE_PATH, f))]

    target_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if ".obj" in name_split and "result" in name_split and name_split_nums[1]==cnot_count and name_split_nums[0]==number_of_qubits: #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            target_files.append(file)
    
    # Load the results
    all_circs_results_sorted=[]
    for file_name in target_files:
        #We're looking at one circuit
        # Load data from each pickle and get the results object for the circuit. 
        with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
                pickle_info=pickle.load(circ_file)
                circ_results=pickle_info["results"]
                assert pickle_info["qubits"]==number_of_qubits, file_name+" does not have the required qubit count."
                assert pickle_info["cx"]==cnot_count, file_name+" does not have the required CNOT count."
        
        # Get the results in order. 
        circ_results_sorted=sorted(circ_results, key=lambda d: d['error_idx'])
        assert len(circ_results_sorted)==len(circ_results), "error when sorting outcomes."
        all_circs_results_sorted.append(circ_results_sorted)

    assert len(all_circs_results_sorted)==len(target_files), "lost some files from sorting."
    
    count_diffs=[]
    for circ_results in all_circs_results_sorted:
        circ_outcomes_temp=[]
        # print(len(circ))
        for idx, rnd in enumerate(circ_results):
            # print(idx)
            circ_outcomes_temp.append(rnd["percent_results_after_postselect"]/rnd["percent_results_before_postselect"])
        count_diffs.append(circ_outcomes_temp)
    
    average_counts_diffs=[sum(x)/len(x) for x in zip(*count_diffs)]
    print(average_counts_diffs)
    std_deviations=[np.std(x) for x in zip(*count_diffs)]
    print(std_deviations)

    plt.clf()
    plt.errorbar(ONE_QUBIT_ERROR_SPACE, average_counts_diffs, yerr=std_deviations, capsize=3.0)
    plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Probability of Zero Measurement on Ancilla")
    plt.xscale("log")
    plt.savefig(f"{BASE_PATH}Qubits_{number_of_qubits}_CNOTS_{cnot_count}_.png")
    # plt.show()
    return (cnot_count, average_counts_diffs, std_deviations)

def create_all_count_percentage_plots(number_of_qubits, cnot_counts):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_count_percentage_plot(number_of_qubits, cnot_count)
        if result:
            all_results.append(result)


    for elem in all_results:
        plt.errorbar(ONE_QUBIT_ERROR_SPACE, elem[1], yerr=elem[2], capsize=3.0, label="CNOT Count: "+str(elem[0]))
    plt.title("Qubits: "+str(number_of_qubits)+ " All CNOTS")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Probability of Zero Measurement on Ancilla")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(BASE_PATH+"Qubits_"+str(number_of_qubits)+ "_allCounts_.png")
    

if __name__ == "__main__":
    # File path stuff
    CODE_DIR=sys.path[0]
    # SUBDIR="/qubits_5_results/"
    SUBDIR="/data/"

    BASE_PATH=CODE_DIR+SUBDIR
    os.chdir(BASE_PATH)
    NUMBER_OF_QUBITS=10
    # CNOT_COUNT=30
    NUMBER_OF_ERROR_POINTS=21
    # NUMBER_OF_ERROR_POINTS=23
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
#     ONE_QUBIT_ERROR_SPACE=[1.00000000e-05, 1.41253754e-05, 1.99526231e-05, 2.81838293e-05,
#  3.98107171e-05, 5.62341325e-05, 7.94328235e-05, 1.12201845e-04,
#  1.58489319e-04, 2.23872114e-04, 3.16227766e-04, 4.46683592e-04,
#  6.30957344e-04, 8.91250938e-04, 1.25892541e-03, 1.77827941e-03,
#  2.51188643e-03, 3.54813389e-03, 5.01187234e-03, 7.07945784e-03,
#  1.00000000e-02, 0.015, 0.02, 0.025, 0.03, 0.035, 0.06]
    cnot_counts=[1,5,10,15,20,25,30,35,40]
    create_all_fidelity_plot_cnots(NUMBER_OF_QUBITS, cnot_counts)
    # Plots the percentages of circuits we keep.
    create_all_count_percentage_plots(NUMBER_OF_QUBITS, cnot_counts)