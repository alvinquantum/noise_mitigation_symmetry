#!/usr/bin/env python3

import os, sys, pickle
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile

def get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits, result_str="result"):
    '''Gets the sorted (by error idx, i.e., smallest single qubit error to largest) test results for 
    the given cnot count and number of qubits.'''
    all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
    # print(all_files)
    target_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names for checking
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        if (
            ".obj" in name_split and result_str in name_split 
            and name_split_nums[1]==cnot_count 
            and name_split_nums[0]==number_of_qubits): #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            target_files.append(file)

    all_circ_results_sorted=[]
    if target_files:
        # Load the results
        for file_name in target_files:
            #We're looking at one circuit
            # Load data from each pickle and get the results object for the circuit. 
            with open(os.path.join(base_path, file_name), "rb") as circ_file:
                    pickle_info=pickle.load(circ_file)
                    circ_result=pickle_info["results"]
                    assert pickle_info["qubits"]==number_of_qubits, f"{file_name} does not have the required qubit count."
                    assert pickle_info["cx"]==cnot_count, f"{file_name} does not have the required CNOT count."
            
            # Get the results in order. 
            circ_result_sorted=sorted(circ_result, key=lambda d: d['error_idx'])
            assert len(circ_result_sorted)==len(circ_result), "error when sorting outcomes."
            all_circ_results_sorted.append(circ_result_sorted)

        assert len(all_circ_results_sorted)==len(target_files), "lost some files from sorting."
    return all_circ_results_sorted

def create_fidelity_plot_cnot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits)
    if all_sorted_circ_results:    
        average_fidelity_gains=[]
        fidelity_gains=[]
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]

            for result in circ_result:
                circ_outcomes_temp.append(result["state_fidelity_with_checks_with_errors"]-result["state_fidelity_no_checks_with_errors"])
            fidelity_gains.append(circ_outcomes_temp)
        
        average_fidelity_gains=[sum(x)/len(x) for x in zip(*fidelity_gains)]
        print(average_fidelity_gains)
        std_deviations=[np.std(x) for x in zip(*fidelity_gains)]
        print(std_deviations)

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_fidelity_gains, yerr=std_deviations, capsize=3.0)
        plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Fidelity Gain")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_fidelity_.png"))
        return (cnot_count, average_fidelity_gains, std_deviations)

def create_all_fidelity_plot_cnots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space):
    '''Wrapper for create_fidelity_plot_cnots(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_fidelity_plot_cnot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)


    for elem in all_results:
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, label="CNOT Count: "+str(elem[0]))
    plt.title("Qubits: "+str(number_of_qubits)+ " All CNOTS")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Fidelity Gain")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_fidelity_.png"))

def create_count_percentage_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits)

    if all_sorted_circ_results:
        count_diffs=[]
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]
            for result in circ_result:
                # In our simulations we have the full density matrix so the percent results before post selecting is actually 1.
                # We do it this way so the graphing can accomidate physical runs where these two values are actually counts.
                circ_outcomes_temp.append(result["percent_results_after_postselect"]/result["percent_results_before_postselect"])
            count_diffs.append(circ_outcomes_temp)
        
        average_counts_diffs=[sum(x)/len(x) for x in zip(*count_diffs)]
        print(average_counts_diffs)
        std_deviations=[np.std(x) for x in zip(*count_diffs)]
        print(std_deviations)

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_counts_diffs, yerr=std_deviations, capsize=3.0)
        plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Probability of Zero Measurement on Ancilla")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_counts_.png"))
        return (cnot_count, average_counts_diffs, std_deviations)

def create_all_count_percentage_plots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_count_percentage_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)

    for elem in all_results:
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, label="CNOT Count: "+str(elem[0]))
    plt.title("Qubits: "+str(number_of_qubits)+ " All CNOTS")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Probability of Zero Measurement on Ancilla")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_counts_.png"))

def create_sso_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits, result_str="resultsso")
    if all_sorted_circ_results:    
        average_fidelity_gains=[]
        fidelity_gains=[]
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]

            for result in circ_result:
                circ_outcomes_temp.append(result["sso_with_checks_with_errors"]-result["sso_no_checks_with_errors"])
            fidelity_gains.append(circ_outcomes_temp)
        
        average_fidelity_gains=[sum(x)/len(x) for x in zip(*fidelity_gains)]
        print(average_fidelity_gains)
        std_deviations=[np.std(x) for x in zip(*fidelity_gains)]
        print(std_deviations)

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_fidelity_gains, yerr=std_deviations, capsize=3.0)
        plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average SSO Gain")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_sso_.png"))
        return (cnot_count, average_fidelity_gains, std_deviations)

def create_all_sso_plots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_sso_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)

    for elem in all_results:
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, label="CNOT Count: "+str(elem[0]))
    plt.title("Qubits: "+str(number_of_qubits)+ " All CNOTS")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average SSO Gain")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_sso_.png"))     

if __name__ == "__main__":
    # File path stuff
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=5
    # SUBDIR=f"qubits_{number_of_qubits}_results"
    SUBDIR="data"
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    os.chdir(BASE_PATH)
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
    cnot_counts=[1,5,10,15,20,25,30,35,40, 80]
    create_all_fidelity_plot_cnots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE)
    # Plots the percentages of circuits we keep.
    create_all_count_percentage_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE)
    create_all_sso_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE)
        