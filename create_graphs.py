#!/usr/bin/env python3

from cProfile import label
import os, sys, pickle
from tkinter import ON
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile
import circtester
import seaborn as sns

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
                    # assert pickle_info["qubits"]==number_of_qubits, f"{file_name} does not have the required qubit count."
                    # assert pickle_info["cx"]==cnot_count, f"{file_name} does not have the required CNOT count."
            
            # Get the results in order. 
            circ_result_sorted=sorted(circ_result, key=lambda d: d['error_idx'])
            assert len(circ_result_sorted)==len(circ_result), "error when sorting outcomes."
            all_circ_results_sorted.append(circ_result_sorted)

        assert len(all_circ_results_sorted)==len(target_files), "lost some files from sorting."
    return all_circ_results_sorted

def get_all_sorted_circ_results_by_layer(base_path, cnot_count, number_of_qubits, layers_count, result_str="result"):
    '''Gets the sorted (by error idx, i.e., smallest single qubit error to largest) test results for 
    the given cnot count and number of qubits.'''
    all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
    # print(all_files)
    target_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names for checking
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        if (f"layers_{layers_count}_" in file and
            ".obj" in name_split and result_str in name_split 
            and name_split_nums[1]==cnot_count 
            and name_split_nums[0]==number_of_qubits): #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            target_files.append(file)

    all_circ_results_sorted=[]
    if target_files:
        print(target_files)
        # Load the results
        for file_name in target_files:
            #We're looking at one circuit
            # Load data from each pickle and get the results object for the circuit. 
            with open(os.path.join(base_path, file_name), "rb") as circ_file:
                    pickle_info=pickle.load(circ_file)
                    circ_result=pickle_info["results"]
                    # assert pickle_info["qubits"]==number_of_qubits, f"{file_name} does not have the required qubit count."
                    # assert pickle_info["cx"]==cnot_count, f"{file_name} does not have the required CNOT count."
            
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
        
        print("raw fidelity", fidelity_gains)
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

def create_fidelity_plot_cnot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count):
    '''Creates plots based on CNOT counts and layers.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results_by_layer(base_path, cnot_count, number_of_qubits, layers_count)
    if all_sorted_circ_results:
        average_fidelity_gains=[]
        fidelity_gains=[]
        # Looking at one circuit.
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]
            # Looking at one error idx for this specific circuit.
            for result in circ_result:
                circ_outcomes_temp.append(result["state_fidelity_with_checks_with_errors"]-result["state_fidelity_no_checks_with_errors"])
            fidelity_gains.append(circ_outcomes_temp)
        
        # print("raw fidelity", fidelity_gains)
        average_fidelity_gains=[sum(x)/len(x) for x in zip(*fidelity_gains)]
        # print("average gains ", average_fidelity_gains)
        print(len(average_fidelity_gains))
        std_deviations=[np.std(x) for x in zip(*fidelity_gains)]
        print(len(std_deviations))
        print(len(one_qubit_error_space))
        # print(std_deviations)

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_fidelity_gains, yerr=std_deviations, capsize=3.0)
        plt.title(f"Compute Qubits: {str(number_of_qubits)}, CNOTS: {cnot_count}, Layers: {layers_count}")
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Fidelity Gain")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_fidelity_layers_{layers_count}_.png"))
        return (cnot_count, average_fidelity_gains, std_deviations)

def create_all_fidelity_plot_cnots_layers(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, layers_count):
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_fidelity_plot_cnot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count)
        if result:
            all_results.append(result)


    for idx, elem in enumerate(all_results):
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    plt.title(f"Compute Qubits: {str(number_of_qubits)}, Layers: {layers_count}")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Fidelity Gain")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_fidelity_gain_layers_{layers_count}_.png"))

def create_all_fidelity_plot_cnots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, NLINES):
    '''Wrapper for create_fidelity_plot_cnots(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_fidelity_plot_cnot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)


    for idx, elem in enumerate(all_results):
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    plt.title(f"Qubits: {number_of_qubits}, Layers: 1")
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
        plt.title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}")
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Probability of Zero Measurement on Ancilla")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_counts_.png"))
        return (cnot_count, average_counts_diffs, std_deviations)

def create_count_percentage_plot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results_by_layer(base_path, cnot_count, number_of_qubits, layers_count = layers_count)

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
        plt.title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}, Layers: {layers_count}")
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Probability of Zero Measurement on Ancilla")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_counts_layers_{layers_count}_.png"))
        return (cnot_count, average_counts_diffs, std_deviations)


def create_all_count_percentage_plots_layers(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, layers_count):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_count_percentage_plot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count=layers_count)
        if result:
            all_results.append(result)

    for idx, elem in enumerate(all_results):
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    plt.title(f"Compute Qubits: {str(number_of_qubits)}, Layers: {layers_count}")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Probability of Zero Measurement on Ancilla")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_counts_layers_{layers_count}_.png"))

def create_all_count_percentage_plots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, NLINES):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_count_percentage_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)

    for idx, elem in enumerate(all_results):
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    plt.title(f"Compute Qubits: {number_of_qubits}, Layers: 1")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Probability of Zero Measurement on Ancilla")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_counts_.png"))

def create_sso_gain_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
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

def create_all_sso_gain_plots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_sso_gain_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
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

def create_raw_sso_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits, result_str="resultsso")
    if all_sorted_circ_results:    
        # average_fidelity_gains=[]
        ssos_with_checks=[]
        ssos_no_checks=[]
        for circ_result in all_sorted_circ_results:
            temp_sso_with_checks=[]
            temp_sso_no_checks=[]
            for result in circ_result:
                temp_sso_with_checks.append(result["sso_with_checks_with_errors"])
                temp_sso_no_checks.append(result["sso_no_checks_with_errors"])
            ssos_with_checks.append(temp_sso_with_checks)
            ssos_no_checks.append(temp_sso_no_checks)

        print(ssos_with_checks)
        average_ssos_with_checks=[sum(x)/len(x) for x in zip(*ssos_with_checks)]
        average_ssos_no_checks=[sum(x)/len(x) for x in zip(*ssos_no_checks)]

        std_deviations_sso_with_checks=[np.std(x) for x in zip(*ssos_with_checks)]
        std_deviations_sso_no_checks=[np.std(x) for x in zip(*ssos_no_checks)]

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_ssos_with_checks, 
            yerr=std_deviations_sso_with_checks, capsize=3.0, label="SSO with Checks")
        plt.errorbar(one_qubit_error_space, average_ssos_no_checks, 
            yerr=std_deviations_sso_no_checks, capsize=3.0, label="SSO No Checks")
        plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        plt.xlabel("Single Qubit Error")
        plt.ylabel("SSO")
        plt.xscale("log")
        plt.legend()
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_rawsso_.png"))
        # return (cnot_count, average_fidelity_gains, std_deviations)

def create_split_layers_fidelity_plot_cnots(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    plt.clf()
    NLINES = len(layers)+1 #plus one for the zero checks.
    # NLINES = 45
    palette = sns.color_palette("viridis",NLINES)
    for layer_idx, layer in enumerate(layers):
        all_sorted_circ_results=get_all_sorted_circ_results_by_layer(base_path, cnot_count, number_of_qubits, layer)
        # print(all_sorted_circ_results)
        if all_sorted_circ_results:
            #Regular plots
            # average_fidelity_gains=[]
            # fidelity_gains=[]

            # Each element contains results for one error value
            fidelity_no_checks=[]
            fidelity_with_checks=[]
            # Looking at one circuit.
            for circ_result in all_sorted_circ_results:
                # circ_outcomes_temp=[]

                # Looking at one error value for that circuit.
                temp_fidelity_no_checks=[]
                temp_fidelity_with_checks=[]
                for result in circ_result:
                    temp_fidelity_no_checks.append(result["state_fidelity_no_checks_with_errors"])
                    temp_fidelity_with_checks.append(result["state_fidelity_with_checks_with_errors"])
                #     circ_outcomes_temp.append(result["state_fidelity_with_checks_with_errors"]-result["state_fidelity_no_checks_with_errors"])
                # fidelity_gains.append(circ_outcomes_temp)
                fidelity_no_checks.append(temp_fidelity_no_checks)
                fidelity_with_checks.append(temp_fidelity_with_checks)
            average_fidelity_with_checks=[sum(x)/len(x) for x in zip(*fidelity_with_checks)]
            # print("raw fidelity", fidelity_gains)
            # average_fidelity_with_checks=[sum(x)/len(x) for x in zip(*fidelity_with_checks)]
            # print(f"average fidelity with checks {average_fidelity_with_checks}")
            std_deviations_with_checks=[np.std(x) for x in zip(*fidelity_with_checks)]
            # print(f"std deviation with checks {std_deviations_with_checks}")

            average_fidelity_no_checks=[sum(x)/len(x) for x in zip(*fidelity_no_checks)]
            # print(f"average fidelity no checks {average_fidelity_no_checks}")
            std_deviations_no_checks=[np.std(x) for x in zip(*fidelity_no_checks)]
            # print(f"std deviation no checks {std_deviations_no_checks}")
            if layer_idx==0:
                plt.errorbar(one_qubit_error_space, average_fidelity_no_checks, std_deviations_no_checks, capsize=3.0, color=palette[0], label="Number of layers: 0")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            plt.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")
    plt.title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Fidelity")
    plt.xscale("log")
    plt.legend()
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_fidelity_.png"))
    # return (cnot_count, average_fidelity_gains, std_deviations)

def create_split_layers_count_percentage_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    plt.clf()
    NLINES = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", NLINES)
    for layer_idx, layer in enumerate(layers):
        all_sorted_circ_results=get_all_sorted_circ_results_by_layer(base_path, cnot_count, number_of_qubits, layer)
        # print(all_sorted_circ_results)
        if all_sorted_circ_results:
            #Regular plots
            # NLINES = 45
            # each element is a circuit.
            counts_after=[]
            counts_before=[]
            for circ_result in all_sorted_circ_results:
                counts_after_temp=[]
                counts_before_temp=[]
                # each element here is one error idx for the circuit.
                for result in circ_result:
                    # In our simulations we have the full density matrix so the percent results before post selecting is actually 1.
                    # We do it this way so the graphing can accomidate physical runs where these two values are actually counts.
                    counts_after_temp.append(result["percent_results_after_postselect"])
                    counts_before_temp.append(result["percent_results_before_postselect"])
                counts_after.append(counts_after_temp)
                counts_before.append(counts_before_temp)
            
            average_counts_after=[sum(x)/len(x) for x in zip(*counts_after)]
            # print(average_counts_diffs)
            std_deviations_after=[np.std(x) for x in zip(*counts_after)]
            # print(std_deviations)

            if layer_idx==0:
                plt.errorbar(one_qubit_error_space, counts_before[0], color=palette[0], capsize=3.0, label="Number of layers: 0")
            # for layer_idx, circ_count_after in enumerate(counts_after):
            #     layer_count=layer_idx+1
            plt.errorbar(one_qubit_error_space, average_counts_after, std_deviations_after,
                capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")
        plt.title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}")
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Probability of Zero Measurement on All Ancillas")
        plt.xscale("log")
        plt.legend()
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_counts_.png"))

def create_baseline_sso_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits, result_str="resultsso")
    if all_sorted_circ_results:    
        # average_fidelity_gains=[]
        ssos_with_checks=[]
        ssos_no_checks=[]
        for circ_result in all_sorted_circ_results:
            temp_sso_with_checks=[]
            temp_sso_no_checks=[]
            for result in circ_result:
                temp_sso_with_checks.append(result["correct_distribution"])
                temp_sso_no_checks.append(result["correct_distribution"])
                # print(result["correct_distribution"])
                sso = 0
                for v in result["correct_distribution"].values():
                    sso += np.sqrt(v / 2**5)
                sso = sso ** 2
                print(f"SSO: {sso}")
        #     ssos_with_checks.append(temp_sso_with_checks)
        #     ssos_no_checks.append(temp_sso_no_checks)

        # print(ssos_with_checks)
        # average_ssos_with_checks=[sum(x)/len(x) for x in zip(*ssos_with_checks)]
        # average_ssos_no_checks=[sum(x)/len(x) for x in zip(*ssos_no_checks)]

        # std_deviations_sso_with_checks=[np.std(x) for x in zip(*ssos_with_checks)]
        # std_deviations_sso_no_checks=[np.std(x) for x in zip(*ssos_no_checks)]

        # plt.clf()
        # plt.errorbar(one_qubit_error_space, average_ssos_with_checks, 
        #     yerr=std_deviations_sso_with_checks, capsize=3.0, label="SSO with Checks")
        # plt.errorbar(one_qubit_error_space, average_ssos_no_checks, 
        #     yerr=std_deviations_sso_no_checks, capsize=3.0, label="SSO No Checks")
        # plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        # plt.xlabel("Single Qubit Error")
        # plt.ylabel("SSO")
        # plt.xscale("log")
        # plt.legend()
        # plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_rawsso_.png"))

if __name__ == "__main__":
    # File path stuff
    CODE_DIR=sys.path[0]
    SUBDIR=os.path.join("qubits_5", "min_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    os.chdir(BASE_PATH)
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    
    # ax.plot(
    #             x,y,
    #             label=f"CNOT count: {count}", 
    #             color=palette[count]
    # )

    #Regular plots
    NUMBER_OF_QUBITS=5
    NLINES = 50
    # NLINES = 45
    palette = sns.color_palette("viridis",NLINES)
    cnot_counts=[1,5,10,15,20,25,30,35,40, 80]
    # cnot_counts=[1,5,10,15,20,25,30,35,40, 80]
    cnot_count=40
    layers_all=[1, 2, 3, 4, 5 , 6]

    # create_all_fidelity_plot_cnots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, NLINES)
    for layers in layers_all:
        create_all_fidelity_plot_cnots_layers(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers)
        create_all_count_percentage_plots_layers(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers)
    create_split_layers_fidelity_plot_cnots(BASE_PATH, NUMBER_OF_QUBITS, cnot_count, ONE_QUBIT_ERROR_SPACE, layers_all)
    # Plots the percentages of circuits we keep.
    # create_all_count_percentage_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, NLINES)

    create_split_layers_count_percentage_plot(BASE_PATH, NUMBER_OF_QUBITS, cnot_count, ONE_QUBIT_ERROR_SPACE, layers_all)

    #SSO stuff
    # create_all_sso_gain_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE)
    # create_raw_sso_plot(BASE_PATH, NUMBER_OF_QUBITS, 30, ONE_QUBIT_ERROR_SPACE)
    # create_baseline_sso_plot(BASE_PATH, NUMBER_OF_QUBITS, 25, ONE_QUBIT_ERROR_SPACE)
    # all_sorted_circ_results=get_all_sorted_circ_results_by_layer(BASE_PATH, cnot_count, NUMBER_OF_QUBITS)
    # print(all_sorted_circ_results)